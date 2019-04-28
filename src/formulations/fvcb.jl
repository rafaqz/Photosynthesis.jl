
# Abstract/mixin implementation of FvCB photosynthesis

"Abstract type for all FvCB photosynthesis"
abstract type AbstractFvCBPhotosynthesis <: AbstractPhotosynthesis end

"Mixin fields for FvCB photosynthesis"
@mix @columns struct MixinFvCBPhoto{V<:AbstractVcJmax,
                                KM<:AbstractCompensation,
                                Ru<:AbstractRubiscoRegen,
                                Re<:Union{Nothing,AbstractRespiration},
                                GS<:AbstractStomatalConductance,
                                SM<:AbstractSoilMethod
                               }
    vcjmax::V         | VcJmax()                       | _ | _ | _ | _
    compensation::KM  | BernacchiCompensation()        | _ | _ | _ | _
    rubisco_regen::Ru | RubiscoRegen()                 | _ | _ | _ | _
    respiration::Re   | Respiration()                  | _ | _ | _ | _       
    gsmodel::GS       | BallBerryStomatalConductance() | _ | _ | _ | _
    soilmethod::SM    | PotentialSoilMethod()          | _ | _ | _ | _
end


check_extremes!(p::AbstractFvCBPhotosynthesis, v) = begin
    if v.jmax <= zero(v.jmax) || v.vcmax <= zero(v.vcmax)
        update_extremes!(p, v)
        return true
    end
    false
end

"""
    photosynthesis!(p::AbstractFvCBPhotosynthesis, v)
Calculates photosynthesis according to the ECOCRAFT
agreed formulation of the Farquharvon Caemmerer (1982) equations.

Farquhar, G.D., S. Caemmerer and J.A. Berry. 1980.
A biochemical model of photosynthetic CO2 assimilation in leaves of C3 species.
Planta. 149:78-90.
"""
function photosynthesis!(p::AbstractFvCBPhotosynthesis, v)
    v.gammastar = co2_compensation_point(p.compensation, v.tleaf) # CO2 compensation point, umol mol-1
    v.km = rubisco_compensation_point(p.compensation, v.tleaf) # Michaelis-Menten for Rubisco, umol mol-1
    v.jmax = max_electron_transport_rate(p.vcjmax, v.tleaf) # Potential electron transport rate, umol m-2 s-1
    v.vcmax = max_rubisco_activity(p.vcjmax, v.tleaf) # Maximum Rubisco activity, umol m-2 s-1
    v.rd = respiration(p.respiration, v.tleaf) # Day leaf respiration, umol m-2 s-1

    soilmoisture_conductance!(p.soilmethod, v)

    v.vj = rubisco_regeneration(p.rubisco_regen, v)

    # Zero values and exit in extreme cases.
    check_extremes!(p, v) && return nothing

    stomatal_conductance!(p, v)
    
    v.ci = v.gs > zero(v.gs) && v.aleaf > zero(v.aleaf) ? v.cs - v.aleaf / v.gs : v.cs
end



# Abstract or concrete implementation of an FvCB energy balance

abstract type AbstractFvCBEnergyBalance <: AbstractEnergyBalance end

@flattenable @default_kw struct FvCBEnergyBalance{
                                    Ra<:AbstractRadiationConductance,
                                    Bo<:AbstractBoundaryConductance,
                                    De<:AbstractDecoupling,
                                    Ev, 
                                    Ph
                                   } <: AbstractFvCBEnergyBalance
    radiation_conductance::Ra | YingPingRadiationConductance()     | true
    boundary_conductance::Bo  | BoundaryConductance()              | true
    decoupling::De            | McNaughtonJarvisDecoupling()       | true
    evapotranspiration::Ev    | PenmanMonteithEvapotranspiration() | true
    photosynthesis::Ph        | BallBerryPhotosynthesis()          | true
    itermax::Int              | 100                                | false
end


enbal_init!(f::AbstractFvCBEnergyBalance, v) = begin
    # Initialise plant with environmental values
    v.tleaf = v.tair
    v.vpdleaf = v.vpd
    v.cs = v.ca

    photo_init!(f.photosynthesis, v)
end 
enbal_update!(f::AbstractFvCBEnergyBalance, v, tleaf) = photo_update!(f.photosynthesis, v, tleaf)


function enbal!(p::AbstractFvCBEnergyBalance, v)
    enbal_init!(p, v)

    # Calculations that don't depend on tleaf
    v.lhv = latent_heat_water_vapour(v.tair)
    v.slope = calc_slope(v.tair)
    v.gradn = radiation_conductance(p.radiation_conductance, v)
    v.gbhu = boundary_conductance_forced(p.boundary_conductance, v)

    converge_tleaf!(p, v) == false && error("leaf temperature convergence failed")

    v.fheat = v.rnet - v.lhv * v.et

    nothing
end

"""
    converge_tleaf!(p, v)
Run energy balance process in a loop to converge on leaf temperature
"""
function converge_tleaf!(p::AbstractFvCBEnergyBalance, v)
    for iter = 1:p.itermax
        photosynthesis!(p.photosynthesis, v)
        vapour_conductance!(p.boundary_conductance, v)

        v.et = penman_monteith(v.pressure, v.slope, v.lhv, v.rnet, v.vpd, v.gh, v.gv)
        v.decoup = calc_decoupling(p.decoupling, v)

        # End of subroutine if no iterations wanted.
        p.itermax == 0 || v.aleaf <= zero(v.aleaf) && return true

        gbc = v.gbh / GBHGBC
        v.cs = v.ca - v.aleaf / gbc
        tleaf = leaftemp(p, v)
        v.vpdleaf = v.et * v.pressure / v.gv

        enbal_update!(p, v, tleaf) # Model-specific var updates

        # Check to see whether convergence was achieved or failed
        if abs(v.tleaf - tleaf) < TOL
            v.tleaf = tleaf
            return true
        end

        v.tleaf = tleaf # Update temperature for another iteration
    end
    false
end

@inline vapour_conductance!(f, v) = begin
    v.gbhf = boundary_conductance_free(f, v)

    # Total boundary layer conductance for heat
    # See Leuning et al (1995) PCE 18:1183-1200 Eqn E5
    v.gbh = v.gbhu + v.gbhf
    # Total conductance for heat - two-sided
    v.gh = 2.0(v.gbh + v.gradn)

    # Total conductance for water vapour
    v.gbv = GBVGBH * v.gbh
    v.gsv = GSVGSC * f.gsc
    # gv = nsides * (gbv * gsv) / (gbv + gsv) # already one-sided value
    v.gv = (v.gbv * v.gsv) / (v.gbv + v.gsv)
end



@MixinEnviroVars @mix struct MixinFvCBVars{μMoMo,kPa,F,WM2,MMoPaJS,μMoM2S,JMo,PaK,MoM2S,MoμMo,μMoMo}
    # shared
    cs::μMoMo        | 400.0       | μmol*mol^-1           | _
    vpdleaf::kPa     | 0.8         | kPa                   | _
    rhleaf::F        | 0.99        | _                     |  "Only in Ball-Berry Stomatal Conductance"
    # energy balance
    fheat::WM2       | 1.0         | W*m^-2                | _
    gbhu::MMoPaJS    | 1.0         | m*mol*Pa*J^-1*s^-1    | _
    gbhf::MMoPaJS    | 1.0         | m*mol*Pa*J^-1*s^-1    | _
    gh::MMoPaJS      | 1.0         | m*mol*Pa*J^-1*s^-1    | _
    gbh::MMoPaJS     | 1.0         | m*mol*Pa*J^-1*s^-1    | _
    gv::MMoPaJS      | 1.0         | m*mol*Pa*J^-1*s^-1    | _
    gradn::MoM2S     | 1.0         | mol*m^-2*s^-1         | _
    lhv::JMo         | 1.0         | J*mol^-1              | _
    et::MoM2S        | 1.0         | mol*m^-2*s^-1         | _
    slope::PaK       | 1.0         | Pa*K^-1               | _
    decoup::F        | 0.0         | _                     | _
    # photosynthesis
    gsdiva::MoμMo    | 1.0         | mol*μmol^-1           | _
    km::μMoMo        | 1.0         | μmol*mol^-1           | _
    ci::μMoMo        | 1.0         | μmol*mol^-1           | _
    gammastar::μMoMo | 1.0         | μmol*mol^-1           | _
    gs::MoM2S        | 1.0         | mol*m^-2*s^-1         | _
    gsv::MoM2S       | 1.0         | mol*m^-2*s^-1         | _
    gbv::MoM2S       | 1.0         | mol*m^-2*s^-1         | _
    jmax::μMoM2S     | 1.0         | μmol*m^-2*s^-1        | _
    vcmax::μMoM2S    | 1.0         | μmol*m^-2*s^-1        | _
    rd::μMoM2S       | 1.0         | μmol*m^-2*s^-1        | _
    ac::μMoM2S       | 1.0         | μmol*m^-2*s^-1        | _
    aleaf::μMoM2S    | 1.0         | μmol*m^-2*s^-1        | _
    vj::μMoM2S       | 1.0         | μmol*m^-2*s^-1        | _
    aj::μMoM2S       | 1.0         | μmol*m^-2*s^-1        | _
    # soil
    fsoil::F         | 1.0         | _                     | _
end

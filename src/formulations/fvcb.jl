
# Abstract/mixin implementation of FvCB photosynthesis

"Abstract type for all FvCB photosynthesis"
abstract type AbstractFvCBPhotosynthesis <: AbstractPhotosynthesis end

"Mixin fields for FvCB photosynthesis"
@columns struct FvCBPhotosynthesis{F,KM,Ru,Re,GS} <: AbstractFvCBPhotosynthesis
    flux::F                  | PotentialModifiedFlux()        | _ | _ | _ | _
    compensation::KM         | BernacchiCompensation()        | _ | _ | _ | _
    rubisco_regen::Ru        | RubiscoRegen()                 | _ | _ | _ | _
    respiration::Re          | Respiration()                  | _ | _ | _ | _       
    stomatal_conductance::GS | BallBerryStomatalConductance() | _ | _ | _ | _
end

check_extremes!(p::AbstractFvCBPhotosynthesis, v) = begin
    if v.jmax <= zero(v.jmax) || v.vcmax <= zero(v.vcmax)
        update_extremes!(p.stomatal_conductance, v)
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
    v.jmax, v.vcmax = flux(p.flux, v) # Potential electron transport rate and maximum Rubisco activity
    v.rd = respiration(p.respiration, v.tleaf) # Day leaf respiration, umol m-2 s-1
    v.vj = rubisco_regeneration(p.rubisco_regen, v)

    # Zero values and exit in extreme cases.
    check_extremes!(p, v) && return nothing

    v.aleaf, v.gs = stomatal_conductance!(p.stomatal_conductance, v)
    
    v.ci = v.gs > zero(v.gs) && v.aleaf > zero(v.aleaf) ? v.cs - v.aleaf / v.gs : v.cs

    v.aleaf
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
    photosynthesis::Ph        | FvCBPhotosynthesis()               | true
    itermax::Int              | 100                                | false
end

enbal_update!(f::AbstractFvCBEnergyBalance, v, tleaf) = nothing

function enbal!(p::AbstractFvCBEnergyBalance, v)

    # Initialise to ambient conditions
    v.tleaf = v.tair
    v.vpdleaf = v.vpd
    v.rhleaf = v.rh
    v.cs = v.ca

    # Calculations that don't depend on tleaf
    v.lhv = latent_heat_water_vapour(v.tair)
    v.slope = calc_slope(v.tair)
    v.gradn = radiation_conductance(p.radiation_conductance, v)
    v.gbhu = boundary_conductance_forced(p.boundary_conductance, v)

    # Converge on leaf temperature
    for iter = 1:p.itermax
        aleaf = photosynthesis!(p.photosynthesis, v)
        v.gv = vapour_conductance!(p.boundary_conductance, v)

        v.et = evapotranspiration(p.evapotranspiration, v)
        v.decoup = calc_decoupling(p.decoupling, v) # only for output?

        # End of subroutine if no iterations wanted.
        p.itermax == 0 || v.aleaf <= zero(v.aleaf) && return true

        gbc = v.gbh / GBHGBC
        v.cs = v.ca - v.aleaf / gbc
        tleaf = leaftemp(p, v)
        v.vpdleaf = v.et * v.pressure / v.gv
        v.rhleaf = 1.0 - v.vpdleaf / saturated_vapour_pressure(tleaf)

        enbal_update!(p, v, tleaf) # Model-specific var updates

        # Check to see whether convergence was achieved or failed
        if abs(v.tleaf - tleaf) < TOL
            v.tleaf = tleaf
            return true
        end

        v.tleaf = tleaf # Update temperature for another iteration

        iter == p.itermax && error("leaf temperature convergence failed")
    end
    # fheat = v.rnet - v.lhv * v.et

    tleaf
end

@MixinEnviroVars @mix struct MixinFvCBVars{μMoMo,kPa,F,WM2,MMoPaJS,μMoM2S,JMo,PaK,MoM2S,MoμMo,μMoMo}
    # shared
    cs::μMoMo        | 400.0       | μmol*mol^-1           | _
    vpdleaf::kPa     | 0.8         | kPa                   | _
    rhleaf::F        | 0.99        | _                     |  "Only in Ball-Berry Stomatal Conductance"
    # energy balance
    fheat::WM2       | 0.0         | W*m^-2                | _
    gbhu::MMoPaJS    | 0.0         | m*mol*Pa*J^-1*s^-1    | _
    gbhf::MMoPaJS    | 0.0         | m*mol*Pa*J^-1*s^-1    | _
    gh::MMoPaJS      | 0.0         | m*mol*Pa*J^-1*s^-1    | _
    gbh::MMoPaJS     | 0.0         | m*mol*Pa*J^-1*s^-1    | _
    gv::MMoPaJS      | 0.0         | m*mol*Pa*J^-1*s^-1    | _
    gradn::MoM2S     | 0.0         | mol*m^-2*s^-1         | _
    lhv::JMo         | 0.0         | J*mol^-1              | _
    et::MoM2S        | 0.0         | mol*m^-2*s^-1         | _
    slope::PaK       | 0.0         | Pa*K^-1               | _
    decoup::F        | 0.0         | _                     | _
    # photosynthesis
    gsdiva::MoμMo    | 0.0         | mol*μmol^-1           | _
    km::μMoMo        | 0.0         | μmol*mol^-1           | _
    ci::μMoMo        | 0.0         | μmol*mol^-1           | _
    gammastar::μMoMo | 0.0         | μmol*mol^-1           | _
    gs::MoM2S        | 0.0         | mol*m^-2*s^-1         | _
    gsv::MoM2S       | 0.0         | mol*m^-2*s^-1         | _
    gbv::MoM2S       | 0.0         | mol*m^-2*s^-1         | _
    jmax::μMoM2S     | 0.0         | μmol*m^-2*s^-1        | _
    vcmax::μMoM2S    | 0.0         | μmol*m^-2*s^-1        | _
    rd::μMoM2S       | 0.0         | μmol*m^-2*s^-1        | _
    ac::μMoM2S       | 0.0         | μmol*m^-2*s^-1        | _
    aleaf::μMoM2S    | 0.0         | μmol*m^-2*s^-1        | _
    vj::μMoM2S       | 0.0         | μmol*m^-2*s^-1        | _
    aj::μMoM2S       | 0.0         | μmol*m^-2*s^-1        | _
    # soil
    fsoil::F         | 0.0         | _                     | _
end

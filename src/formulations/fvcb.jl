
"""
Abstract supertype for all Farquhar/von Caemmerer/Berry derived photosynthesis models
"""
abstract type AbstractFvCBPhotosynthesis <: AbstractPhotosynthesis end

flux_model(p::AbstractFvCBPhotosynthesis) = p.flux_model
compensation_model(p::AbstractFvCBPhotosynthesis) = p.compensation_model
rubisco_regen_model(p::AbstractFvCBPhotosynthesis) = p.rubisco_regen_model
respiration_model(p::AbstractFvCBPhotosynthesis) = p.respiration_model
stomatal_conductance_model(p::AbstractFvCBPhotosynthesis) = p.stomatal_conductance_model


"""
    FvCBPhotosynthesis(flux, compensation, rubisco_regen, respiration, stomatal_conductance)

General Farquhar von Caemmerer Berry model of photosynthesis. Organised as a modular
components of electron flux, CO2 and rubusco compensation, rubisco regeneration, 
respiration and stomatal conductance.

Calculates photosynthesis according to the ECOCRAFT
agreed formulation of the Farquharvon Caemmerer (1982) equations.

Farquhar, G.D., S. Caemmerer and J.A. Berry. 1980.
A biochemical model of photosynthetic CO2 assimilation in leaves of C3 species.
Planta. 149:78-90.
"""
@default_kw struct FvCBPhotosynthesis{F,KM,Ru,Re,GS} <: AbstractFvCBPhotosynthesis
    flux_model::F                  | PotentialModifiedFlux()
    compensation_model::KM         | BernacchiCompensation()
    rubisco_regen_model::Ru        | RubiscoRegen()
    respiration_model::Re          | Respiration()
    stomatal_conductance_model::GS | BallBerryStomatalConductance()
end

check_extremes!(v, p::AbstractFvCBPhotosynthesis) = begin
    if v.jmax <= zero(v.jmax) || v.vcmax <= zero(v.vcmax)
        update_extremes!(stomatal_conductance_model(p), v)
        return true
    end
    false
end

update_extremes!(v, m::AbstractStomatalConductance) = begin
    v.aleaf = -v.rd
    v.gs = g0(m)
end

function photosynthesis!(v, p::AbstractFvCBPhotosynthesis)
    v.gammastar = co2_compensation_point(compensation_model(p), v.tleaf) # CO2 compensation point, umol mol-1
    v.km = rubisco_compensation_point(compensation_model(p), v.tleaf) # Michaelis-Menten for Rubisco, umol mol-1
    v.jmax, v.vcmax = flux(flux_model(p), v) # Potential electron transport rate and maximum Rubisco activity
    v.rd = respiration(respiration_model(p), v.tleaf) # Day leaf respiration, umol m-2 s-1
    v.vj = rubisco_regeneration(rubisco_regen_model(p), v)

    # Zero values and exit in extreme cases.
    check_extremes!(v, p) && return nothing

    v.aleaf, v.gs = stomatal_conductance!(v, stomatal_conductance_model(p))

    v.ci = v.gs > zero(v.gs) && v.aleaf > zero(v.aleaf) ? v.cs - v.aleaf / v.gs : v.cs

    v.aleaf
end


"""
Energy balance models based on genral FvCB photosynthesis, derived from Maespa/Maestra models.
"""
abstract type AbstractFvCBEnergyBalance <: AbstractEnergyBalance end

enbal_update!(v, f::AbstractFvCBEnergyBalance, tleaf) = nothing

"""
    FvCBEnergyBalance(radiation_conductance, 
                      boundary_conductance, 
                      decoupling, 
                      evapotranspiration, 
                      photosynthesis, 
                      max_itererations)

Energy-balance model composed of submodels from radiation conductance, 
boundary conductance, decoupling, evapotranspiration
photosynthesis.

`max_itererations` determines the maximum number of iterations to perform to determine 
flux and temperature.
"""
@default_kw @flattenable struct FvCBEnergyBalance{
                                    Ra<:AbstractRadiationConductance,
                                    Bo<:AbstractBoundaryConductance,
                                    De<:AbstractDecoupling,
                                    Ev,
                                    Ph
                                   } <: AbstractFvCBEnergyBalance
    radiation_conductance_model::Ra | WangRadiationConductance()     | true
    boundary_conductance_model::Bo  | BoundaryConductance()              | true
    decoupling_model::De            | McNaughtonJarvisDecoupling()       | true
    evapotranspiration_model::Ev    | PenmanMonteithEvapotranspiration() | true
    photosynthesis_model::Ph        | FvCBPhotosynthesis()               | true
    max_itererations::Int           | 100                                | false
end

radiation_conductance_model(p::AbstractFvCBEnergyBalance) = p.radiation_conductance_model 
boundary_conductance_model(p::AbstractFvCBEnergyBalance) = p.boundary_conductance_model
decoupling_model(p::AbstractFvCBEnergyBalance) = p.decoupling_model
evapotranspiration_model(p::AbstractFvCBEnergyBalance) = p.evapotranspiration_model
photosynthesis_model(p::AbstractFvCBEnergyBalance) = p.photosynthesis_model
max_itererations(p::AbstractFvCBEnergyBalance) = p.max_itererations

function enbal!(v, p::AbstractFvCBEnergyBalance)

    # Initialise to ambient conditions
    v.tleaf = v.tair
    v.vpdleaf = v.vpd
    v.rhleaf = v.rh
    v.cs = v.ca

    # Calculations that don't depend on tleaf
    v.lhv = latent_heat_water_vapour(v.tair)
    v.slope = slope(v.tair)
    v.gradn = radiation_conductance(radiation_conductance_model(p), v)
    v.gbhu = boundary_conductance_forced(boundary_conductance_model(p), v)

    # Converge on leaf temperature
    for iter = 1:max_itererations(p)
        aleaf = photosynthesis!(v, photosynthesis_model(p))
        v.gv = vapour_conductance!(v, boundary_conductance_model(p))

        v.et = evapotranspiration(evapotranspiration_model(p), v)
        v.decoup = decoupling(decoupling_model(p), v) # only for output?

        # End of subroutine if no iterations wanted.
        max_itererations(p) == 0 || v.aleaf <= zero(v.aleaf) && return true

        gbc = v.gbh / GBHGBC
        v.cs = v.ca - v.aleaf / gbc
        tleaf = leaftemp(p, v)
        v.vpdleaf = v.et * v.pressure / v.gv
        v.rhleaf = 1.0 - v.vpdleaf / saturated_vapour_pressure(tleaf)

        enbal_update!(v, p, tleaf) # Model-specific var updates

        # Check to see whether convergence was achieved or failed
        if abs(v.tleaf - tleaf) < TOL
            v.tleaf = tleaf
            return true
        end

        v.tleaf = tleaf # Update temperature for another iteration

        iter == max_itererations(p) && error("leaf temperature convergence failed")
    end
    # fheat = v.rnet - v.lhv * v.et

    tleaf
end

"""
    @MixinEnviroVars

Mixin variables for [`FvCBEnergyBalance`](@ref) variables objects.
"""
@MixinEnviroVars @mix struct MixinFvCBVars{μMoMo,kPa,F,WM2,MMoPaJS,μMoM2S,JMo,PaK,MoM2S,MoμMo,μMoMo}
    # shared
    cs::μMoMo        | 400.0       | μmol*mol^-1           | _
    vpdleaf::kPa     | 0.8         | kPa                   | _
    rhleaf::F        | 0.99        | _                     | "Only in Ball-Berry Stomatal Conductance"
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
    gs_div_a::MoμMo  | 0.0         | mol*μmol^-1           | _
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

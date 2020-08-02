
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

$(FIELDDOCTABLE)
"""
@default_kw struct FvCBPhotosynthesis{F,KM,Ru,Re,GS} <: AbstractFvCBPhotosynthesis
    flux_model::F                  | PotentialModifiedFlux()
    compensation_model::KM         | BernacchiCompensation()
    rubisco_regen_model::Ru        | RubiscoRegen()
    respiration_model::Re          | Respiration()
    stomatal_conductance_model::GS | BallBerryStomatalConductance()
end

check_extremes!(v, p::AbstractFvCBPhotosynthesis) =
    if v.jmax <= zero(v.jmax) || v.vcmax <= zero(v.vcmax)
        update_extremes!(stomatal_conductance_model(p), v)
        true
    else
        false
    end

update_extremes!(v, m::AbstractStomatalConductance) = begin
    v.aleaf = -v.rd
    v.gs = g0(m)
end

function photosynthesis!(v, m::AbstractFvCBPhotosynthesis)
    v.gammastar = co2_compensation_point(compensation_model(m), v.tleaf) # CO2 compensation point, umol mol-1
    v.km = rubisco_compensation_point(compensation_model(m), v.tleaf) # Michaelis-Menten for Rubisco, umol mol-1
    v.jmax, v.vcmax = flux(flux_model(m), v) # Potential electron transport rate and maximum Rubisco activity
    v.rd = respiration(respiration_model(m), v.tleaf) # Day leaf respiration, umol m-2 s-1
    v.vj = rubisco_regeneration(rubisco_regen_model(m), v)

    # Zero values and exit in extreme cases.
    check_extremes!(v, m) && return

    v.aleaf, v.gs = stomatal_conductance!(v, stomatal_conductance_model(m))

    v.ci = if v.gs > zero(v.gs) && v.aleaf > zero(v.aleaf) 
        v.cs - v.aleaf / v.gs
    else
        v.cs
    end
    return
end


"""
Energy balance models based on genral FvCB photosynthesis, derived from Maespa/Maestra models.
"""
abstract type AbstractFvCBEnergyBalance <: AbstractEnergyBalance end

"""
    FvCBEnergyBalance(radiation_conductance, 
                      boundary_conductance, 
                      decoupling, 
                      evapotranspiration, 
                      photosynthesis, 
                      max_iterations)

Energy-balance model composed of submodels from radiation conductance, 
boundary conductance, decoupling, evapotranspiration
photosynthesis.

`max_iterations` determines the maximum number of iterations to perform to determine 
flux and temperature.

$(FIELDDOCTABLE)
"""
@default_kw @flattenable struct FvCBEnergyBalance{
        Ra<:AbstractRadiationConductance,
        Bo<:AbstractBoundaryConductance,
        De<:AbstractDecoupling,Ev,Ph,I,A} <: AbstractFvCBEnergyBalance
    radiation_conductance_model::Ra | WangRadiationConductance()         | true
    boundary_conductance_model::Bo  | BoundaryConductance()              | true
    decoupling_model::De            | McNaughtonJarvisDecoupling()       | true
    evapotranspiration_model::Ev    | PenmanMonteithEvapotranspiration() | true
    photosynthesis_model::Ph        | FvCBPhotosynthesis()               | true
    max_iterations::I               | 100                                | false
    atol::A                         | TOL                                | false
end

radiation_conductance_model(p::AbstractFvCBEnergyBalance) = p.radiation_conductance_model 
boundary_conductance_model(p::AbstractFvCBEnergyBalance) = p.boundary_conductance_model
decoupling_model(p::AbstractFvCBEnergyBalance) = p.decoupling_model
evapotranspiration_model(p::AbstractFvCBEnergyBalance) = p.evapotranspiration_model
photosynthesis_model(p::AbstractFvCBEnergyBalance) = p.photosynthesis_model
max_iterations(p::AbstractFvCBEnergyBalance) = p.max_iterations
atol(p::AbstractFvCBEnergyBalance) = p.atol

"""
    enbal!(v, m::AbstractFvCBEnergyBalance)

Calculates leaf photosynthesis and transpiration for an `AbstractEnergyBalance`
model `m` and variables `v`.

Results are written to v.

These may be calculated by:

1. Assuming leaf temperature is the same as air temperature, 
   and stomatal carbon has the same conentration as in the air. 
2. Using iterative scheme of Leuning et al (1995) (PCE 18:1183-1200) 
   to calculate leaf temperature and stomatal carbon concentration.

Setting `max_iterations=0` gives 1, max_iterations > 0 (default 100) gives 2.
"""
function enbal!(v, m::AbstractFvCBEnergyBalance)

    # Initialise to ambient conditions
    v.tleaf = v.tair
    v.vpdleaf = v.vpd
    v.rhleaf = v.rh
    v.cs = v.ca

    # Calculations that don't depend on tleaf
    v.lhv = latent_heat_water_vapour(v.tair)
    v.slope = slope(v.tair)
    v.gradn = radiation_conductance(radiation_conductance_model(m), v)
    v.gbhu = boundary_conductance_forced(boundary_conductance_model(m), v)

    iter = 1
    # Converge on leaf temperature
    while true 
        photosynthesis!(v, photosynthesis_model(m))
        conductance!(v, m)
        # This isn't actually used - it's only for output
        v.decoup = decoupling(decoupling_model(m), v)

        # End of subroutine if no iterations wanted.
        (max_iterations(m) == 0 || v.aleaf <= zero(v.aleaf)) && return true

        gbc = v.gbh / GBHGBC
        v.cs = v.ca - v.aleaf / gbc
        tleaf = leaftemp(m, v)

        # Recalculate
        conductance!(v, m)

        v.vpdleaf = v.et * v.pressure / v.gv
        v.rhleaf = 1 - v.vpdleaf / saturated_vapour_pressure(tleaf)

        # Check to see whether convergence has occurred
        if abs(v.tleaf - tleaf) < atol(m)
            v.tleaf = tleaf
            return true
        end

        v.tleaf = tleaf # Update temperature for another iteration

        iter >= max_iterations(m) && break
        iter += 1
    end

    @warn "leaf temperature convergence failed"
    return false
end

function conductance!(v, m)
    # Total boundary layer conductance for heat
    # See Leuning et al (1995) PCE 18:1183-1200 Eqn E5
    v.gbhf = boundary_conductance_free(boundary_conductance_model(m), v)
    v.gbh = v.gbhu + v.gbhf
    # Total conductance for heat: two-sided
    v.gh = 2.0(v.gbh + v.gradn)
    # Total conductance for water vapour
    gbv = GBVGBH * v.gbh
    gsv = GSVGSC * v.gs
    v.gv = (gbv * gsv) / (gbv + gsv)
    v.et = evapotranspiration(evapotranspiration_model(m), v)
end

"""
    @MixinEnviroVars

Mixin variables for [`FvCBEnergyBalance`](@ref) variables objects.
"""
@MixinEnviroVars @mix struct MixinFvCBVars{μMoMo,kPa,F,WM2,μMoM2S,JMo,PaK,MoM2S,MoμMo,μMoMo}
    # shared
    cs::μMoMo        | 400.0  | μmol*mol^-1    | _
    vpdleaf::kPa     | 0.8    | kPa            | _
    rhleaf::F        | 0.99   | _              | "Only in Ball-Berry Stomatal Conductance"
    # energy balance
    fheat::WM2       | 0.0    | W*m^-2         | _
    gbhu::MoM2S      | 0.0    | mol*m^-2*s^-1  | _
    gbhf::MoM2S      | 0.0    | mol*m^-2*s^-1  | _
    gh::MoM2S        | 0.0    | mol*m^-2*s^-1  | _
    gbh::MoM2S       | 0.0    | mol*m^-2*s^-1  | _
    gv::MoM2S        | 0.0    | mol*m^-2*s^-1  | _
    gradn::MoM2S     | 0.0    | mol*m^-2*s^-1  | _
    lhv::JMo         | 0.0    | J*mol^-1       | _
    et::MoM2S        | 0.0    | mol*m^-2*s^-1  | _
    slope::PaK       | 0.0    | Pa*K^-1        | _
    decoup::F        | 0.0    | _              | _
    # photosynthesis
    gs_div_a::MoμMo  | 0.0    | mol*μmol^-1    | _
    km::μMoMo        | 0.0    | μmol*mol^-1    | _
    ci::μMoMo        | 0.0    | μmol*mol^-1    | _
    gammastar::μMoMo | 0.0    | μmol*mol^-1    | _
    gs::MoM2S        | 0.0    | mol*m^-2*s^-1  | _
    gsv::MoM2S       | 0.0    | mol*m^-2*s^-1  | _
    gbv::MoM2S       | 0.0    | mol*m^-2*s^-1  | _
    jmax::μMoM2S     | 0.0    | μmol*m^-2*s^-1 | _
    vcmax::μMoM2S    | 0.0    | μmol*m^-2*s^-1 | _
    rd::μMoM2S       | 0.0    | μmol*m^-2*s^-1 | _
    ac::μMoM2S       | 0.0    | μmol*m^-2*s^-1 | _
    aleaf::μMoM2S    | 0.0    | μmol*m^-2*s^-1 | _
    vj::μMoM2S       | 0.0    | μmol*m^-2*s^-1 | _
    aj::μMoM2S       | 0.0    | μmol*m^-2*s^-1 | _
    # soil
    fsoil::F         | 1.0    | _              | _
end

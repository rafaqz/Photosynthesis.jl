
"""
Energy balance models derived from MAESPA/MAESTRA models.
"""
abstract type AbstractMaespaEnergyBalance <: AbstractEnergyBalance end


"""
    MaespaEnergyBalance(radiation_conductance, 
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
@flattenable @default_kw @description struct MaespaEnergyBalance{
        Ra<:AbstractRadiationConductance,
        Bo<:AbstractBoundaryConductance,
        De<:AbstractDecoupling,Ev,Ph,I,A} <: AbstractMaespaEnergyBalance
    radiation_conductance_model::Ra  | true  | WangRadiationConductance()         | "Radiation condutance model"
    boundary_conductance_model::Bo   | true  | BoundaryConductance()              | "Boundary conductance model"
    decoupling_model::De             | true  | McNaughtonJarvisDecoupling()       | ""
    evapotranspiration_model::Ev     | true  | PenmanMonteithEvapotranspiration() | "Evapotranspiration model"
    photosynthesis_model::Ph         | true  | FvCBPhotosynthesis()               | "Photosynthesis model"
    max_iterations::I                | false | 100                                | "Number of iterations used to find leaf temperature"
    atol::A                          | false | 0.02K                              | "Tolerance for difference with previous leaf temperature"
end

radiation_conductance_model(p::AbstractMaespaEnergyBalance) = p.radiation_conductance_model 
boundary_conductance_model(p::AbstractMaespaEnergyBalance) = p.boundary_conductance_model
decoupling_model(p::AbstractMaespaEnergyBalance) = p.decoupling_model
evapotranspiration_model(p::AbstractMaespaEnergyBalance) = p.evapotranspiration_model
photosynthesis_model(p::AbstractMaespaEnergyBalance) = p.photosynthesis_model
max_iterations(p::AbstractMaespaEnergyBalance) = p.max_iterations
atol(p::AbstractMaespaEnergyBalance) = p.atol

"""
    enbal!(v, m::AbstractMaespaEnergyBalance)

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
function enbal!(v, m::AbstractMaespaEnergyBalance)

    # Initialise to ambient conditions
    v.tleaf = v.tair
    v.vpdleaf = v.vpd
    v.rhleaf = v.rh
    v.cs = v.ca

    # Calculations that don't depend on tleaf
    v.lhv = latent_heat_water_vapour(v.tair)
    v.gradn = radiation_conductance(radiation_conductance_model(m), v)
    v.gbhu = boundary_conductance_forced(boundary_conductance_model(m), v)
    # Slope is precalculated as it's expensive
    v.slope = vapour_pressure_slope(evapotranspiration_model(m), v)

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
    v.gbv = GBVGBH * v.gbh
    v.gsv = GSVGSC * v.gs
    v.gv = (v.gbv * v.gsv) / (v.gbv + v.gsv)
    v.et = evapotranspiration(evapotranspiration_model(m), v)
end


"""
    @MixinMaespaVars

Mixin variables for [`AbstractMaespaEnergyBalance`](@ref) variables objects.
"""
@mix @vars struct MixinMaespaVars{TL,VPDL,RHL,CD,FHi,GBHU,GBHF,GH,GSV,GBV,GBH,GV,GR,LHV,ET,LS,DC}
    # shared
    tleaf::TL        | 298.15 | K              | _
    vpdleaf::VPDL    | 0.0    | Pa             | _
    rhleaf::RHL      | 0.0    | _              | _
    cs::CD           | 0.0    | μmol*mol^-1    | _
    # energy balance
    fheat::FHi       | 0.0    | W*m^-2         | _
    gbhu::GBHU       | 0.0    | mol*m^-2*s^-1  | _
    gbhf::GBHF       | 0.0    | mol*m^-2*s^-1  | _
    gh::GH           | 0.0    | mol*m^-2*s^-1  | _
    gsv::GSV         | 0.0    | mol*m^-2*s^-1  | _
    gbv::GBV         | 0.0    | mol*m^-2*s^-1  | _
    gbh::GBH         | 0.0    | mol*m^-2*s^-1  | _
    gv::GV           | 0.0    | mol*m^-2*s^-1  | _
    gradn::GR        | 0.0    | mol*m^-2*s^-1  | _
    lhv::LHV         | 0.0    | J*mol^-1       | _
    et::ET           | 0.0    | mol*m^-2*s^-1  | _
    slope::LS        | 0.0    | Pa*K^-1        | _
    decoup::DC       | 0.0    | _              | _
end

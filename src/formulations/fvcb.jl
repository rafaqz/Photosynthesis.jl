abstract type AbstractFvCBEnergyBalance <: AbstractEnergyBalance end


@flattenable @default_kw struct FvCBEnergyBalance{Ph,
                                    Ra<:AbstractRadiationConductance,
                                    Bo<:AbstractBoundaryConductance,
                                    De<:AbstractDecoupling} <: AbstractFvCBEnergyBalance
    radiation_conductance::Ra | YingPingRadiationConductance() | true
    boundary_conductance::Bo  | BoundaryConductance()          | true
    decoupling::De            | McNaughtonJarvisDecoupling()   | true
    photo::Ph                 | FvCBPhoto()                    | true
    itermax::Int              | 100                            | false
end

run_enbal!(v, p::AbstractFvCBEnergyBalance) = run_enbal!(p, v)

"""
    enbal!(v, p)
This subroutine calculates leaf photosynthesis and transpiration.

These may be calculated by:
(1) assuming leaf temperature = air temperature, cs = ca and ds = da
(2) using iterative scheme of Leuning et al (1995) (PCE 18:1183-1200) to calculate leaf temp, CsCa.

Setting itermax = 0 gives (1); itermax > 0 (suggest 100) gives (2).
"""
function enbal!(p::AbstractFvCBEnergyBalance, v)
    enbal_init!(p, v)

    # Initialise with environmental values
    v.tleaf = v.tair
    v.vpdleaf = v.vpd
    v.cs = v.ca

    # Calculations that do not depend on tleaf
    v.lhv = latent_heat_water_vapour(v.tair)
    v.slope = calc_slope(v.tair)
    v.gradn = radiation_conductance(p.radiation_conductance, v, p)
    v.gbhu = boundary_conductance_forced(p.boundary_conductance, v, p)

    converge_tleaf!(f, v, p) == false && error("leaf temperature convergence failed")

    v.fheat = v.rnet - v.lhv * v.et

    nothing
end

"""
    converge_tleaf!(v, p)
Run energy balance process in a loop to converge on leaf temperature
"""
function converge_tleaf!(p::AbstractFvCBEnergyBalance, v)
    for iter = 1:p.itermax
        photosynthesis!(p.photo, v)
        v.gbhf = boundary_conductance_free(p.boundary_conductance, v, p)

        # Total boundary layer conductance for heat
        # See Leuning et al (1995) PCE 18:1183-1200 Eqn E5
        gbh = v.gbhu + v.gbhf
        # Total conductance for heat - two-sided
        v.gh = 2.0(gbh + v.gradn)

        # Total conductance for water vapour
        gbv = GBVGBH * gbh
        gsv = GSVGSC * v.gs # TODO check this is right?? was: * p.gsc
        # gv = nsides * (gbv * gsv) / (gbv + gsv) # already one-sided value
        gv = (gbv * gsv) / (gbv + gsv)

        v.et = penman_monteith(v.pressure, v.slope, v.lhv, v.rnet, v.vpd, v.gh, gv)
        v.decoup = calc_decoupling(p.decoupling, v, p, gbv, gsv)

        # End of subroutine if no iterations wanted.
        p.itermax == 0 || v.aleaf <= zero(v.aleaf) && return true

        gbc = gbh / GBHGBC
        v.cs = v.ca - v.aleaf / gbc # TODO this value is way too low
        tleaf1 = leaftemp(v, p)
        v.vpdleaf = v.et * v.pressure / gv # TODO and this seems too high?
        uconvert(kPa, v.vpdleaf)

        model_update!(v, f, p, tleaf1) # Model-specific var updates

        # Check to see whether convergence achieved or failed
        if abs(v.tleaf - tleaf1) < TOL
            v.tleaf = tleaf1
            return true
        end

        v.tleaf = tleaf1 # Update temperature for another iteration
    end
    false
end

leaftemp(v, p) = v.tair + (v.rnet - v.et * v.lhv) / (CPAIR * AIRMA * v.gh)

enbal_init!(f::AbstractFvCBEnergyBalance, v, p) = photo_init!(f.photo, v, p) 
enbal_update!(f::AbstractFvCBEnergyBalance, v, p, tleaf) = photo_update!(f.photo, v, p, tleaf)

photo_init!(f::AbstractFvCBPhoto, v, p) = photo_init!(f.photo.model, v, p) 
photo_update!(f::AbstractFvCBPhoto, v, p, tleaf) = photo_update!(f.photo.model, v, p, tleaf)

"""
Fixed parameters for photosynthesis.

THe majority of these are "composed" from submodels that both hold
the parameters and use their own specific methods during the photosynthesis
routines.

Calling `PhotoParams()` will give the default values for all of these submodels.
Any parameters and submodels can be overridden with keyword arguments:

`PhotoParams(model=TuzetModel, ca= 450 | μmol*mol^-1)`
"""
@default_kw struct FvCBPhoto{F<:AbstractPhotoModel,
                             V<:AbstractVcJmax,
                             KM<:AbstractCompensation,
                             Ru<:AbstractRubiscoRegen,
                             Re<:Union{Nothing,AbstractRespiration}}
    model::F          | BallBerryModel()       
    vcjmax::V         | VcJmax()               
    compensation::KM  | BernacchiCompensation()
    rubisco_regen::Ru | RubiscoRegen()         
    respiration::Re   | nothing          
end

"""
    photosynthesis!(v, p)
Calculates photosynthesis according to the ECOCRAFT
agreed formulation of the Farquharvon Caemmerer (1982) equations.

Farquhar, G.D., S. Caemmerer and J.A. Berry. 1980. 
A biochemical model of photosynthetic CO2 assimilation in leaves of C3 species. 
Planta. 149:78-90. 
"""
function photosynthesis!(p::AbstractFvCBPhoto, v)
    v.gammastar = co2_compensation_point(p.compensation, v, p) # CO2 compensation point, umol mol-1
    v.km = rubisco_compensation_point(p.compensation, v, p) # Michaelis-Menten for Rubisco, umol mol-1
    v.jmax = max_electron_transport_rate(p.vcjmax, v, p) # Potential electron transport rate, umol m-2 s-1
    v.vcmax = max_rubisco_activity(p.vcjmax, v, p) # Maximum Rubisco activity, umol m-2 s-1

    v.rd = respiration(p.respiration, v, p, v.rd) # Day leaf respiration, umol m-2 s-1

    soilmoisture_conductance!(p.model.soilmethod, v, p)

    v.vj = rubisco_regeneration(p.rubisco_regen, v, p)

    extremes!(v, p) && return nothing
    stomatal_conductance!(p.model, v, p)

    nothing
end

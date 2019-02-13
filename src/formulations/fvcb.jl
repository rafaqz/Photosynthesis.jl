abstract type AbstractFvCBEnergyBalance <: AbstractEnergyBalance end


@flattenable @default_kw struct FvCBEnergyBalance{Ph,
                                    Ra<:AbstractRadiationConductance,
                                    Bo<:AbstractBoundaryConductance,
                                    De<:AbstractDecoupling} <: AbstractFvCBEnergyBalance
    radiation_conductance::Ra | YingPingRadiationConductance() | true
    boundary_conductance::Bo  | BoundaryConductance()          | true
    decoupling::De            | McNaughtonJarvisDecoupling()   | true
    photo::Ph                 | BallBerryModel()               | true
    itermax::Int              | 100                            | false
end

run_enbal!(p::AbstractFvCBEnergyBalance, v) = enbal!(p, v)

"""
    enbal!(p, v)
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
        photosynthesis!(p.photo, v)
        vapour_conductance!(p, v)

        v.et = penman_monteith(v.pressure, v.slope, v.lhv, v.rnet, v.vpd, v.gh, gv)
        v.decoup = calc_decoupling(p.decoupling, v, gbv, gsv)

        # End of subroutine if no iterations wanted.
        p.itermax == 0 || v.aleaf <= zero(v.aleaf) && return true

        v.gbc = gbh / GBHGBC
        v.cs = v.ca - v.aleaf / v.gbc # TODO this value is way too low
        tleaf = leaftemp(p, v)
        v.vpdleaf = v.et * v.pressure / v.gv # TODO and this seems too high?

        enbal_update!(p, v, tleaf) # Model-specific var updates

        # Check to see whether convergence achieved or failed
        if abs(v.tleaf - tleaf) < TOL
            v.tleaf = tleaf
            return true
        end

        v.tleaf = tleaf1 # Update temperature for another iteration
    end
    false
end

@inline vapour_conductance(p, v) = begin
    v.gbhf = boundary_conductance_free(p, v)

    # Total boundary layer conductance for heat
    # See Leuning et al (1995) PCE 18:1183-1200 Eqn E5
    v.gbh = v.gbhu + v.gbhf
    # Total conductance for heat - two-sided
    v.gh = 2.0(gbh + v.gradn)

    # Total conductance for water vapour
    gbv = GBVGBH * v.gbh
    gsv = GSVGSC * p.gsc
    # gv = nsides * (gbv * gsv) / (gbv + gsv) # already one-sided value
    v.gv = (gbv * gsv) / (gbv + gsv)
end

leaftemp(p, v) = v.tair + (v.rnet - v.et * v.lhv) / (CPAIR * AIRMA * v.gh)

enbal_init!(f::AbstractFvCBEnergyBalance, v) = photo_init!(f.photo, v)
enbal_update!(f::AbstractFvCBEnergyBalance, v, tleaf) = photo_update!(f.photo, v, tleaf)

abstract type AbstractFvCBPhotosynthesis <: AbstractPhotosynthesis end

check_extremes!(p, v) = begin
    if v.jmax <= zero(v.jmax) || v.vcmax <= zero(v.vcmax)
        println("extreme")
        update_extremes!(p, v)
        return true
    end
    false
end

"""
    photosynthesis!(p, v)
Calculates photosynthesis according to the ECOCRAFT
agreed formulation of the Farquharvon Caemmerer (1982) equations.

Farquhar, G.D., S. Caemmerer and J.A. Berry. 1980.
A biochemical model of photosynthetic CO2 assimilation in leaves of C3 species.
Planta. 149:78-90.
"""
function photosynthesis!(p, v)
    v.gammastar = co2_compensation_point(p.compensation, v) # CO2 compensation point, umol mol-1
    v.km = rubisco_compensation_point(p.compensation, v) # Michaelis-Menten for Rubisco, umol mol-1
    v.jmax = max_electron_transport_rate(p.vcjmax, v) # Potential electron transport rate, umol m-2 s-1
    v.vcmax = max_rubisco_activity(p.vcjmax, v) # Maximum Rubisco activity, umol m-2 s-1

    v.rd = respiration(p.respiration, v) # Day leaf respiration, umol m-2 s-1

    soilmoisture_conductance!(p.soilmethod, v)

    v.vj = rubisco_regeneration(p.rubisco_regen, v)

    check_extremes!(p, v) && return nothing
    stomatal_conductance!(p, v)

    nothing
end

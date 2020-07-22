"""
Evapotranspiration models, define an [`evapotranspiration`](@ref) method.
"""
abstract type AbstractEvapotranspiration end

"""
    evapotranspiration(f::AbstractEvapotranspiration, v)

Calculate the rate of evapotranspiration in u"mol*m^-2*s^-1".

Where v is a variables object.
"""
function evapotranspiration end

"""
    PenmanMonteithEvapotranspiration()

Penman-Monteith evapotranspiration model
"""
struct PenmanMonteithEvapotranspiration <: AbstractEvapotranspiration end

evapotranspiration(f::PenmanMonteithEvapotranspiration, v) =
    penman_monteith(v.pressure, v.slope, v.lhv, v.rnet, v.vpd, v.gh, v.gv)

"""
    penman_monteith(pressure, slope, lhv, rnet, vpd, gh, gv)

This subroutine calculates evapotranspiration by leaves using the Penman - Monteith equation.

Inputs:

ρ atmospheric pressure, Pa
Δ slope of VPD / T curve, Pa K-1
lhv latent heat of water at air T, J mol-1
Rn net radiation, J m-2 s-1
Da vapour pressure deficit of air, Pa
gh boundary layer conductance to heat (freeforcedradiative components), mol m-2 s-1
gv conductance to water vapour (stomatalbdry layer components), mol*m^-2*s^-1

Result in mol*m^-2*s^-1
"""
penman_monteith(ρa, Δ, lhv, Rn, Da, gh, gv) = begin
    #FIXME gv <= zero(gv) && return zero(Rn / lhv)

    γ = CPAIR * ρa * AIRMA / lhv
    ET = (Δ * Rn + CPAIR * gh * Da * AIRMA) / (Δ + γ * gh * (1/gv))

    return ET / lhv
    # if (penmon < 0.0) penmon = 0.0 end # BM 12 / 05 Should not be negative
end


""" 
    slope(tair)

Calculate vapour pressure change with temperature -
Slope `s` for Penman-Monteith equation, in Pa K^-1

this is somewhat expensive with two calls to exp(), so it is good to precalculate instead of repeating
"""
slope(tair) = (saturated_vapour_pressure(tair + 0.1oneunit(tair)) -
               saturated_vapour_pressure(tair)) / 0.1oneunit(tair)

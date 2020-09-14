"""
Evapotranspiration models, define an [`evapotranspiration`](@ref) method.
"""
abstract type AbstractEvapotranspiration end

"""
    evapotranspiration(m::AbstractEvapotranspiration, v)

Calculate the rate of leaf evapotranspiration in `u"mol*m^-2*s^-1"` 
for some evapotranspiration model `m`, where `v` is a variables object.
"""
function evapotranspiration end

"""
    PenmanMonteithEvapotranspiration()

Calculates leaf evapotranspiration using the Penman-Monteith equation.
"""
struct PenmanMonteithEvapotranspiration <: AbstractEvapotranspiration end

evapotranspiration(f::PenmanMonteithEvapotranspiration, v) =
    penman_monteith_evapotranspiration(v.pressure, v.slope, v.lhv, v.rnet, v.vpd, v.gh, v.gv)

"""
    penman_monteith_evapotranspiration(pressure, slope, lhv, rnet, vpd, gh, gv)

Calculates leaf evapotranspiration using the Penman-Monteith equation.

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
penman_monteith_evapotranspiration(ρa, Δ, lhv, Rn, Da, gh, gv) =
    if gv > zero(gv) 
        γ = CPAIR * ρa * AIRMA / lhv
        (Δ * Rn + CPAIR * gh * Da * AIRMA) / (Δ + γ * gh / gv) / lhv
    else
        zero(Rn / lhv)
    end


""" 
    slope(tair)

Calculate vapour pressure change with temperature,
`slope` or `Δ` for the Penman-Monteith equation. 

Returns value in `u"Pa*K^-1"`.

Expensive to calculate so should be precomputed separately to 
evapotranspiration where possible.
"""
vapour_pressure_slope(m::PenmanMonteithEvapotranspiration, v) = 
    (saturated_vapour_pressure(v.tair + 0.1oneunit(v.tair)) -
     saturated_vapour_pressure(v.tair)) / 0.1oneunit(v.tair)

"""
    saturated_vapour_pressure(tair)

Calculate saturated water vapour pressure in `u"Pa"` 
at air temperature `tair` in `u"K"`.

From Jones 1992 p 110 (note error in a - wrong units)

TODO: name the magic numbers
      explain the °C multiplication: this is an empirical hack
"""
@inline saturated_vapour_pressure(tair) = begin
    T = ustrip(°C, tair)
    613.75Pa * exp(17.502 * T / (240.97 + T))
end

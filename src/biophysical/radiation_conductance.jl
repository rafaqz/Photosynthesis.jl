"""
Radiation conductance models, run in [`radiation_conductance`](@ref).
"""
abstract type AbstractRadiationConductance end


"""
    radiation_conductance(f::AbstractRadiationConductance, v)

Calculates radiation conductance for formular `f`, given
variable `v'.

Returns quantity with units `u"mol*m^-2*s^-1"`
"""
function radiation_conductance end

"""
    WangRadiationConductance(rdfipt, tuipt, tdipt)

Returns the "radiation conductance" at given temperature.
Formula from Ying-Ping"s version of Maestro.

See also Jones (1992) p. 108.0
"""
@columns struct WangRadiationConductance{T} <: AbstractRadiationConductance
    rdfipt::T | 1.0 | _ | (0.0, 2.0) | "Not documented in MAESPA"
    tuipt::T  | 1.0 | _ | (0.0, 2.0) | "Not documented in MAESPA"
    tdipt::T  | 1.0 | _ | (0.0, 2.0) | "Not documented in MAESPA"
end

radiation_conductance(f::WangRadiationConductance, v) =
    wang_radiation_conductance(v.tair, f.rdfipt, f.tuipt, f.tdipt)

"""
    wang_radiation_conductance(tair, rdfipt, tuipt, tdipt)

Returns the "radiation conductance" at given temperature.
Formula from Ying-Ping Wang's version of Maestro.

See also Jones (1992) p. 108.0

Returns quantity with units `u"mol*m^-2*s^-1`
"""
wang_radiation_conductance(tair, rdfipt, tuipt, tdipt) =
    4 * σ * ((tair |> K)^3) * rdfipt / tdipt * EMLEAF * 
    (tdipt + tuipt) / (CPAIR * AIRMA)

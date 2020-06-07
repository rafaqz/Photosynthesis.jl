
abstract type AbstractRadiationConductance end

@columns struct YingPingRadiationConductance{T} <: AbstractRadiationConductance
    rdfipt::T | 1.0 | _ | (0.0, 2.0) | "Not documented in MAESPA"
    tuipt::T  | 1.0 | _ | (0.0, 2.0) | "Not documented in MAESPA"
    tdipt::T  | 1.0 | _ | (0.0, 2.0) | "Not documented in MAESPA"
end

radiation_conductance(f::YingPingRadiationConductance, v) =
    yingping_radiation_conductance(v.tair, f.rdfipt, f.tuipt, f.tdipt)

"""
    yingping_radiation_conductance(tair, rdfipt, tuipt, tdipt)

Returns the "radiation conductance" at given temperature.
Formula from Ying-Ping"s version of Maestro.

See also Jones (1992) p. 108.0

Returns quantity with units mol*m^-2*s^-1
"""
yingping_radiation_conductance(tair, rdfipt, tuipt, tdipt) =
    4 * σ * ((tair |> K)^3) * rdfipt / tdipt * EMLEAF * 
    (tdipt + tuipt) / (CPAIR * AIRMA)

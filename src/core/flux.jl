"""
Electron flux formulation, calculated in [`max_electron_transport_rate`](@ref).
"""
abstract type AbstractJmax end

"""
    Jmax(jmax25, delsj, eavj, edvj)

Standard formulation for maximum electron transport rate.

$(FIELDDOCTABLE)
"""
@columns struct Jmax{μMoM2S,JMoK,JMo} <: AbstractJmax
    jmax25::μMoM2S  | 184.0    | μmol*m^-2*s^-1 | (0.0, 1000.0)    | "Maximum rate of electron transport at 25° C"
    delsj::JMoK     | 640.02   | J*mol^-1*K^-1  | (0.0, 1000.0)    | "DELTAS in Medlyn et al. (2002)"
    eavj::JMo       | 37259.0  | J*mol^-1       | (0.0, 100000.0)  | "Ha in Medlyn et al. (2002)"
    edvj::JMo       | 200000.0 | J*mol^-1       | (0.0, 1000000.0) | "Hd in Medlyn et al. (2002)"
end

"""
    max_electron_transport_rate(f::Jmax, tleaf)

Calculates the potential max_electron transport rate (Jmax) at the leaf temperature
"""
max_electron_transport_rate(f::Jmax, tleaf) = begin
    tleafK = tleaf |> K
    K25 = 25.0°C |> K
    f.jmax25 * exp((tleafK - K25) * f.eavj / (R * tleafK * K25)) *
    (1 + exp((f.delsj * K25 - f.edvj) / (R * K25))) /
    (1 + exp((f.delsj * tleafK - f.edvj) / (R * tleafK)))
end


"""
Maximum rubisco activity formulations, 
calculated in [`max_rubisco_activity`](@ref).
"""
abstract type AbstractVcmax end

" Mixin fields for maximum rubisco transport rate "
@mix @columns struct Vcmax{μMoM2S,JMo}
    vcmax25::μMoM2S | 110.0    | μmol*m^-2*s^-1 | (0.0, 200.0)    | "Maximumrate rate of rubisco activity at 25° C"
    eavc::JMo       | 47590.0  | J*mol^-1       | (0.0, 100000.0) | "Ha Medlyn et al. (2002)"
end

"""
    max_rubisco_activity(formulation::AbstractVcmax, tleaf)

Calculates the maximum Rubisco activity (Vcmax) at the leaf temperature `tleaf`.

There is still disagreement as to whether this function has an optimum or not.
Both versions here are well-behaved for `tleaf < 0.0`
"""
function max_rubisco_activity end

"""
    NoOptimumVcmax(vcmax25, eavc)

Formulation for maximum Rubisco activity with no optimum.

$(FIELDDOCTABLE)
"""
@Vcmax struct NoOptimumVcmax{} <: AbstractVcmax end

max_rubisco_activity(f::NoOptimumVcmax, tleaf) = begin
    tleafK = K(tleaf)
    K25 = K(25°C)
    f.vcmax25 * exp((f.eavc * (tleaf - K25)) / (K25 * R * tleafK))
end

"""
    OptimumVcmax(edvc, delsc, vcmax25, eavc)

Formulation for maximum Rubisco activity with an optimum.

$(FIELDDOCTABLE)
"""
@Vcmax struct OptimumVcmax{JMo,JMoK} <: AbstractVcmax
    edvc::JMo       | 1.0      | J*mol^-1       | (0.0, 10.0)    | "Hd in Medlyn et al. (2002)"
    delsc::JMoK     | 629.26   | J*mol^-1*K^-1  | (0.0, 2000.0)  | "DELTAS in Medlyn et al. (2002)"
end

max_rubisco_activity(f::OptimumVcmax, tleaf) = begin
    tleafK = K(tleaf)
    K25 = K(25°C)
    f.vcmax25 * exp((tleaf - K25) * f.eavc / (R * tleafK * K25)) *
    (1.0 + exp((f.delsc * K25 - f.edvc) / (R * K25))) /
    (1.0 + exp((f.delsc * tleafK - f.edvc) / (R * tleafK)))
end


"""
Abstract supertype for flux models.

Electron flux Jmax and Rubisco activity Vcmax are often modified 
by the same function, so we group them.
"""
abstract type AbstractFlux end

"""
    flux(f::AbstractFlux, v)

Run jamax and vcmax formulations and any modifications, returning a 2-tuple
"""
function flux end

"""
    Flux(jmaxformulation, vcmaxformulation)

Formulation grouping jmax and vcmax formultions

$(FIELDDOCTABLE)
"""
@udefault_kw struct Flux{J,V} <: AbstractFlux
    jmaxformulation::J  | Jmax()
    vcmaxformulation::V | NoOptimumVcmax()
end

flux_model(x::Flux) = x
flux_model(x) = flux_model(x)

flux(f::Flux, v) =
    max_electron_transport_rate(f.jmaxformulation, v.tleaf),
    max_rubisco_activity(f.vcmaxformulation, v.tleaf)

"""
    DukeFlux(flux_model, tvjup, tvjdn)

Flux model modified that allow Jmax and Vcmax to be forced linearly to zero at low T.

$(FIELDDOCTABLE)
"""
@columns struct DukeFlux{F,K} <: AbstractFlux
    flux_model::F       | Flux() | _  | _              | _
    tvjup::K            | 283.15 | K  | (250.0, 350.0) | _
    tvjdn::K            | 273.15 | K  | (250.0, 350.0) | _
end

flux_model(m::DukeFlux) = m.flux_model

flux(f::DukeFlux, v) = begin
    jmax, vcmax = flux(flux_model(f), v)
    v.tleaf < f.tvjdn && return zero.((jmax, vcmax))

    if v.tleaf < f.tvjup
        (jmax, vcmax) .* ((v.tleaf - f.tvjdn) / (f.tvjup - f.tvjdn))
    else
        jmax, vcmax
    end
end


"""
    PotentialModifiedFlux(flux, potential_model)

Flux model modified by non-stomatal potential dependence.

Mdifying both electron flux and rubisco activity by the result of 
[`non_stomatal_potential_dependence`](@ref) for `potential_model`.

$(FIELDDOCTABLE)
"""
@default_kw struct PotentialModifiedFlux{F,P} <: AbstractFlux
    flux_model::F       | Flux()
    potential_model::P  | ZhouPotentialDependence()
end

flux_model(m::PotentialModifiedFlux) = m.flux_model

flux(f::PotentialModifiedFlux, v) = begin
    jmax, vcmax = flux(flux_model(f), v)
    pd = non_stomatal_potential_dependence(f.potential_model, v.swp)
    jmax * pd, vcmax * pd
end


"Electron flux parameters"
abstract type AbstractJmax end
 
@columns struct Jmax{μMoM2S,JMoK,JMo} <: AbstractJmax
    jmax25::μMoM2S  | 184.0    | μmol*m^-2*s^-1 | Gamma(100, 184/100)    | (0.0, 1000.0)    | _
    delsj::JMoK     | 640.02   | J*mol^-1*K^-1  | Gamma(100, 640.02/100) | (0.0, 1000.0)    | "DELTAS in Medlyn et al. (2002)"
    eavj::JMo       | 37259.0  | J*mol^-1       | Gamma(100, 37259/100)  | (0.0, 100000.0)  | "Ha     in Medlyn et al. (2002)"
    edvj::JMo       | 200000.0 | J*mol^-1       | Gamma(100, 200000/100) | (0.0, 1000000.0) | "Hd     in Medlyn et al. (2002)"
end

abstract type AbstractVcmax end

"@Vcmax mixin macro adds vcmax fields to a struct"
@mix @columns struct Vcmax{μMoM2S,JMo}
    vcmax25::μMoM2S | 110.0    | μmol*m^-2*s^-1 | Gamma(100, 114/100)   | (0.0, 200.0)    | _
    eavc::JMo       | 47590.0  | J*mol^-1       | Gamma(100, 47590/100) | (0.0, 100000.0) | "Ha    in Medlyn et al. (2002)"
end

@Vcmax struct NoOptimumVcmax{} <: AbstractVcmax end
@Vcmax struct OptimumVcmax{JMo,JMoK} <: AbstractVcmax
    edvc::JMo       | 1.0      | J*mol^-1       | Gamma(10, 1/10)        | (0.0, 10.0)    | "Hd in Medlyn et al. (2002)"
    delsc::JMoK     | 629.26   | J*mol^-1*K^-1  | Gamma(100, 629.26/100) | (0.0, 2000.0)  | "DELTAS in Medlyn et al. (2002)"
end

abstract type AbstractVcJmax end

@default_kw struct VcJmax{J,V} <: AbstractVcJmax where {J<:AbstractVcmax,V<:AbstractJmax}
    jmaxformulation::J  | Jmax()
    vcmaxformulation::V | NoOptimumVcmax()
end
@columns struct DukeVcJmax{J,V,K} <: AbstractVcJmax where {J<:AbstractVcmax,V<:AbstractJmax}
    jmaxformulation::J  | Jmax()           | _  | _ | _ | _
    vcmaxformulation::V | NoOptimumVcmax() | _  | _ | _ | _
    tvjup::K            | 283.15           | K  | _ | (250.0, 350.0) | _
    tvjdn::K            | 273.15           | K  | _ | (250.0, 350.0) | _
end


""" 
    max_electron_transport_rate(f::DukeVcJmax, v, p)
Wrapper to `max_electron_transport_rate()` allowing Vcmax to be forced linearly to zero at low T """
max_electron_transport_rate(f::DukeVcJmax, v, p) = begin
    v.tleaf < f.tvjdn && return zero(f.jmax25)
    jmax = max_electron_transport_rate(f.jmaxformulation, v, p)
    v.tleaf < f.tvjup && (v.tleaf - f.tvjdn) / (f.tvjup - f.tvjdn) * jmax
    jmax
end

""" 
    max_electron_transport_rate(f::VcJmax, v, p)
Wrapper for `max_electron_transport_rate() that simply runs `max_electron_transport_rate` for `jmaxformulation`
"""
max_electron_transport_rate(f::VcJmax, v, p) = max_electron_transport_rate(f.jmaxformulation, v, p)

""" 
    max_electron_transport_rate(f::Jmax, v, p)
Calculates the potential max_electron transport rate (Jmax) at the leaf temperature 
"""
max_electron_transport_rate(f::Jmax, v, p) = begin
    tleafK = v.tleaf |> K
    K25 = 25.0°C |> K
    f.jmax25 * exp((tleafK - K25) * f.eavj / (R * tleafK * K25)) *
    (1 + exp((f.delsj * K25 - f.edvj) / (R * K25))) /
    (1 + exp((f.delsj * tleafK - f.edvj) / (R * tleafK)))
end

"""
Calculates the maximum Rubisco activity (Vcmax) at the leaf temperature.

There is still disagreement as to whether this function has an optimum or not. 
Both versions are well-behaved for tleaf < 0.0
"""
function max_rubisco_activity() end

""" Function allowing Vcmax to be forced linearly to zero at low T.  
Introduced for Duke data. """
function max_rubisco_activity(f::DukeVcJmax, v, p)
    v.tleaf < f.tvjdn && return zero(f.vcmax25)
    vcmax = max_rubisco_activity(f.vcmaxformulation, v, p)
    v.tleaf < f.tvjup && (v.tleaf - f.tvjdn) / (f.tvjup - f.tvjdn) * vcmax
    vcmax
end

""" Wrapper with no alterations to vcmax
Runs max_rubisco_activity for `vcmaxformlation`"""
max_rubisco_activity(f::VcJmax, v, p) = max_rubisco_activity(f.vcmaxformulation, v, p)

""" Vcmax forulation with no optimum"""
function max_rubisco_activity(f::NoOptimumVcmax, v, p)
    tleafK = v.tleaf |> K
    K25 = K(25°C)
    f.vcmax25 * exp((f.eavc * (v.tleaf - K25)) / (K25 * R * tleafK))
end
""" Vcmax formulation with optimum"""
function max_rubisco_activity(f::OptimumVcmax, v, p)
    tleafK = v.tleaf |> K
    K25 = K(25°C)
    f.vcmax25 * exp((v.tleaf - K25) * f.eavc / (R * tleafK * K25)) *
    (1.0 + exp((f.delsc * K25 - f.edvc) / (R * K25))) /
    (1.0 + exp((f.delsc * tleafK - f.edvc) / (R * tleafK)))
end


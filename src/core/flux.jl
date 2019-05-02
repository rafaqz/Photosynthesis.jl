"Electron flux parameters"
abstract type AbstractJmax end
 
@columns struct Jmax{μMoM2S,JMoK,JMo} <: AbstractJmax
    jmax25::μMoM2S  | 184.0    | μmol*m^-2*s^-1 | Gamma(100, 184/100)    | (0.0, 1000.0)    | "Maximum rate of electron transport at 25° C"
    delsj::JMoK     | 640.02   | J*mol^-1*K^-1  | Gamma(100, 640.02/100) | (0.0, 1000.0)    | "DELTAS in Medlyn et al. (2002)"
    eavj::JMo       | 37259.0  | J*mol^-1       | Gamma(100, 37259/100)  | (0.0, 100000.0)  | "Ha in Medlyn et al. (2002)"
    edvj::JMo       | 200000.0 | J*mol^-1       | Gamma(100, 200000/100) | (0.0, 1000000.0) | "Hd in Medlyn et al. (2002)"
end

abstract type AbstractVcmax end

@mix @columns struct Vcmax{μMoM2S,JMo}
    vcmax25::μMoM2S | 110.0    | μmol*m^-2*s^-1 | Gamma(100, 114/100)   | (0.0, 200.0)    | "Maximumrate rate of rubisco activity at 25° C"
    eavc::JMo       | 47590.0  | J*mol^-1       | Gamma(100, 47590/100) | (0.0, 100000.0) | "Ha Medlyn et al. (2002)"
end

@Vcmax struct NoOptimumVcmax{} <: AbstractVcmax end
@Vcmax struct OptimumVcmax{JMo,JMoK} <: AbstractVcmax
    edvc::JMo       | 1.0      | J*mol^-1       | Gamma(10, 1/10)        | (0.0, 10.0)    | "Hd in Medlyn et al. (2002)"
    delsc::JMoK     | 629.26   | J*mol^-1*K^-1  | Gamma(100, 629.26/100) | (0.0, 2000.0)  | "DELTAS in Medlyn et al. (2002)"
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
Calculates the maximum Rubisco activity (Vcmax) at the leaf temperature.

There is still disagreement as to whether this function has an optimum or not. 
Both versions are well-behaved for tleaf < 0.0
"""
function max_rubisco_activity end

""" Vcmax forulation with no optimum"""
function max_rubisco_activity(f::NoOptimumVcmax, tleaf)
    tleafK = tleaf |> K
    K25 = K(25°C)
    f.vcmax25 * exp((f.eavc * (tleaf - K25)) / (K25 * R * tleafK))
end
""" Vcmax formulation with optimum"""
function max_rubisco_activity(f::OptimumVcmax, tleaf)
    tleafK = tleaf |> K
    K25 = K(25°C)
    f.vcmax25 * exp((tleaf - K25) * f.eavc / (R * tleafK * K25)) *
    (1.0 + exp((f.delsc * K25 - f.edvc) / (R * K25))) /
    (1.0 + exp((f.delsc * tleafK - f.edvc) / (R * tleafK)))
end


"""
Jmax and Vcmax are often modified by the same function, so we group them 
inside AbstractFlux wrappers. 
"""
abstract type AbstractFlux end

@udefault_kw struct Flux{J,V} <: AbstractFlux 
    jmaxformulation::J  | Jmax()
    vcmaxformulation::V | NoOptimumVcmax()
end
@columns struct DukeFlux{F,K} <: AbstractFlux
    flux::F             | Flux()           | _  | _ | _              | _
    tvjup::K            | 283.15           | K  | _ | (250.0, 350.0) | _
    tvjdn::K            | 273.15           | K  | _ | (250.0, 350.0) | _
end

@default_kw struct PotentialModifiedFlux{F,P} <: AbstractFlux
    flux::F             | Flux()           
    potential_model::P  | ZhouPotentialDependence()
end

"""
    flux(f, v)
Run jamx and vcmax formulations, returning a 2-tuple
"""
flux(f::Flux, v) = 
    max_electron_transport_rate(f.jmaxformulation, v.tleaf),  
    max_rubisco_activity(f.vcmaxformulation, v.tleaf)

""" 
    flux(f::DukeFlux, v)
Wrapper to `flux()` allowing Jmax and Vcmax to be forced linearly to zero at low T 
"""
flux(f::DukeFlux, v) = begin
    jmax, vcmax = flux(f.flux, v)
    v.tleaf < f.tvjdn && return zero.((jmax, vcmax))

    if v.tleaf < f.tvjup
        (jmax, vcmax) .* ((v.tleaf - f.tvjdn) / (f.tvjup - f.tvjdn))
    else
        jmax, vcmax
    end
end

flux(f::PotentialModifiedFlux, v) = begin
    jmax, vcmax = flux(f.flux, v)
    pd = non_stomatal_potential_dependence(f.potential_model, v.swp)
    # println((pd, v.swp))
    jmax * pd, vcmax * pd
end


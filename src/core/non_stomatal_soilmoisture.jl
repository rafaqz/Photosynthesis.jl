
abstract type AbstractPotentialDependence end

"""
    non_stomatal_potential_dependence (f::AbstractPotentialDependence, swp)

Calculate dependence on soil water potential outside of stomatal effects.
"""
function non_stomatal_potential_dependence end

"""
    NoPotentialDependence()

Parameterless model where soil moisture has no non-stomatal effects.
"""
struct NoPotentialDependence{} <: AbstractPotentialDependence end

"""
    non_stomatal_potential_dependence(f::NoPotentialDependence, swp)

Simply returns 1.0
"""
non_stomatal_potential_dependence(f::NoPotentialDependence, swp) = kneunit(swp) / oneunit(swp)

"""
    LinearPotentialDependence(vpara, vparb)

Simple linear dependance. `vpara` is the swp where vcmax and 
jmax are set to zero, while at `vparb` fluxes are at their maximum rate.
"""
@columns struct LinearPotentialDependence{Pa} <: AbstractPotentialDependence
    vpara::Pa | -300.0 | kPa | (0.0, -2000.0) | _
    vparb::Pa | -100.0 | kPa | (0.0, -2000.0) | _
end

"""
    non_stomatal_potential_dependence(f::LinearPotentialDependence, swp)

Applies LinearPotentialDependence, returning a Float64
"""
non_stomatal_potential_dependence(f::LinearPotentialDependence, swp) =
    if soilwaterpotential < f.vpara
        zero(oneunit(f.vpara) / oneunit(f.vpara))
    elseif soilwaterpotential > f.vparb
        oneunit(f.vpara) / oneunit(f.vpara)
    else
        (soilwaterpotential - f.vpara) / (f.vparb - f.vpara)
    end

"""
    ZhouPotentialDependence(s, ψ)

Parameters following Zhou, et al. Agricultural and Forest Meteorology, 2013.
"""
@columns struct ZhouPotentialDependence{S,Ψ} <: AbstractPotentialDependence
    s::S | 2.0  | MPa^-1 | (0.4, 12.0)  | "Sensitivity parameter indicating the steepness of the decline"
    ψ::Ψ | -1.0 | MPa    | (-0.1, -4.0) | "The water potential at which f(Ψpd) decreases to half of its maximum value"
end

"""
    non_stomatal_potential_dependence(f::ZhouPotentialDependence, swp)

Zhou, et al. Agricultural and Forest Meteorology. 2013.0
"""
non_stomatal_potential_dependence(f::ZhouPotentialDependence, swp) =
    (1 + exp(f.s * f.ψ)) / (1 + exp(f.s * (f.ψ - swp)))

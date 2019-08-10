
abstract type AbstractPotentialDependence end

@columns struct LinearPotentialDependence{Pa} <: AbstractPotentialDependence
    vpara::Pa | -300.0 | kPa | (0.0, -2000.0) | _
    vparb::Pa | -100.0 | kPa | (0.0, -2000.0) | _
end

@columns struct ZhouPotentialDependence{S,Ψ} <: AbstractPotentialDependence
    s::S | 2.0  | MPa^-1 | (0.4, 12.0)  | "Sensitivity parameter indicating the steepness of the decline"
    ψ::Ψ | -1.0 | MPa    | (-0.1, -4.0) | "The water potential at which f(Ψpd) decreases to half of its maximum value"
end

struct NoPotentialDependence{} <: AbstractPotentialDependence end


"""
    non_stomatal_potential_dependence(f::NoPotentialDependence, v)
Returns 1.0
"""
non_stomatal_potential_dependence(f::NoPotentialDependence, swp) = 1.0

"""
    non_stomatal_potential_dependence(f::LinearPotentialDependence, v, p)
Simple linear dependance. At vpara is the swp where vcmax and 
jmax are set to zero, at vparb fluxes are at their maximum rate.
"""
non_stomatal_potential_dependence(f::LinearPotentialDependence, swp) =
    if soilwaterpotential < f.vpara
        0.0
    elseif soilwaterpotential > f.vparb
        1.0
    else
        (soilwaterpotential - f.vpara) / (f.vparb - f.vpara)
    end

"""
    non_stomatal_potential_dependence(f::ZhouPotentialDependence, v, p)
Zhou, et al. Agricultural and Forest Meteorology. 2013.0
"""
non_stomatal_potential_dependence(f::ZhouPotentialDependence, swp) =
    (1 + exp(f.s * f.ψ)) / (1 + exp(f.s * (f.ψ - swp)))

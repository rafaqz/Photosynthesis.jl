"""
Rubisco regeneration models. 
These are run in the [`rubisco_regeneration`](@ref) method.
"""
abstract type AbstractRubiscoRegen end

""" 
    rubisco_regeneration(f::RubiscoRegen, v)

Returns RuBP regeneration rate, in `u"umol*m-2*s-1"`
for formulation `f` given variables `v`.
"""
function rubisco_regeneration end

"""
    RubiscoRegen(theta, ajq)

Rubisco regeneration model. 

TODO: specify origin of formulation

$(FIELDDOCTABLE)
"""
@columns struct RubiscoRegen{} <: AbstractRubiscoRegen
    theta::Float64 | 0.4    | _ | (0.0, 1.0) | "Shape parameter of the non-rectangular hyperbola"
    ajq::Float64   | 0.324  | _ | (0.0, 1.0) | "Quantum yield of electron transport"
end

rubisco_regeneration(f::RubiscoRegen, v) = begin
    a = f.theta
    b = -(f.ajq * v.par + v.jmax)
    c = f.ajq * v.par * v.jmax
    j = quad(Lower(), a, b, c) # Actual e- transport rate, umol m-2 s-1
    return j / 4
end

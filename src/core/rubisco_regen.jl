abstract type AbstractRubiscoRegen end
@columns struct RubiscoRegen{} <: AbstractRubiscoRegen
    theta::Float64  | 0.4    | _              | Beta(8, 12)  | (0.0, 1.0) | "Shape parameter of the non-rectangular hyperbola."
    ajq::Float64    | 0.324  | _              | Beta(4, 8.3) | (0.0, 1.0) | "Quantum yield of electron transport."
end

""" 
    rubisco_regeneration(f::RubiscoRegen, v)
RuBP-regen rate, in umol m-2 s-1 """
function rubisco_regeneration(f::RubiscoRegen, v)
    a = f.theta
    b = -(f.ajq * v.par + v.jmax)
    c = f.ajq * v.par * v.jmax
    j = quad(Lower(), a, b, c) # Actual e- transport rate, umol m-2 s-1
    return j / 4
end

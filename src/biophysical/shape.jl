abstract type AbstractGSShape end

struct HardMinimumGS <: AbstractGSShape end

@columns struct HyperbolicMinimumGS <: AbstractGSShape
    hmshape::Float64 | 0.999 | _ | Beta(20, 1.01) | (0.9, 1.0) | _
end

function shape_gs(f::HardMinimumGS, v, p) 
    gs = p.g0 + v.gsdiva * v.aleaf
    max(p.g0, gs)
end

""" 
Hyperbolic minimum; use for optimization to avoid discontinuity.
`hmshape` determines how smooth the transition is.
"""
function shape_gs(f::HyperbolicMinimumGS, v, p) 
    aleafhypmin = (v.ac + v.aj - sqrt((v.ac + v.aj)^2 - 4f.hmshape * v.ac * v.aj)) /
              (2 * f.hmshape) - v.rd
    p.g0 + v.gsdiva * aleafhypmin
end

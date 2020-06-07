abstract type StomatalConductanceShape end

struct HardMinimum <: StomatalConductanceShape end

@columns struct HyperbolicMinimum <: StomatalConductanceShape
    hmshape::Float64 | 0.999 | _ | (0.9, 1.0) | _
end

function shape_gs(f::HardMinimum, v, p) 
    gs = p.g0 + v.gsdiva * v.aleaf
    max(p.g0, gs)
end

""" 
Hyperbolic minimum; use for optimization to avoid discontinuity.
`hmshape` determines how smooth the transition is.
"""
function shape_gs(f::HyperbolicMinimum, v, p) 
    aleafhypmin = (v.ac + v.aj - sqrt((v.ac + v.aj)^2 - 4f.hmshape * v.ac * v.aj)) /
              (2 * f.hmshape) - v.rd
    p.g0 + v.gsdiva * aleafhypmin
end

"""
Determines the shape of stomatal conductance as it approaches zero,
run in [`shape_gs`](@ref) method.
""" 
abstract type StomatalConductanceShape end

"""
    shape_gs(f::StomatalConductanceShape, v, p)

Function to determine the hyperbolic minimum shape,
for a formulation `f` in [`StomatalConductanceShape`](@ref),
fiven variables `v` and stomatal conductance parameters `p`.

Returns the stomatal conductance `gs` in `u"mol*m^-2*s^-1"`.
"""
function shape_gs end

""" 
    HardMinimum()

Return stomatal conductance value that abrutly
switches to the minimum value `g0`.
"""
struct HardMinimum <: StomatalConductanceShape end

shape_gs(f::HardMinimum, v, p) = begin
    gs = p.g0 + v.gs_div_a * v.aleaf
    max(p.g0, gs)
end

""" 
    HyperbolicMinimum()

Return stomatal conductance value that switches to the minimum value `g0`
using a hyperbolic curve, used for optimization to avoid discontinuity. 
`hmshape` determines how smooth the transition is.

$(FIELDDOCTABLE)
"""
@columns struct HyperbolicMinimum <: StomatalConductanceShape
    hmshape::Float64 | 0.999 | _ | (0.9, 1.0) | _
end

shape_gs(f::HyperbolicMinimum, v, p) = begin
    aleafhypmin = (v.ac + v.aj - sqrt((v.ac + v.aj)^2 - 4f.hmshape * v.ac * v.aj)) /
              (2 * f.hmshape) - v.rd
    p.g0 + v.gs_div_a * aleafhypmin
end

abstract type AbstractStomatalConductance end

" @Gs mixin macro adds stomatal conductance fields to a struct"
@mix @columns struct Gs{μMoMo,F}
    gamma::μMoMo    | 0.0    | μmol*mol^-1    | Gamma(1, 2)     | (0.0, 10.0) | "Gamma for all Ball-Berry type models"
    g1::F           | 7.0    | _              | Gamma(10, 7/10) | (0.0, 10.0) | "Slope parameter"
end

"""
    stomatal_conductance!(v, ::AbstractPhotoModel, p)
Stomatal conductance and intercellular CO2 partial pressure calculations.

v.aleaf is NET leaf photosynthesis.
"""
function stomatal_conductance!(f, v, p) end


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
    aleafhypmin = (ac + aj - sqrt((v.ac + v.aj)^2 - 4f.hmshape * v.ac * v.aj)) /
              (2 * f.hmshape) - v.rd
    p.g0 + v.gsdiva * aleafhypmin
end


function assimilation!(v, p)
    v.ac = rubisco_limited_rate(p.model, v, p)
    v.aj = transport_limited_rate(p.model, v, p)
    v.aleaf = min(v.ac, v.aj) - v.rd
end

function rubisco_limited_rate end

function transport_limited_rate end

"""
    gsdiva(::AbstractModelGS, v, p)
Formulation-specific component for the Ball-Berry family of stomatal conductance models.
"""
function gsdiva end


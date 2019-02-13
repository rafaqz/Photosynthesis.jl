
""" Three parameter Ball-Berry stomatal conductance formulation parameters
Has the extra `vpdmin` paramater in Pa and gk scalar parameter.
(modelgs = 5 in maestra)
"""
@Gs struct ThreeParStomatalConductance{F,Pa} <: AbstractStomatalConductance
    gk::F      | 0.3    | _  | Gamma(2, 0.3/2)    | _ | _
    vpdmin::Pa | 1500.0 | Pa | Gamma(10, 1500/10) | _ | _
end

"""
    gsdiva(::ThreeParStomatalConductance, v) """
gsdiva(f::ThreeParStomatalConductance, v) =
    f.g1 / (v.cs - f.gamma) / (max(f.vpdmin, v.vpdleaf)Pa^-1)^f.gk

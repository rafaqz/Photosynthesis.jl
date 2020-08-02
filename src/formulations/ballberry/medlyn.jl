
"""
    MedlynStomatalConductanceSubModel(vpdmin, gk, gamma, g1)

Medlyn stomatal conductance formulation parameters
Has the extra `vpdmin` paramater in Pa.
(modelgs = 4 in maestra)

$(FIELDDOCTABLE)
"""
@MixinBallBerryStomatalConductanceSubModel struct MedlynStomatalConductanceSubModel{Pa} <: AbstractBallBerryStomatalConductanceSubModel
    vpdmin::Pa | 1500.0 | kPa | _ | _
    gk::F      | 0.3    | _   | _ | _
end

gs_div_a(f::MedlynStomatalConductanceSubModel, v) =
    (oneunit(f.g1) + f.g1 / (max(f.vpdmin, v.vpdleaf)/MPa)^(1 - f.gk)) / v.cs

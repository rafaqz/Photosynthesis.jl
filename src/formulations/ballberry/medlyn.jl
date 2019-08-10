
"""
Medlyn stomatal conductance formulation parameters
Has the extra `vpdmin` paramater in Pa.
(modelgs = 4 in maestra)
"""
@MixinBallBerryGs struct MedlynGSsubModel{Pa} <: AbstractGSsubModel
    vpdmin::Pa | 1500.0 | kPa | _ | _
    gk::F      | 0.3    | _   | _ | _
end

"""
    gsdiva(::MedlynStomatalConductance, v) """
gsdiva(f::MedlynGSsubModel, v) =
    (oneunit(f.g1) + f.g1) / (max(f.vpdmin, v.vpdleaf)/kPa)^(1 - f.gk) / v.cs

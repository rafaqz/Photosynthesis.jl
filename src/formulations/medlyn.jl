
"""
Medlyn stomatal conductance formulation parameters
Has the extra `vpdmin` paramater in Pa.
(modelgs = 4 in maestra)
"""
@Gs struct MedlynStomatalConductance{Pa} <: AbstractStomatalConductance
    vpdmin::Pa | 1500.0 | Pa | Gamma(10, 1500/10) | _ | _
end

"""
    gsdiva(::MedlynStomatalConductance, v, p) """
gsdiva(f::MedlynStomatalConductance, v, p) =
    (1.0 + (f.g1 * v.fsoil) / sqrt(max(f.vpdmin, v.vpdleaf)Pa^-1)) / v.cs

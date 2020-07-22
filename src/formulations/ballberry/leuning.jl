"""
    LeuningStomatalConductanceSubModel(d0l, gamma, g1)

Leuning stomatal conductance formulation with extra `d0l` paramater.

From R. Leuning, A critical appraisal of a combined stomatal-photosynthesis
model for C3 plants. Plant, Celt and Environment (1995) 18, 339-355

$(FIELDDOCTABLE)
"""
@MixinBallBerryStomatalConductanceSubModel struct LeuningStomatalConductanceSubModel{Pa} <: AbstractBallBerryStomatalConductanceSubModel
    d0l::Pa    | 1500.0 | Pa | (0.0, 2000.0) | _
end

gs_div_a(f::LeuningStomatalConductanceSubModel, v) = 
    f.g1 / ( v.cs - f.gamma) / (1.0 + v.vpdleaf / f.d0l)

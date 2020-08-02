"""
    LeuningStomatalConductanceSubModel(D0, gamma, g1)

Leuning stomatal conductance formulation with extra `D0` paramater.

From R. Leuning, A critical appraisal of a combined stomatal-photosynthesis
model for C3 plants. Plant, Celt and Environment (1995) 18, 339-355

$(FIELDDOCTABLE)
"""
@MixinBallBerryStomatalConductanceSubModel struct LeuningStomatalConductanceSubModel{D0} <: AbstractBallBerryStomatalConductanceSubModel
    D0::D0    | 1500.0 | kPa | (0.0, 2000.0) | _
end

gs_div_a(f::LeuningStomatalConductanceSubModel, v) = 
    f.g1 / (v.cs - f.gamma) / (1 + v.vpdleaf / f.D0)

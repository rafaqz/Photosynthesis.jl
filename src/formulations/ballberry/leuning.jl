"""
Leuning stomatal conductance formulation
Has the extra `d0l` paramater in Pa.
"""
@MixinBallBerryGs struct LeuningGSsubModel{Pa} <: AbstractGSsubModel
    d0l::Pa    | 1500.0 | Pa | Gamma(10, 1500/10) | (0.0, 2000.0) | _
end

"""
    gsdiva(::LeuningGSsubModel, v) 
From R. Leuning, A critical appraisal of a combined stomatal-photosynthesis
model for C3 plants. Plant, Celt and Environment (1995) 18, 339-355
"""
gsdiva(f::LeuningGSsubModel, v) =
    f.g1 / ( v.cs - f.gamma) / (1.0 + v.vpdleaf / f.d0l)

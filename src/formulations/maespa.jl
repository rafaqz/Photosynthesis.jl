abstract type AbstractMaespaModel{GS,SM} <: AbstractBallBerryModel{GS,SM} end

" @Maespa mixin macro adds base maespa fields to a struct"
@BB @mix struct Maespa{mMoM2S,M2SMPaMo,SH}
    plantk::mMoM2S       | 3.0             | mmol*m^-2*s^-1*MPa^-1 | (0.0, 10.0) | _ | _
    totsoilres::M2SMPaMo | 0.5             | m^2*s^1*MPa^1*mmol^-1 | (0.0, 10.0) | _ | _
    gsshape::SH          | HardMinimumGS() | _                     | _ | _ | _
end

@mix @default_kw struct MaespaSoil{T<:AbstractPotentialDependence}
    soildata::T | LinearPotentialDependence()
end

"@MaespaVars mixin macro"
@Vars @mix struct MaespaVars{mMoM2SMPa,kPa}
    ktot::mMoM2SMPa  | 2.0         | mmol*m^-2*s^-1*MPa^-1 | _
    weightedswp::kPa | 0.0         | kPa                   | _
    psil::kPa        | -111.0      | kPa                   | _
end

photo_init!(f::AbstractMaespaModel, v) = v.ktot = 10 / (f.totsoilres + 1.0 / f.plantk)


function psil_gs!(f, v, p) end

"""
    stomatal_conductance!(f::TuzetModel, v, p)
Stomatal conductance calculations for the Jarvis model
"""
function stomatal_conductance!(f::AbstractMaespaModel, v, p)
    v.gsdiva = gsdiva(p.model.gsmodel, v, p)
    assimilation!(v, p)

    v.gs = shape_gs(f.gsshape, v, f)
    psil_gs!(f, v, p)

    if v.gs > zero(v.gs) && v.aleaf > zero(v.aleaf)
        v.ci = v.cs - v.aleaf / v.gs
    else
        v.ci = v.cs
    end

    nothing
end


"""
Dependance of assimilation on soil water potential (SWP).
Modifier function (return value between 0.0 and 1.0) 
"""
function vjmax_water end

"""
    vjmax_water(f::NoPotentialDependence, v)
Returns 1.0
"""
vjmax_water(f::NoPotentialDependence, v) = 1.0

"""
    vjmax_water(f::LinearPotentialDependence, v, p)
Simple linear dependance.
vpara is swp where vcmax is zero, vparb is swp where vcmax is one.
"""
vjmax_water(f::LinearPotentialDependence, v) =
    if v.weightedswp < f.vpara 
        0.0
    elseif v.weightedswp > f.vparb 
        1.0
    else 
        (v.weightedswp - f.vpara) / (f.vparb - f.vpara)
    end

"""
    vjmaxw(f::ZhouPotentialDependence, v, p)
Zhou, et al. Agricultural and Forest Meteorology. 2013.0
"""
vjmax_water(f::ZhouPotentialDependence, v) = 
    (1 + exp(f.vpara * f.vparb)) / (1 + exp(f.vpara * (f.vparb - v.weightedswp)))

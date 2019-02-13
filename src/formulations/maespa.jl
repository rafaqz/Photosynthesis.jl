# Abstract Maespa model. Tuzet and Emax models inherit from this

abstract type AbstractMaespaModel <: AbstractBallBerryModel end

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


function psil_gs!(f, v) end

"""
    stomatal_conductance!(f::TuzetModel, v, p)
Stomatal conductance calculations for the Jarvis model
"""
function stomatal_conductance!(f::AbstractMaespaModel, v)
    v.gsdiva = gsdiva(f.gsmodel, v)
    assimilation!(v, f)

    v.gs = shape_gs(f.gsshape, v, f)
    psil_gs!(f, v)

    if v.gs > zero(v.gs) && v.aleaf > zero(v.aleaf)
        v.ci = v.cs - v.aleaf / v.gs
    else
        v.ci = v.cs
    end

    nothing
end

# Maespa abstract/mixin implementations. Tuzet and Emax models inherit from these.


# Maespa abstract/mixin implementation of photosynthesis

abstract type AbstractMaespaPhotosynthesis <: AbstractBallBerryPhotosynthesis end

" @Maespa mixin macro adds base maespa fields to a struct"
@MixinBallBerryPhoto @mix struct MixinMaespaPhoto{mMoM2S,M2SMPaMo,SH}
    plantk::mMoM2S       | 3.0             | mmol*m^-2*s^-1*MPa^-1 | (0.0, 10.0) | _ | _
    totsoilres::M2SMPaMo | 0.5             | m^2*s^1*MPa^1*mmol^-1 | (0.0, 10.0) | _ | _
    gsshape::SH          | HardMinimumGS() | _                     | _ | _ | _
end

@mix @udefault_kw struct MixinMaespaSoil{T,NS}
    soilmethod::T    | ConstantSoilMethod()
    non_stomatal::NS | ZhouPotentialDependence()
end

photo_init!(f::AbstractMaespaPhotosynthesis, v) = v.ktot = 10 / (f.totsoilres + 1.0 / f.plantk)


"@MaespaVars mixin macro"
@MixinFvCBVars @mix struct MixinMaespaVars{mMoM2SMPa,kPa}
    ktot::mMoM2SMPa  | 2.0         | mmol*m^-2*s^-1*MPa^-1 | _
    weightedswp::kPa | 0.0         | kPa                   | _
    psil::kPa        | -111.0      | kPa                   | _
end

function psil_gs!(f, v) end

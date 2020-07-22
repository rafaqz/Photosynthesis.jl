"""
Abstract supertype for component to specify the source data for soil moisture.
"""
abstract type AbstractSoilData end

"""
Soil water data is measure as moisture-deficit.
"""
abstract type AbstractDeficitSoilData <: AbstractSoilData end

" @Deficit mixin macro adds water deficit fields to a struct"
@mix @default_kw struct MixinDeficit{T}
    smd1::T   | 1.0
    smd2::T   | 1.0
end


"""
Soil water data is measure by volumetric content.
"""
abstract type AbstractContentSoilData <: AbstractDeficitSoilData end

"@Content mixin macro adds water content fields to a struct"
@MixinDeficit @mix struct MixinContent{W}
    swmin::W  | 0.0
    swmax::W  | 1.0
end

"""
    DeficitSoilData(smd1, smd2)

Soil water is measured as deficit
"""
@MixinDeficit struct DeficitSoilData{} <: AbstractDeficitSoilData end

"""
    ContentSoilData(swmin, swmax, smd1, smd2)

Soil water is measured as volumetric content
"""
@MixinContent struct ContentSoilData{} <: AbstractContentSoilData end

"""
    SimulatedSoilData(swmin, swmax, smd1, smd2)

Simulated soil volumetric content.
"""
@MixinContent struct SimulatedSoilData{} <: AbstractContentSoilData end



"""
Methods for calculating the effect of soil moisture on stomatal conductance. 
"""
abstract type AbstractSoilMethod end

"""
    soilmoisture_conductance(soilmethod::AbstractSoilMethod, vars)

Calculate the effect of soil water content on stomatal conductance

These methods dispatch on AbstractSoilMethod or AbstractSoilData types.
"""
function soilmoisture_conductance end


"""
    VolumetricSoilMethod(soildata, wc1, wc2, soilroot, soildepth)

Soil method where soil water is handled volumetrically.
"""
@default_kw struct VolumetricSoilMethod{WC,R,D,T<:AbstractContentSoilData} <: AbstractSoilMethod
    soildata::T  | ContentSoilData()
    wc1::WC      | 0.5
    wc2::WC      | 0.5
    soilroot::R  | 0.5
    soildepth::D | 0.5
end

soilmoisture_conductance(f::VolumetricSoilMethod{<:SimulatedSoilData}, v) =
    v.soilmoist = f.soilroot / f.soildepth # TODO fix this, where is fsoil set??

soilmoisture_conductance(f::VolumetricSoilMethod{<:AbstractContentSoilData}, v) = begin
    # TODO: check in constructor f.wc1 > f.wc2 && error("error: wc1 needs to be smaller than wc2.")
    fsoil = -f.wc1 / (f.wc2 - f.wc1) + v.soilmoist / (f.wc2 - f.wc1)
    max(min(fsoil, 1.0), 0.0)
end


"""
    NoSoilMethod()

No soil method is used, soil moisture has no effect.
"""
struct NoSoilMethod <: AbstractSoilMethod end

soilmoisture_conductance(f::NoSoilMethod, v) = 1.0

"""
Soil method where soil water is held constant.
"""
struct ConstantSoilMethod <: AbstractSoilMethod end

"""
    PotentialSoilMethod(soildata, swpexp)

Soil mmodel where soil water is measured as water potential.
"""
@udefault_kw @units @bounds @description struct PotentialSoilMethod{E} <: AbstractSoilMethod
    swpexp::E | 0.5 | kPa^-1 | (0.0, 10.0) | "Exponent for soil water-potential response of stomata"
end

soilmoisture_conductance(f::PotentialSoilMethod, v) = exp(f.swpexp * v.swp)


"""
    DeficitSoilMethod(soildata)

Soil method where soil water is handled as a moisture deficit.
"""
@default_kw struct DeficitSoilMethod{T<:AbstractDeficitSoilData} <: AbstractSoilMethod
    soildata::T | ContentSoilData()
end

"Convert volumetric data to deficit"
soilmoisture_conductance(f::DeficitSoilMethod{<:AbstractContentSoilData}, v) = begin
    s = f.soildata
    soilmoist = (s.swmax - v.soilmoist) / (s.swmax - s.swmin)
    soilmoisture_conductance(s, soilmoist, v)
end

soilmoisture_conductance(f::DeficitSoilMethod{<:DeficitSoilData}, v) =
    soilmoisture_conductance(f.soildata, v.soilmoist, v)

"""
    soilmoisture_conductance(soil::AbstractDeficitSoilData, soilmoist, v)

GranierLoustau 1994 Fs = 1 - a exp(b SMD) where SMD is soil
moisture deficit, dimless

TODO: move smd1/smd2 fields to SoilMethod
"""
soilmoisture_conductance(soil::AbstractDeficitSoilData, soilmoist, v) = begin
    effect = 1.0

    # Exponential relationship with deficit: params SMD1, SMD2
    if soil.smd1 > zero(soil.smd1)
        effect = 1.0 - soil.smd1 / oneunit(soil.smd1) * exp(soil.smd2 * soilmoist / oneunit(soilmoist)^2)
    # Linear decline with increasing deficit: pUT SMD1 < 0, param SMD2
    elseif soil.smd2 > zero(soil.smd2)
        if oneunit(soilmoist) - soilmoist < soil.smd2
            effect = (oneunit(soilmoist) - soilmoist) / soil.smd2
        end
    end
    return max(effect, zero(effect))
end

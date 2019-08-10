abstract type AbstractSoilData end
abstract type AbstractDeficitSoilData <: AbstractSoilData end
abstract type AbstractContentSoilData <: AbstractDeficitSoilData end

" @Deficit mixin macro adds water deficit fields to a struct"
@mix @default_kw struct MixinDeficit{T}
    smd1::T   | 1.0
    smd2::T   | 1.0
end

"@Content mixin macro adds water content fields to a struct"
@MixinDeficit @mix struct MixinContent{W}
    swmax::W  | 1.0
    swmin::W  | 0.0
end

"Soil water is measured as deficit"
@MixinDeficit struct DeficitSoilData{} <: AbstractDeficitSoilData end
"Soil water is measured as volumetric content"
@MixinContent struct ContentSoilData{} <: AbstractContentSoilData end
"Simulated soil volumetric content"
@MixinContent struct SimulatedSoilData{} <: AbstractContentSoilData end
"Soil water is measured as water potential"


abstract type AbstractSoilMethod end

struct NoSoilMethod <: AbstractSoilMethod end

@default_kw struct VolumetricSoilMethod{WC,R,D,T<:AbstractContentSoilData} <: AbstractSoilMethod
    soildata::T  | ContentSoilData()
    wc1::WC      | 0.5
    wc2::WC      | 0.5
    soilroot::R  | 0.5
    soildepth::D | 0.5
end

struct ConstantSoilMethod <: AbstractSoilMethod end

@default_kw struct DeficitSoilMethod{T<:AbstractDeficitSoilData} <: AbstractSoilMethod
    soildata::T | ContentSoilData()
end

@description @bounds @units @udefault_kw struct PotentialSoilMethod{E} <: AbstractSoilMethod
    swpexp::E | 0.5 | kPa^-1 | [0.0, 10.0] | "Exponent for soil water-potential response of stomata"
end



"""
Calculate the effect of soil water content on stomatal conductance

These methods dispatch on AbstractSoilMethod or AbstractSoilData types.
"""
function soilmoisture_conductance end

soilmoisture_conductance(f::NoSoilMethod, v) = 1.0

soilmoisture_conductance(f::PotentialSoilMethod, v) = exp(f.swpexp * v.swp)

"Convert volumetric data to deficit"
soilmoisture_conductance(f::DeficitSoilMethod{<:AbstractContentSoilData}, v) = begin
    s = f.soildata
    soilmoist = (s.swmax - v.soilmoist) / (s.swmax - s.swmin)
    soilmoisture_conductance(s, soilmoist, v)
end

soilmoisture_conductance(f::DeficitSoilMethod{<:DeficitSoilData}, v) =
    soilmoisture_conductance(f.soildata, v.soilmoist, v)

soilmoisture_conductance(f::VolumetricSoilMethod{<:SimulatedSoilData}, v) =
    v.soilmoist = f.soilroot / f.soildepth # TODO fix this, where is fsoil set??

soilmoisture_conductance(f::VolumetricSoilMethod{<:AbstractContentSoilData}, v) = begin
    # TODO: check in constructor f.wc1 > f.wc2 && error("error: wc1 needs to be smaller than wc2.")
    fsoil = -f.wc1 / (f.wc2 - f.wc1) + v.soilmoist / (f.wc2 - f.wc1)
    max(min(fsoil, 1.0), 0.0)
end

"""
    soilmoisture_conductance(soil::AbstractDeficitSoilData, soilmoist, v)
GranierLoustau 1994 Fs = 1 - a exp(b SMD) where SMD is soil
moisture deficit, dimless """
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

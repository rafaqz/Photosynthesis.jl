# Jarvis parameters and functions are kep separately here
# to simplify the rest of the module.
# There are a lot of them, and they may not be used regularly.

abstract type AbstractJarvisLight end

@with_kw mutable struct JarvisLight{T} <: AbstractJarvisLight 
    i0::T = 1.0u"mol*m^-2*s^-1"
end

abstract type AbstractJarvisCO2 end

struct JarvisNoCO2 <: AbstractJarvisCO2 end
@with_kw mutable struct JarvisLinearCO2{T} <: AbstractJarvisCO2 
    gsja::T = 1.0u"μmol^-1*mol"
end
@with_kw mutable struct JarvisNonlinearCO2{T} <: AbstractJarvisCO2 
    gsjb::T = 1.0u"μmol*mol^-1"
end


abstract type AbstractJarvisTemp end

@mix @with_kw struct JarTemp{T}
    tmax::T = 40.0u"°C"
    tref::T = 25.0u"°C"
    t0::T = 0.0u"°C"
end

struct JarvisNoTemp <: AbstractJarvisTemp end
@JarTemp mutable struct JarvisTemp1{} <: AbstractJarvisTemp end
@JarTemp mutable struct JarvisTemp2{} <: AbstractJarvisTemp end


abstract type AbstractJarvisVPD end

"""
Hyperbolic response to vapour pressure deficit. 
Parameters vk1 and vk2 are the dimensionless scalar and exponent.
"""
@with_kw mutable struct JarvisHyperbolicVPD{T} <: AbstractJarvisVPD 
    vk1::T = 1.0
    vk2::T = 1.0
end

"""
Non-linear Lohammer response to vapour pressure deficit.
Parameters vpd1 and vpd2 are in Pascals.
"""
@with_kw mutable struct JarvisLohammerVPD{T} <: AbstractJarvisVPD 
    vpd1::T = 1.0u"Pa"
    vpd2::T = 1.0u"Pa"
end
@with_kw mutable struct JarvisFractionDeficitVPD{T} <: AbstractJarvisVPD 
    vmfd0::T = 1.0u"mmol*mol^-1"
end
@with_kw mutable struct JarvisLinearDeclineVPD{T} <: AbstractJarvisVPD 
    d0::T = 1.0u"Pa"
end


"""
Calculate stomatal conductance gs according to the Jarvis model.

This model calculates gs by multiplying together factors for several
environmental variables: light, VPD, CO2 and temperature.
"""
function factor_conductance(f, v, p)
    flight = min(max(calc_flight(f.lightmethod, v, p), 1.0), 0.0)
    fvpd = min(max(calc_fvpd(f.vpdmethod, v, p), 1.0), 0.0)
    fco2 = min(max(calc_fco2(f.co2method, v, p), 1.0), 0.0)
    ftemp = min(max(calc_ftemp(f.tempmethod, v, p), 1.0), 0.0)

    return (f.gsref - f.gsmin) * flight * fvpd * fco2 * ftemp * v.fsoil + f.gsmin
end


"""
Jarvis stomatal conductance model

Combines factors from soilmethod, co2 method, vpdmethod and tempmethod
to gain an overall stomatal conductance.
"""
@with_kw mutable struct JarvisModel{S<:AbstractSoilMethod, C<:AbstractJarvisCO2, 
                                    V<:AbstractJarvisVPD, L<:AbstractJarvisLight,
                                    T<:AbstractJarvisTemp, MoMeS, mMoMoS
                                   } <: AbstractPhotoModel
    soilmethod::S  = DeficitSoilMethod()
    co2method::C   = JarvisNonlinearCO2()
    vpdmethod::V   = JarvisLohammerVPD()
    lightmethod::L = JarvisLight()
    tempmethod::T  = JarvisTemp2()
    gsmin::MoMeS   = 1.0u"mol*m^-2*s^-1"
    gsref::MoMeS   = 1.0u"mol*m^-2*s^-1"
    vmfd::mMoMoS   = 1.0u"mmol*mol^-1"
end

""" Response to incident radiation in umol m^-2 s^-1 """
calc_flight(f::JarvisLight, v, p) = f.i0 > zero(f.i0) ? v.par / (v.par + f.i0) : 1.0

""" Hyperbolic decline with VPD (VPD in Pa) * BRAY (ALEX BOSC) """
calc_fvpd(f::JarvisHyperbolicVPD, v, p) =
    v.vpdleaf > zero(v.vpdleaf) ? 1.0 / (f.vk1 * v.vpdleaf^f.vk2) : 0.0
""" Lohammer response to VPD (VPDD0 in Pa) * ECOCRAFT """
calc_fvpd(f::JarvisLohammerVPD, v, p) =
    v.vpdleaf >= f.vpd1 ? 1.0 - (v.vpdleaf - f.vpd1) / (f.vpd2 - f.vpd1) : 1.0
""" Mole fraction deficit (VMFD in mmol mol - 1) * MARK RAYMENT """
calc_fvpd(f::JarvisFractionDeficitVPD, v, p) = 
    f.vmfd0 > 0.0 ? 1.0 - v.vmleaf / f.vmfd0 : 1.0
""" Linear decline with VPD (VPD1, VPD2 in Pa) * GROMIT (TIM RANDLE) """
calc_fvpd(f::JarvisLinearDeclineVPD, v, p) = f.d0 > 0.0 ? 1.0 / (1.0 + v.vpdleaf / f.d0) : 1.0

""" No response to temperature in Jarvis stomatal conductance """
calc_ftemp(f::JarvisNoTemp, v, p) = 1.0
calc_ftemp(f::JarvisTemp1, v, p) = begin
    p = (f.tmax - f.tref) / (f.tref - f.t0)
    (v.tleaf - f.t0) * ((f.tmax - v.tleaf)^p) / ((f.tref - f.t0) * ((f.tmax - tref)^p))
end
calc_ftemp(f::JarvisTemp2, v, p) = 
    (v.tleaf - f.t0) * (2 * f.tmax - f.t0 - v.tleaf) / ((f.tref - f.t0) * (2 * f.tmax - f.t0 - f.tref))

""" No influence from CO2 for Jarvis stomatal conductance"""
calc_fco2(f::JarvisNoCO2, v, p) = 1.0
""" Linear response to CO2 for Jarvis stomatal conductance"""
calc_fco2(f::JarvisLinearCO2, v, p) = 
    f.gsja != 0.0 ? 1 - f.gsja * (v.cs - 350.0) : 1.0
""" Non-linear response to CO2 for Jarvis stomatal conductance"""
calc_fco2(f::JarvisNonlinearCO2, v, p) =
    f.gsjb != 0.0 ? (f.gsjb + 350.0) / (f.gsjb + v.cs) : 1.0

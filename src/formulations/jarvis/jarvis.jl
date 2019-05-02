# Jarvis model
# Uses factor combination to detaermine stomatal conductance.

abstract type AbstractJarvisLight end
@columns mutable struct JarvisLight{T} <: AbstractJarvisLight
    i0::T | 1.0 | mol*m^-2*s^-1 | _ | _ | _
end

abstract type AbstractJarvisCO2 end
struct JarvisNoCO2 <: AbstractJarvisCO2 end
@columns mutable struct JarvisLinearCO2{T} <: AbstractJarvisCO2
    gsja::T | 1.0 | μmol^-1*mol | _ | _ | _
end
@columns mutable struct JarvisNonlinearCO2{T} <: AbstractJarvisCO2
    gsjb::T | 1.0 | μmol*mol^-1 | _ | _ | _
end


abstract type AbstractJarvisTemp end

@mix @columns struct JarTemp{T}
    tmax::T | (273.15 + 40.0) | K | _ | _ | _
    tref::T | (273.15 + 25.0) | K | _ | _ | _
    t0::T   | 273.15          | K | _ | _ | _
end

struct JarvisNoTemp <: AbstractJarvisTemp end
@JarTemp mutable struct JarvisTemp1{} <: AbstractJarvisTemp end
@JarTemp mutable struct JarvisTemp2{} <: AbstractJarvisTemp end


abstract type AbstractJarvisVPD end

"""
Hyperbolic response to vapour pressure deficit.
Parameters vk1 and vk2 are the dimensionless scalar and exponent.
"""
@columns mutable struct JarvisHyperbolicVPD{T} <: AbstractJarvisVPD
    vk1::T | 1.0 | _ | _ | _ | _
    vk2::T | 1.0 | _ | _ | _ | _
end

"""
Non-linear Lohammer response to vapour pressure deficit.
Parameters vpd1 and vpd2 are in Pascals.
"""
@columns mutable struct JarvisLohammerVPD{T} <: AbstractJarvisVPD
    vpd1::T  | 1.0 | Pa | _ | _ | _
    vpd2::T  | 1.0 | Pa | _ | _ | _
end
@columns mutable struct JarvisFractionDeficitVPD{T} <: AbstractJarvisVPD
    vmfd0::T | 1.0 | mmol*mol^-1 | _ | _ | _
end
@columns mutable struct JarvisLinearDeclineVPD{T} <: AbstractJarvisVPD
    d0::T    | 1.0 | Pa | _ | _ | _
end


abstract type AbstractJarvisStomatalConductance <: AbstractStomatalConductance end

"""
Jarvis stomatal conductance model

Combines factors from soilmethod, co2 method, vpdmethod and tempmethod
to gain an overall stomatal conductance.
"""
@columns struct JarvisStomatalConductance{JC<:AbstractJarvisCO2,
                               JV<:AbstractJarvisVPD, JL<:AbstractJarvisLight,
                               JT<:AbstractJarvisTemp, MoMeS, mMoMoS
                              } <: AbstractJarvisStomatalConductance
    co2method::JC   | JarvisNonlinearCO2() | _             | _ | _ | _
    vpdmethod::JV   | JarvisLohammerVPD()  | _             | _ | _ | _
    lightmethod::JL | JarvisLight()        | _             | _ | _ | _
    tempmethod::JT  | JarvisTemp2()        | _             | _ | _ | _
    gsmin::MoMeS    | 1.0                  | mol*m^-2*s^-1 | _ | _ | _
    gsref::MoMeS    | 1.0                  | mol*m^-2*s^-1 | _ | _ | _
    vmfd::mMoMoS    | 1.0                  | mmol*mol^-1   | _ | _ | _
end

@MixinFvCBVars mutable struct JarvisVars{mMoMo}
    vmleaf::mMoMo  | 1.0 | mmol*mol^-1 | _
end


photo_init!(f::AbstractJarvisStomatalConductance, v) = v.vmleaf = f.vmfd
photo_update!(::AbstractJarvisStomatalConductance, v, tleaf1) = v.vmleaf = v.vpdleaf / v.pressure

update_extremes!(f::AbstractJarvisStomatalConductance, v) = begin
    v.aleaf = -v.rd
    v.gs = f.gsmin
end


"""
    stomatal_conductance!(f::AbstractJarvisStomatalConductance, v, p)
Stomatal conductance and for the Jarvis model
"""
function stomatal_conductance!(f::AbstractJarvisStomatalConductance, v)
    v.gs = factor_conductance(f, v)
    v.ac = rubisco_limited_rate(f, v)
    v.aj = transport_limited_rate(f, v)
    v.aleaf = min(v.ac, v.aj) - v.rd
    nothing
end


"""
    rubisco_limited_rate(v, p)
Solution when Rubisco activity is limiting for the Jarvis model """
@inline function rubisco_limited_rate(f::AbstractJarvisStomatalConductance, v)
    a = 1.0 / v.gs
    b = (v.rd - v.vcmax) / v.gs - v.cs - v.km
    c = v.vcmax * (v.cs - v.gammastar) - v.rd * (v.cs + v.km)
    quad(Lower(), a, b, c)
end

"""
    transport_limited_rate(f::AbstractJarvisStomatalConductance, v, p)
Solution when electron transport rate is limiting for the Jarvis model
"""
@inline function transport_limited_rate(f::AbstractJarvisStomatalConductance, v)
    a = 1.0 / v.gs
    b = (v.rd - v.vj) / v.gs - v.cs - 2v.gammastar
    c = v.vj * (v.cs - v.gammastar) - v.rd * (v.cs + 2v.gammastar)
    quad(Lower(), a, b, c)
end

"""
Calculate stomatal conductance gs according to the Jarvis model.

This model calculates gs by multiplying together factors for several
environmental variables: light, VPD, CO2 and temperature.
"""
function factor_conductance(f, v)
    flight = min(max(calc_flight(f.lightmethod, v, f), 1.0), 0.0)
    fvpd = min(max(calc_fvpd(f.vpdmethod, v, f), 1.0), 0.0)
    fco2 = min(max(calc_fco2(f.co2method, v, f), 1.0), 0.0)
    ftemp = min(max(calc_ftemp(f.tempmethod, v, f), 1.0), 0.0)

    return (f.gsref - f.gsmin) * flight * fvpd * fco2 * ftemp * v.fsoil + f.gsmin
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

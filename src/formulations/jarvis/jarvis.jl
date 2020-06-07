# Jarvis model
# Uses factor combination to detaermine stomatal conductance.

abstract type AbstractJarvisLight end

@columns mutable struct JarvisLight{T} <: AbstractJarvisLight
    i0::T   | 1.0 | mol*m^-2*s^-1 | _ | _
end

abstract type AbstractJarvisCO2 end

struct JarvisNoCO2 <: AbstractJarvisCO2 end

@columns mutable struct JarvisLinearCO2{T} <: AbstractJarvisCO2
    gsja::T | 1.0 | μmol^-1*mol | _ | _
end
@columns mutable struct JarvisNonlinearCO2{T} <: AbstractJarvisCO2
    gsjb::T | 1.0 | μmol*mol^-1 | _ | _
end


abstract type AbstractJarvisTemp end

@mix @columns struct JarTemp{T}
    tmax::T | (273.15 + 40.0) | K | _ | _
    tref::T | (273.15 + 25.0) | K | _ | _
    t0::T   | 273.15          | K | _ | _
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
    vk1::T | 1.0 | _ | _ | _
    vk2::T | 1.0 | _ | _ | _
end

"""
Non-linear Lohammer response to vapour pressure deficit.
Parameters vpd1 and vpd2 are in Pascals.
"""
@columns mutable struct JarvisLohammerVPD{T} <: AbstractJarvisVPD
    vpd1::T  | 1.0 | Pa | _ | _
    vpd2::T  | 1.0 | Pa | _ | _
end
@columns mutable struct JarvisFractionDeficitVPD{T} <: AbstractJarvisVPD
    vmfd0::T | 1.0 | mmol*mol^-1 | _ | _
end
@columns mutable struct JarvisLinearDeclineVPD{T} <: AbstractJarvisVPD
    d0::T    | 1.0 | Pa | _ | _
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
    co2method::JC   | JarvisNonlinearCO2() | _             | _ | _
    vpdmethod::JV   | JarvisLohammerVPD()  | _             | _ | _
    lightmethod::JL | JarvisLight()        | _             | _ | _
    tempmethod::JT  | JarvisTemp2()        | _             | _ | _
    gsmin::MoMeS    | 1.0                  | mol*m^-2*s^-1 | _ | _
    gsref::MoMeS    | 1.0                  | mol*m^-2*s^-1 | _ | _
    vmfd::mMoMoS    | 1.0                  | mmol*mol^-1   | _ | _
end

@MixinFvCBVars mutable struct JarvisVars{mMoMo}
    vmleaf::mMoMo  | 1.0 | mmol*mol^-1 | _
end


photo_init!(v, m::AbstractJarvisStomatalConductance) = v.vmleaf = m.vmfd
photo_update!(v, m::AbstractJarvisStomatalConductance, tleaf1) = v.vmleaf = v.vpdleaf / v.pressure

update_extremes!(v, m::AbstractJarvisStomatalConductance) = begin
    v.aleaf = -v.rd
    v.gs = m.gsmin
end


"""
    stomatal_conductance!(v, m::AbstractJarvisStomatalConductance, p)

Stomatal conductance and for the Jarvis model
"""
function stomatal_conductance!(v, m::AbstractJarvisStomatalConductance)
    v.gs = factor_conductance(m, v)
    v.ac = rubisco_limited_rate(m, v)
    v.aj = transport_limited_rate(m, v)
    v.aleaf = min(v.ac, v.aj) - v.rd
    v.aleaf, v.gs
end


"""
    rubisco_limited_rate(v, p)
Solution when Rubisco activity is limiting for the Jarvis model """
@inline function rubisco_limited_rate(m::AbstractJarvisStomatalConductance, v)
    a = 1.0 / v.gs
    b = (v.rd - v.vcmax) / v.gs - v.cs - v.km
    c = v.vcmax * (v.cs - v.gammastar) - v.rd * (v.cs + v.km)
    quad(Lower(), a, b, c)
end

"""
    transport_limited_rate(m::AbstractJarvisStomatalConductance, v, p)
Solution when electron transport rate is limiting for the Jarvis model
"""
@inline function transport_limited_rate(m::AbstractJarvisStomatalConductance, v)
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
function factor_conductance(m, v)
    flight = min(max(calc_flight(m.lightmethod, v), 1.0), 0.0)
    fvpd = min(max(calc_fvpd(m.vpdmethod, v), 1.0), 0.0)
    fco2 = min(max(calc_fco2(m.co2method, v), 1.0), 0.0)
    ftemp = min(max(calc_ftemp(m.tempmethod, v), 1.0), 0.0)

    return (m.gsref - m.gsmin) * flight * fvpd * fco2 * ftemp * v.fsoil + m.gsmin
end

""" Response to incident radiation in umol m^-2 s^-1 """
calc_flight(m::JarvisLight, v) = m.i0 > zero(m.i0) ? v.par / (v.par + m.i0) : 1.0


""" Hyperbolic decline with VPD (VPD in Pa) * BRAY (ALEX BOSC) """
calc_fvpd(m::JarvisHyperbolicVPD, v) =
    v.vpdleaf > zero(v.vpdleaf) ? 1.0 / (m.vk1 * v.vpdleaf^m.vk2) : 0.0
""" Lohammer response to VPD (VPDD0 in Pa) * ECOCRAFT """
calc_fvpd(m::JarvisLohammerVPD, v) =
    v.vpdleaf >= m.vpd1 ? 1.0 - (v.vpdleaf - m.vpd1) / (m.vpd2 - m.vpd1) : 1.0
""" Mole fraction deficit (VMFD in mmol mol - 1) * MARK RAYMENT """
calc_fvpd(m::JarvisFractionDeficitVPD, v) =
    m.vmfd0 > 0.0 ? 1.0 - v.vmleaf / m.vmfd0 : 1.0
""" Linear decline with VPD (VPD1, VPD2 in Pa) * GROMIT (TIM RANDLE) """
calc_fvpd(m::JarvisLinearDeclineVPD, v) = m.d0 > 0.0 ? 1.0 / (1.0 + v.vpdleaf / m.d0) : 1.0


""" No response to temperature in Jarvis stomatal conductance """
calc_ftemp(m::JarvisNoTemp, v) = 1.0
calc_ftemp(m::JarvisTemp1, v) = begin
    p = (m.tmax - m.tref) / (m.tref - m.t0)
    (v.tleaf - m.t0) * ((m.tmax - v.tleaf)^p) / ((m.tref - m.t0) * ((m.tmax - tref)^p))
end
calc_ftemp(m::JarvisTemp2, v) =
    (v.tleaf - m.t0) * (2 * m.tmax - m.t0 - v.tleaf) / ((m.tref - m.t0) * (2 * m.tmax - m.t0 - m.tref))


""" No influence from CO2 for Jarvis stomatal conductance"""
calc_fco2(m::JarvisNoCO2, v) = 1.0
""" Linear response to CO2 for Jarvis stomatal conductance"""
calc_fco2(m::JarvisLinearCO2, v) =
    m.gsja != 0.0 ? 1 - m.gsja * (v.cs - 350.0) : 1.0
""" Non-linear response to CO2 for Jarvis stomatal conductance"""
calc_fco2(m::JarvisNonlinearCO2, v) =
    m.gsjb != 0.0 ? (m.gsjb + 350.0) / (m.gsjb + v.cs) : 1.0

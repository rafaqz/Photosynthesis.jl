
abstract type AbstractJarvisStomatalConductance <: AbstractStomatalConductance end

"""
    JarvisStomatalConductance(co2method, vpdmethod, lightmethod, tempmethod, g0, gsref, vmfd)

Jarvis stomatal conductance model.

Combines factors from soilmethod, co2 method, vpdmethod and tempmethod
to gain an overall stomatal conductance.

$(FIELDDOCTABLE)
"""
@columns struct JarvisStomatalConductance{JC,JV,JL,JT,G0,GS,V} <: AbstractJarvisStomatalConductance
    co2method::JC   | JarvisNonlinearCO2() | _             | _ | _
    vpdmethod::JV   | JarvisLohammerVPD()  | _             | _ | _
    lightmethod::JL | JarvisLight()        | _             | _ | _
    tempmethod::JT  | JarvisTemp2()        | _             | _ | _
    g0::G0          | 1.0                  | mol*m^-2*s^-1 | _ | _
    gsref::GS       | 1.0                  | mol*m^-2*s^-1 | _ | _
    vmfd::V         | 1.0                  | mmol*mol^-1   | _ | _
end

g0(f::JarvisStomatalConductance) = f.g0

@MixinFvCBVars mutable struct JarvisVars{mMoMo}
    vmleaf::mMoMo  | 1.0 | mmol*mol^-1 | _
end

gs_init!(v, m::AbstractJarvisStomatalConductance) = v.vmleaf = m.vmfd
gs_update!(v, m::AbstractJarvisStomatalConductance, tleaf1) = v.vmleaf = v.vpdleaf / v.pressure

stomatal_conductance!(v, m::AbstractJarvisStomatalConductance) = begin
    v.gs = factor_conductance(m, v)
    v.ac = rubisco_limited_rate(m, v)
    v.aj = transport_limited_rate(m, v)
    v.aleaf = min(v.ac, v.aj) - v.rd
    v.aleaf, v.gs
end

@inline function rubisco_limited_rate(m::AbstractJarvisStomatalConductance, v)
    a = 1.0 / v.gs
    b = (v.rd - v.vcmax) / v.gs - v.cs - v.km
    c = v.vcmax * (v.cs - v.gammastar) - v.rd * (v.cs + v.km)
    quad(Lower(), a, b, c)
end

@inline function transport_limited_rate(m::AbstractJarvisStomatalConductance, v)
    a = 1.0 / v.gs
    b = (v.rd - v.vj) / v.gs - v.cs - 2v.gammastar
    c = v.vj * (v.cs - v.gammastar) - v.rd * (v.cs + 2v.gammastar)
    quad(Lower(), a, b, c)
end

"""
    factor_conductance(m, v)

Calculate stomatal conductance `gs` according to the
Jarvis model `m`, given variables `v`.

This model calculates `gs` by multiplying together factors for several
environmental variables: light, VPD, CO2 and temperature.
"""
function factor_conductance(m, v)
    flight = min(max(light_factor(m.lightmethod, v), 1.0), 0.0)
    fvpd = min(max(vpd_factor(m.vpdmethod, v), 1.0), 0.0)
    fco2 = min(max(co2_factor(m.co2method, v), 1.0), 0.0)
    ftemp = min(max(temp_factor(m.tempmethod, v), 1.0), 0.0)

    return (m.gsref - m.g0) * flight * fvpd * fco2 * ftemp * v.fsoil + m.g0
end



"""
Jarvis light response factor formulations
Run in [`light_factor`](@ref) methods.
"""
abstract type AbstractJarvisLight end

"""
    light_factor(m::AbstractJarvisLight, v)

Calculates the light related conductance factor.

Returns a value between 0.0 and 1.0
"""
function light_factor end

"""
    JarvisLight(i0)i0

Factor response to incident radiation

$(FIELDDOCTABLE)
"""
@columns mutable struct JarvisLight{T} <: AbstractJarvisLight
    i0::T   | 1.0 | mol*m^-2*s^-1 | _ | _
end

light_factor(m::JarvisLight, v) = m.i0 > zero(m.i0) ? v.par / (v.par + m.i0) : oneunit(m.i0)


"""
Jarvis CO2 response factor formulations.
Run in [`co2_factor`](@ref) methods.
"""
abstract type AbstractJarvisCO2 end

"""
    co2_factor(m::AbstractJarvisLight, v)

Calculates the CO2 related conductance factor.

Returns a value between 0.0 and 1.0
"""
function co2_factor end

"""
    JarvisNoCO2()

No influence from CO2 for Jarvis stomatal conductance
"""
struct JarvisNoCO2 <: AbstractJarvisCO2 end

co2_factor(m::JarvisNoCO2, v) = 1.0

"""
    JarvisLinearDeclineVPD(gsja)

Linear response to CO2 for Jarvis stomatal conductance

$(FIELDDOCTABLE)
"""
@columns mutable struct JarvisLinearCO2{T} <: AbstractJarvisCO2
    gsja::T | 1.0 | μmol^-1*mol | _ | _
end

co2_factor(m::JarvisLinearCO2, v) =
    m.gsja != 0.0 ? 1 - m.gsja * (v.cs - 350.0) : 1.0

"""
    JarvisNonlinearCO2(gsjb)

Non-linear response to CO2 for Jarvis stomatal conductance

$(FIELDDOCTABLE)
"""
@columns mutable struct JarvisNonlinearCO2{T} <: AbstractJarvisCO2
    gsjb::T | 1.0 | μmol*mol^-1 | _ | _
end

co2_factor(m::JarvisNonlinearCO2, v) =
    m.gsjb != 0.0 ? (m.gsjb + 350.0) / (m.gsjb + v.cs) : 1.0



"""
Temperature response factor formulations
Run in [`temp_factor`](@ref) methods.
"""
abstract type AbstractJarvisTemp end

@mix @columns struct JarTemp{T}
    tmax::T | (273.15 + 40.0) | K | _ | _
    tref::T | (273.15 + 25.0) | K | _ | _
    t0::T   | 273.15          | K | _ | _
end

"""
    temp_factor(m::AbstractJarvisLight, v)

Calculates the temperature related conductance factor.

Returns a value between 0.0 and 1.0
"""
function temp_factor end

"""
    JarvisNoTemp()

No response to temperature in Jarvis stomatal conductance
"""
struct JarvisNoTemp <: AbstractJarvisTemp end

temp_factor(m::JarvisNoTemp, v) = 1.0

"""
    JarvisTemp1(tmax, tref, t0)

$(FIELDDOCTABLE)
"""
@JarTemp mutable struct JarvisTemp1{} <: AbstractJarvisTemp end

temp_factor(m::JarvisTemp1, v) = begin
    p = (m.tmax - m.tref) / (m.tref - m.t0)
    (v.tleaf - m.t0) * ((m.tmax - v.tleaf)^p) / ((m.tref - m.t0) * ((m.tmax - tref)^p))
end

"""
    JarvisTemp2(tmax, tref, t0)

$(FIELDDOCTABLE)
"""
@JarTemp mutable struct JarvisTemp2{} <: AbstractJarvisTemp end

temp_factor(m::JarvisTemp2, v) =
    (v.tleaf - m.t0) * (2 * m.tmax - m.t0 - v.tleaf) / 
    ((m.tref - m.t0) * (2 * m.tmax - m.t0 - m.tref))



"""
Vapour pressure deficit response-factor formulations
Run in [`vpd_factor`](@ref) methods.
"""
abstract type AbstractJarvisVPD end

"""
    vpd_factor(m::AbstractJarvisLight, v)

Calculates the vapour-pressure-deficit related conductance factor.

Returns a value between 0.0 and 1.0
"""
function vpd_factor end

"""
    JarvisHyperbolicVPD(vk1, vk2)

Hyperbolic decline with VPD. BRAY (ALEX BOSC)
Parameters vk1 and vk2 are the dimensionless scalar and exponent.

$(FIELDDOCTABLE)
"""
@columns mutable struct JarvisHyperbolicVPD{T} <: AbstractJarvisVPD
    vk1::T | 1.0 | _ | _ | _
    vk2::T | 1.0 | _ | _ | _
end

vpd_factor(m::JarvisHyperbolicVPD, v) =
    v.vpdleaf > zero(v.vpdleaf) ? 1.0 / (m.vk1 * v.vpdleaf^m.vk2) : 0.0

"""
    JarvisLohammerVPD(vpd1, vpd2)

Non-linear Lohammer response to vapour pressure deficit.
Parameters vpd1 and vpd2 are in Pascals.

$(FIELDDOCTABLE)
"""
@columns mutable struct JarvisLohammerVPD{T} <: AbstractJarvisVPD
    vpd1::T  | 1.0 | Pa | _ | _
    vpd2::T  | 1.0 | Pa | _ | _
end

vpd_factor(m::JarvisLohammerVPD, v) =
    v.vpdleaf >= m.vpd1 ? 1.0 - (v.vpdleaf - m.vpd1) / (m.vpd2 - m.vpd1) : 1.0

"""
    JarvisFractionDeficitVPD(vmdf0)

Mole fraction deficit MARK RAYMENT

$(FIELDDOCTABLE)
"""
@columns mutable struct JarvisFractionDeficitVPD{T} <: AbstractJarvisVPD
    vmfd0::T | 1.0 | mmol*mol^-1 | _ | _
end

vpd_factor(m::JarvisFractionDeficitVPD, v) =
    m.vmfd0 > 0.0 ? 1.0 - v.vmleaf / m.vmfd0 : 1.0

"""
    JarvisLinearDeclineVPD(d0)

Linear decline with VPD (VPD1, VPD2 in Pa) * GROMIT (TIM RANDLE)

$(FIELDDOCTABLE)
"""
@columns mutable struct JarvisLinearDeclineVPD{T} <: AbstractJarvisVPD
    d0::T    | 1.0 | Pa | _ | _
end

vpd_factor(m::JarvisLinearDeclineVPD, v) =
    m.d0 > 0.0 ? 1.0 / (1.0 + v.vpdleaf / m.d0) : 1.0


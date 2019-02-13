"Physiological constants for CO2 and rubisco compensation points"
abstract type AbstractCompensation end

@columns struct BadgerCollatzCompensation{μMoMo,kJMo,K} <: AbstractCompensation
    Kc25::μMoMo     | 404.0    | μmol*mol^-1    | Gamma(100, 404/100)    | (0.0, 1000.0) | "MM coefft of Rubisco for CO2"
    Ko25::μMoMo     | 248000.0 | μmol*mol^-1    | Gamma(100, 248000/100) | (0.0, 10000000.0) | "MM coefft of Rubisco for O2"
    ΔHa_Kc::kJMo    | 59.4     | kJ*mol^-1      | Gamma(100, 59.4/100)   | (0.0, 200.0) | "Temp. response of Kc"
    ΔHa_Ko::kJMo    | 36.0     | kJ*mol^-1      | Gamma(100, 36/100)     | (0.0, 200.0) | "Temp. response of Ko"
    tref::K         | 298.15   | K              | _                      | (250.0, 350.0) | "Temperature reference, usually 25.0° C"
end

"""
Bernacchi et al 2001 PCE 24: 253 - 260
Extra deactivation terms may be required above 40°C.
"""
@columns struct BernacchiCompensation{μMoMo,kJMo,K} <: AbstractCompensation
    Kc25::μMoMo     | 404.9    | μmol*mol^-1    | Gamma(10, 404.9/10)    | (1.0, 1000.0) | "MM coefft of Rubisco for CO2"
    Ko25::μMoMo     | 278400.0 | μmol*mol^-1    | Gamma(10, 278400/10)   | (1.0, 10000000.0) | "MM coefft of Rubisco for O2"
    Γ☆25::μMoMo     | 42.75    | μmol*mol^-1    | Gamma(10, 42.75/10)    | (1.0, 100.0)   | _
    ΔHa_Kc::kJMo    | 79.43    | kJ*mol^-1      | Gamma(10, 79.43/10)    | (1.0, 100.0)   | "Temp. response of Kc"
    ΔHa_Ko::kJMo    | 36.38    | kJ*mol^-1      | Gamma(10, 36.38/10)    | (1.0, 100.0)   | "Temp. response of Ko"
    ΔHa_Γ☆::kJMo    | 37.83    | kJ*mol^-1      | Gamma(10, 37.83/10)    | (1.0, 100.0)   | _
    tref::K         | 298.15   | K              | _                      | (250.0, 350.0) | _
end

""" Calculates Γ*, or the CO2 compensation point in the absence of
non-photorespiratory respiration """
function co2_compensation_point end

""" 
    co2_compensation_point(::BadgerCollatzCompensation, v) 
"""
co2_compensation_point(f::BadgerCollatzCompensation, v) = begin
    # if tleaf < -1.0  calculate gamma for t = -1 (quadratic not applicable)
    tleaf = max(K(-1.0°C), v.tleaf)
    # TODO these numbers should be parameters?
    36.9μmol*mol^-1 + 1.88μmol*mol^-1*K^-1 *
    (tleaf - f.tref) + 0.036μmol*mol^-1*K^-2 * (tleaf - f.tref)^2
end

""" 
    co2_compensation_point(::BernacchiCompensation, v)
Bernacchi et al 2001 PCE 24: 253-260 
"""
co2_compensation_point(f::BernacchiCompensation, v) =
    arrhenius(f.Γ☆25, f.ΔHa_Γ☆, v.tleaf, f.tref)

""" 
    rubisco_compensation(::BadgerCollatzCompensation, v)
Calculates Km, or the effective Michaelis-Menten coefficient of Rubisco activity.
"""
rubisco_compensation_point(f::AbstractCompensation, v) = begin
    Kc = arrhenius(f.Kc25, f.ΔHa_Kc, v.tleaf, f.tref)
    Ko = arrhenius(f.Ko25, f.ΔHa_Ko, v.tleaf, f.tref)
    Kc * (1 + OI / Ko)
end

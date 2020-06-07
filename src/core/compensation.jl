"Physiological constants for CO2 and rubisco compensation points"
abstract type Compensation end

""" 
    co2_compensation_point(formulation::Compensation, vars)

Calculates Γ*, or the CO2 compensation point in the absence of
non-photorespiratory respiration 
"""
function co2_compensation_point end


@flattenable @columns struct BadgerCollatzCompensation{μMoMo,kJMo,K} <: Compensation
#   Field        | Flat   | Default  | Units          | Bounds            | Description
    Kc25::μMoMo  | false  | 404.0    | μmol*mol^-1    | (0.0, 1000.0)     | "MM coefft of Rubisco for CO2"
    Ko25::μMoMo  | false  | 248000.0 | μmol*mol^-1    | (0.0, 10000000.0) | "MM coefft of Rubisco for O2"
    ΔHa_Kc::kJMo | false  | 59.4     | kJ*mol^-1      | (0.0, 200.0)      | "Temp. response of Kc"
    ΔHa_Ko::kJMo | false  | 36.0     | kJ*mol^-1      | (0.0, 200.0)      | "Temp. response of Ko"
    tref::K      | false  | 298.15   | K              | (250.0, 350.0)    | "Temperature reference, usually 25.0° C"
end

""" 
    co2_compensation_point(::BadgerCollatzCompensation, v) 
"""
co2_compensation_point(f::BadgerCollatzCompensation, tleaf) = begin
    # if tleaf < -1.0  calculate gamma for t = -1 (quadratic not applicable)
    tleaf = max(K(-1.0°C), tleaf)
    # TODO these numbers should be parameters?
    36.9μmol*mol^-1 + 1.88μmol*mol^-1*K^-1 *
    (tleaf - f.tref) + 0.036μmol*mol^-1*K^-2 * (tleaf - f.tref)^2
end


"""
    BernacchiCompensation(Kc25, Ko25, Γ☆25, ΔHa_Kc, ΔHa_Ko, ΔHa_Γ☆, tref)

Bernacchi et al 2001 PCE 24: 253 - 260. Extra deactivation terms may be required above 40°C.
"""
@flattenable @columns struct BernacchiCompensation{μMoMo,kJMo,K} <: Compensation
#   Field        | Flat  | Default  | Units          | Bounds            | Description
    Kc25::μMoMo  | false | 404.9    | μmol*mol^-1    | (1.0, 1000.0)     | "MM coefft of Rubisco for CO2"
    Ko25::μMoMo  | false | 278400.0 | μmol*mol^-1    | (1.0, 10000000.0) | "MM coefft of Rubisco for O2"
    Γ☆25::μMoMo  | false | 42.75    | μmol*mol^-1    | (1.0, 100.0)      | _
    ΔHa_Kc::kJMo | false | 79.43    | kJ*mol^-1      | (1.0, 100.0)      | "Temp. response of Kc"
    ΔHa_Ko::kJMo | false | 36.38    | kJ*mol^-1      | (1.0, 100.0)      | "Temp. response of Ko"
    ΔHa_Γ☆::kJMo | false | 37.83    | kJ*mol^-1      | (1.0, 100.0)      | _
    tref::K      | false | 298.15   | K              | (250.0, 350.0)    | _
end

""" 
    co2_compensation_point(::BernacchiCompensation, v)
Bernacchi et al 2001 PCE 24: 253-260 
"""
co2_compensation_point(f::BernacchiCompensation, tleaf) =
    arrhenius(f.Γ☆25, f.ΔHa_Γ☆, tleaf, f.tref)

""" 
    rubisco_compensation(::Compensation, v)

Calculates Km, or the effective Michaelis-Menten coefficient of Rubisco activity.
"""
rubisco_compensation_point(f::Compensation, tleaf) = begin
    Kc = arrhenius(f.Kc25, f.ΔHa_Kc, tleaf, f.tref)
    Ko = arrhenius(f.Ko25, f.ΔHa_Ko, tleaf, f.tref)
    Kc * (1 + OI / Ko)
end

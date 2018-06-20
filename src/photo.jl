"""
    photosynthesis!(v, p::FvCBPhoto)
Calculates photosynthesis according to the ECOCRAFT
agreed formulation of the Farquharvon Caemmerer (1982) equations.

Farquhar, G.D., S. Caemmerer and J.A. Berry. 1980. 
A biochemical model of photosynthetic CO2 assimilation in leaves of C3 species. 
Planta. 149:78-90. 
"""
function photosynthesis!(v, p::FvCBPhoto)
    v.gammastar = co2_compensation_point(p.compensation, v, p) # CO2 compensation point, umol mol-1
    v.km = rubisco_compensation_point(p.compensation, v, p) # Michaelis-Menten for Rubisco, umol mol-1
    v.jmax = max_electron_transport_rate(p.vcjmax, v, p) # Potential electron transport rate, umol m-2 s-1
    v.vcmax = max_rubisco_activity(p.vcjmax, v, p) # Maximum Rubisco activity, umol m-2 s-1
    v.rd = respiration(p.respiration, v, p) # Day leaf respiration, umol m-2 s-1

    # soil_water_conductance!(v, p.model.soilmethod, p)

    v.vj = rubisco_regeneration(p.rubisco_regen, v, p)
    extremes!(v, p) && return nothing
    stomatal_conductance!(v, p.model, p)

    nothing
end


""" Calculates Γ*, or the CO2 compensation point in the absence of
non-photorespiratory respiration """
function co2_compensation_point end

""" 
    co2_compensation_point(::BadgerCollatzCompensation, v, p) 
"""
co2_compensation_point(f::BadgerCollatzCompensation, v, p) = begin
    # if tleaf < -1.0  calculate gamma for t = -1 (quadratic not applicable)
    tleaf1 = max(-1.0u"°C", v.tleaf)
    # TODO these numbers should be parameters
    36.9u"μmol*mol^-1" + 1.88u"μmol*mol^-1*°C^-1" *
    (tleaf1 - f.tref) + 0.036u"μmol*mol^-1*°C^-2" * (tleaf1 - f.tref)^2
end

""" 
    co2_compensation_point(::BernacchiCompensation, v, p)
Bernacchi et al 2001 PCE 24: 253-260 
"""
co2_compensation_point(f::BernacchiCompensation, v, p) =
    arrhenius(f.Γ☆25, f.ΔHa_Γ☆, v.tleaf, f.tref)

""" 
    rubisco_compensation(::BadgerCollatzCompensation, v, p)
Calculates Km, or the effective Michaelis-Menten coefficient of Rubisco activity.
"""
rubisco_compensation_point(f::AbstractCompensation, v, p) = begin
    Kc = arrhenius(f.Kc25, f.ΔHa_Kc, v.tleaf, f.tref)
    Ko = arrhenius(f.Ko25, f.ΔHa_Ko, v.tleaf, f.tref)
    Kc * (1 + OI / Ko)
end

""" 
    max_electron_transport_rate(f::DukeVcJmax, v, p)
Wrapper to `max_electron_transport_rate()` allowing Vcmax to be forced linearly to zero at low T """
max_electron_transport_rate(f::DukeVcJmax, v, p) = begin
    v.tleaf < f.tvjdn && return zero(f.jmax25)
    jmax = max_electron_transport_rate(f.jmaxformulation, v, p)
    v.tleaf < f.tvjup && (v.tleaf - f.tvjdn) / (f.tvjup - f.tvjdn) * jmax
    jmax
end

""" 
    max_electron_transport_rate(f::VcJmax, v, p)
Wrapper for `max_electron_transport_rate() that simply runs `max_electron_transport_rate` for `jmaxformulation`
"""
max_electron_transport_rate(f::VcJmax, v, p) = max_electron_transport_rate(f.jmaxformulation, v, p)

""" 
    max_electron_transport_rate(f::Jmax, v, p)
Calculates the potential max_electron transport rate (Jmax) at the leaf temperature 
"""
max_electron_transport_rate(f::Jmax, v, p) = begin
    tleafK = v.tleaf |> u"K"
    K25 = 25.0u"°C" |> u"K"
    f.jmax25 * e^((tleafK - K25) * f.eavj / (u"R" * tleafK * K25)) *
    (1 + e^((f.delsj * K25 - f.edvj) / (u"R" * K25))) /
    (1 + e^((f.delsj * tleafK - f.edvj) / (u"R" * tleafK)))
end

"""
Calculates the maximum Rubisco activity (Vcmax) at the leaf temperature.

There is still disagreement as to whether this function has an optimum or not. 
Both versions are well-behaved for tleaf < 0.0
"""
function max_rubisco_activity end

""" Function allowing Vcmax to be forced linearly to zero at low T.  
Introduced for Duke data. """
function max_rubisco_activity(f::DukeVcJmax, v, p)
    v.tleaf < f.tvjdn && return zero(f.vcmax25)
    vcmax = max_rubisco_activity(f.vcmaxformulation, v, p)
    v.tleaf < f.tvjup && (v.tleaf - f.tvjdn) / (f.tvjup - f.tvjdn) * vcmax
    vcmax
end

""" Wrapper with no alterations to vcmax
Simply runs max_rubisco_activity again for `vcmaxformlation`"""
max_rubisco_activity(f::VcJmax, v, p) = max_rubisco_activity(f.vcmaxformulation, v, p)

""" Vcmax forulation with no optimum"""
function max_rubisco_activity(f::NoOptimumVcmax, v, p)
    tleafK = v.tleaf |> u"K"
    K25 = 25u"°C" |> u"K"
    f.vcmax25 * e^((f.eavc * (v.tleaf - 25u"°C")) / (K25 * u"R" * tleafK))
end
""" Vcmax formulation with optimum"""
function max_rubisco_activity(f::OptimumVcmax, v, p)
    tleafK = v.tleaf |> u"K"
    K25 = 25u"°C" |> u"K"
    f.vcmax25 * e^((v.tleaf - 25u"°C") * f.eavc / (u"R" * tleafK * K25)) *
    (1.0 + e^((f.delsc * K25 - f.edvc) / (u"R" * K25))) /
    (1.0 + e^((f.delsc * tleafK - f.edvc) / (u"R" * tleafK)))
end

""" Calculates respiration from temperature using a Q10 (exponential) formulation """
function respiration(f::Respiration, v, p)
    v.tleaf < f.tbelow && return zero(f.rd0 * f.dayresp)
    # Make sure light suppression of dark respiration only occurs when it is light.
    # See Atkin et al. 1998 (or 2001?). From yplantmc
    resplightfrac = v.par < 100oneunit(v.par) ? oneunit(f.dayresp) : f.dayresp

    f.rd0 * e^(f.q10f * (v.tleaf - f.rtemp)/10) * resplightfrac 
end

function respiration(f::AcclimatizedRespiration, v, p)
    v.tleaf < f.tbelow && return zero(f.rd0 * f.dayresp)
    resplightfrac = v.par < 100oneunit(v.par) ? oneunit(f.dayresp) : f.dayresp

    rd0acc = f.rd0 * exp(f.k10f * (f.tmove - f.rtemp))
    rd0acc * exp(f.q10f * (v.tleaf - f.rtemp)) * resplightfrac 
end 

""" 
    rubisco_regeneration(f::RubiscoRegen, v, p)
RuBP-regen rate, in umol m-2 s-1 """
function rubisco_regeneration(f::RubiscoRegen, v, p)
    a = f.theta
    b = -(f.ajq * v.par + v.jmax)
    c = f.ajq * v.par * v.jmax
    j = quad(Val{:lower}, a, b, c) # Actual e- transport rate, umol m-2 s-1
    return j / 4.0
end

function jmaxN()
  jmaxa * leafn + jmaxb
end

function vcmax()
    vcmaxa * leafn + vcmaxb
end

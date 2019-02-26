abstract type AbstractEnergyBalance end

abstract type AbstractBoundaryConductance end

@columns struct BoundaryConductance{M,MolMS} <: AbstractBoundaryConductance
    leafwidth::M | 0.05 | m             | Gamma(2, 0.05/2) | (0.0, 1.0) | "Mean width of leaves"
    gsc::MolMS   | 1.0  | mol*m^-2*s^-1 | Gamma(2, 2)      | (0.0, 1.0) | "??"
end


abstract type AbstractRadiationConductance end

@columns struct YingPingRadiationConductance{T} <: AbstractRadiationConductance
    rdfipt::T | 1.0 | _ | _ | (0.0, 2.0) | _
    tuipt::T  | 1.0 | _ | _ | (0.0, 2.0) | _
    tdipt::T  | 1.0 | _ | _ | (0.0, 2.0) | _
end

abstract type AbstractDecoupling end

struct McNaughtonJarvisDecoupling <: AbstractDecoupling end
struct NoDecoupling <: AbstractDecoupling end


"""
run_enbal!(p, v)
Potosynthesis model runner for all formulations.
"""
function run_enbal! end

function enbal! end

"""
    model_init!(v, f)
Runs any model initialisation that needs to happen at the start of energy balance
"""
function enbal_init!(f, v) end

"""
    model_update!(v, model, tleaf1)
Runs any model specific variable updates that need to happen at the end of
the leaf temperature convergene loop
"""
function enbal_update!(f, v, tleaf) end

"""
    calc_decoupling(f::McNaughtonJarvisDecoupling, v, gbv, gsv)
Calculate decoupling coefficient (McNaughton and Jarvis 1986)
"""
calc_decoupling(f::McNaughtonJarvisDecoupling, v) = begin
    γc = CPAIR * AIRMA * v.pressure / v.lhv
    epsilon = ustrip(v.slope / γc) # TODO why is ustrip needed here?
    (1.0 + epsilon) / (1.0 + epsilon + v.gbv / v.gsv)
end
"""
    calc_decoupling(f::NoDecoupling, v, gbv, gsv)
Don't calculate decoupling
"""
calc_decoupling(f::NoDecoupling, v) = 0.0

""" Calculate vapour pressure change with temperature -
Slope `s` for Penman-Monteith equation, in Pa K^-1
"""
calc_slope(tair) = (saturated_vapour_pressure(tair + 0.1oneunit(tair)) -
                  saturated_vapour_pressure(tair)) / 0.1oneunit(tair)

"""
Boundary layer conductance for heat - single sided, free convection
"""
function boundary_conductance_free(f::AbstractBoundaryConductance, v)
    gb = free_boundary_conductance(DHEAT, v.tleaf, v.tair, f.leafwidth)
    # Convert from m s-1 to mol m-2 s-1
    gb * cmolar(v.pressure, v.tair)
end

"""
Boundary layer conductance for heat - single sided, forced convection
"""
function boundary_conductance_forced(f::AbstractBoundaryConductance, v)
    gb = forced_boundary_conductance(v.windspeed, f.leafwidth)
    # Convert from m s-1 to mol m-2 s-1
    gb * cmolar(v.pressure, v.tair)
end

radiation_conductance(f::YingPingRadiationConductance, v) =
    yingping_radiation_conductance(v.tair, f.rdfipt, f.tuipt, f.tdipt)


"""
    yingping_radiation_conductance(tair, rdfipt, tuipt, tdipt)
Returns the "radiation conductance" at given temperature.
Formula from Ying-Ping"s version of Maestro.
See also Jones (1992) p. 108.0

Return mol*m^-2*s^-1
"""
yingping_radiation_conductance(tair, rdfipt, tuipt, tdipt) =
    4.0 * σ * ((tair |> K)^3.0) * rdfipt / tdipt * EMLEAF * 
    (tdipt + tuipt) / (CPAIR * AIRMA)


"""
    boundary_conductance_forced(Ta, ρ, U, w)
Leaf boundary layer conductance for heat - single sided, forced convection
Ta is air temperature
ρ is pressure
U is wind speed
w is leaf width

See Leuning et al (1995) PCE 18:1183-1200 Eqn E1
"""
forced_boundary_conductance(wind, width) = 0.003m*s^-1 * sqrt(ustrip(wind / width))

"""
    boundary_conductance_free(α, Tl, Ta, w)
Leaf boundary layer conductance for heat  - single sided, free convection 
Dh is the molecular diffusivity to heat
Ta is air temperature
Tl is leaf temperature
w is leaf width

See Leuning et al (1995) PCE 18:1183-1200 Eqn E3
"""
free_boundary_conductance(Dh, Tl, Ta, w) = 0.5Dh * (grashof_number(Tl, Ta, w)^(1/4))/w

"""
    grashof_number(Ts, Ta, d)
Calculates the Grashof number given leaf temperature, air tempereature
and leaf width.

Ts: Object temperature
T: background temperature
a: coefficient of thermal expansion.

return: dimensionless

See Leuning et al (1995) PCE 18:1183-1200 Eqn E4
"""
grashof_number(Tl, Ta, w) = 1.6e8m^-3*K^-1 * w^3 * abs(Tl - Ta)


"""
    latent_heat_water_vapour(tair)

Caculates the late hear of water vapour from the air temperature.
"""
latent_heat_water_vapour(Ta) = (H2OLV0 - 2.365e3J*kg^-1*K^-1 * Ta) * H2OMW


"""
    arrhenius(kt, ea, t, tref)
The Arrhenius function.
kT is the value at tref deg # 
Ea the activation energy (j mol - 1) and 
T the temp (deg #).
Standard form and temperature difference form.
"""
arrhenius(A, Ea, T::typeof(1.0°C)) = arrhenius(A, Ea, T |> K)
arrhenius(A, Ea, T) = A * exp(Ea / (R * T))
arrhenius(kref, Ea, T::typeof(1.0°C), Tref::typeof(1.0°C)) = arrhenius(kref, Ea, T |> K, Tref |> K)
arrhenius(kref, Ea, T::typeof(1.0K), Tref::typeof(1.0K)) = kref * exp(Ea * (T - Tref) / (R * T * Tref))


"""
# Calculate saturated water vapour pressure (Pa) at temperature Ta (Celsius)
# from Jones 1992 p 110 (note error in a - wrong units)
"""
saturated_vapour_pressure(Ta) = 613.75Pa * exp(17.502 * Ta / (K(240.97°C) + Ta))


"""
    penman_monteith(pressure, slope, lhv, rnet, vpd, gh, gv)
This subroutine calculates evapotranspiration by leaves using the Penman - Monteith equation.

Inputs:

ρ atmospheric pressure, Pa
Δ slope of VPD / T curve, Pa K-1
lhv latent heat of water at air T, J mol-1
Rn net radiation, J m-2 s-1
Da vapour pressure deficit of air, Pa
gh boundary layer conductance to heat (freeforcedradiative components), mol m-2 s-1
gv conductance to water vapour (stomatalbdry layer components), mol*m^-2*s^-1

Result in mol*m^-2*s^-1
"""
penman_monteith(ρa, Δ, lhv, Rn, Da, gh, gv) = begin
    #FIXME gv <= zero(gv) && return zero(Rn / lhv)

    γ = CPAIR * ρa * AIRMA / lhv
    ET = (Δ * Rn + CPAIR * gh * Da * AIRMA) / (Δ + γ * gh * (1/gv))

    return ET / lhv
    # if (penmon < 0.0) penmon = 0.0 end # BM 12 / 05 Should not be negative
end

"""
    cmolar(ρ, Ta)
Convert from m*s^-1 to mol*m^-2*s^-1
"""
cmolar(pressure, tair) = pressure / (R * (tair |> K))

"""
run_phototrans!(v, p)
Potosynthesis model runner for all formulations.
"""
run_phototrans!(v, p) = run_phototrans!(v, p, p.photo.model)

"""
    run_phototrans!(v, pm::TuzetModel)
Potosynthesis runner for `TuzetModel` runs phototranspiration inside a root
finding algorithm to calculate leaf water potential
"""
function run_phototrans!(v, p, m::TuzetModel)
    v.psilin = -0.1oneunit(v.psilin)
    v.psil = -0.1oneunit(v.psil)
    bracket = -100.0oneunit(v.psil), zero(v.psil)
    tolz = 1e-03oneunit(v.psil)

    v.psilin = find_zero(x -> leaf_water_potential_finder(x, v, p), bracket, FalsePosition(), atol=tolz)
    nothing
end

"""
    run_phototrans!(v, p, m::AbstractPhotoModel)
Potosynthesis runner for most models. Identical to calling
`phototranspiration` directly
"""
function run_phototrans!(v, p, m::AbstractPhotoModel)
    phototranspiration!(v, p)
    nothing
end

"""
    leaf_water_potential_finder(psilin, v, p)
PSIL finder. A wrapper function that
return the squared difference in PSILIN and PSIL
"""
function leaf_water_potential_finder(psilin, v, p)
    v.ci = zero(v.ci)
    phototranspiration!(v, p)
    v.psilin - v.psil
end

"""
    phototranspiration!(v, p)
This subroutine calculates leaf photosynthesis and transpiration.

These may be calculated by:
(1) assuming leaf temperature = air temperature, cs = ca and ds = da
(2) using iterative scheme of Leuning et al (1995) (PCE 18:1183-1200) to calculate leaf temp, CsCa.

Setting itermax = 0 gives (1); itermax > 0 (suggest 100) gives (2).
"""
function phototranspiration!(v, p)
    model_init!(v, p.photo.model, p)

    # Initialise with environmental values
    v.tleaf = v.tair
    v.vpdleaf = v.vpd
    v.cs = v.ca

    # Calculations that do not depend on tleaf
    v.lhv = latent_heat_water_vapour(v.tair)
    v.slope = calc_slope(v.tair)
    #v.gradn = radiation_conductance(p.radiation_conductance, v, p)
    v.gbhu = boundary_conductance_forced(p.boundary_conductance, v, p)

    converge_tleaf!(v, p) == false && error("leaf temperature convergence failed")

    v.fheat = v.rnet - v.lhv * v.et

    nothing
end

"""
    converge_tleaf!(v, p)
Run phototranpiration process in a loop to converge on leaf temperature
"""
function converge_tleaf!(v, p)
    for iter = 1:p.itermax
        photosynthesis!(v, p.photo)
        v.gbhf = # boundary_conductance_free(p.boundary_conductance, v, p)

        # Total boundary layer conductance for heat
        # See Leuning et al (1995) PCE 18:1183-1200 Eqn E5
        gbh = v.gbhu + v.gbhf
        # Total conductance for heat - two-sided
        v.gh = 2.0(gbh + v.gradn)

        # Total conductance for water vapour
        gbv = GBVGBH * gbh
        gsv = GSVGSC * v.gs # TODO check this is right?? was: * p.gsc
        # gv = nsides * (gbv * gsv) / (gbv + gsv) # already one-sided value
        gv = (gbv * gsv) / (gbv + gsv)

        v.et = penman_monteith(v.pressure, v.slope, v.lhv, v.rnet, v.vpd, v.gh, gv)
        v.decoup = calc_decoupling(p.decoupling, v, p, gbv, gsv)

        # End of subroutine if no iterations wanted.
        p.itermax == 0 || v.aleaf <= zero(v.aleaf) && return true

        gbc = gbh / GBHGBC
        v.cs = v.ca - v.aleaf / gbc # TODO this value is way too low
        tleaf1 = leaftemp(v, p)
        v.vpdleaf = v.et * v.pressure / gv # TODO and this seems too high?
        uconvert(kPa, v.vpdleaf)

        model_update!(v, p.photo.model, p, tleaf1) # Model-specific var updates

        # Check to see whether convergence achieved or failed
        if abs(v.tleaf - tleaf1) < TOL
            v.tleaf = tleaf1
            return true
        end

        v.tleaf = tleaf1 # Update temperature for another iteration
    end
    false
end

leaftemp(v, p) = v.tair + (v.rnet - v.et * v.lhv) / (CPAIR * AIRMA * v.gh)

"""
    model_init!(v, f, p)
Runs any model initialisation that needs to happen at the start of phototranspiration
"""
model_init!(v, f, p) = nothing
model_init!(v, f::AbstractBallBerryModel, p) = v.rhleaf = v.rh
model_init!(v, f::JarvisModel, p) = v.vmleaf = f.vmfd
model_init!(v, f::AbstractMaespaModel, p) = begin
    v.ktot = 10 / (f.totsoilres + 1.0 / f.plantk)
end

"""
    model_update!(v, model, p, tleaf1)
Runs any model specific variable updates that need to happen at the end of
the leaf temperature convergene loop
"""
model_update!(v, model, p, tleaf1) = nothing
model_update!(v, model::JarvisModel, p, tleaf1) = v.vmleaf = v.vpdleaf / v.pressure
model_update!(v, model::AbstractBallBerryModel, p, tleaf1) = begin
    v.rhleaf = 1.0 - v.vpdleaf / saturated_vapour_pressure(tleaf1)
end

"""
    calc_decoupling(f::McNaughtonJarvisDecoupling, v, p, gbv, gsv)
Calculate decoupling coefficient (McNaughton and Jarvis 1986)
"""
calc_decoupling(f::McNaughtonJarvisDecoupling, v, p, gbv, gsv) = begin
    γc = CPAIR * AIRMA * v.pressure / v.lhv
    epsilon = ustrip(v.slope / γc) # TODO why is ustrip needed here?
    (1.0 + epsilon) / (1.0 + epsilon + gbv / gsv)
end
"""
    calc_decoupling(f::NoDecoupling, v, p, gbv, gsv)
Don't calculate decoupling
"""
calc_decoupling(f::NoDecoupling, v, p, gbv, gsv) = 0.0

""" Calculate vapour pressure change with temperature -
slope `s` for Penman-Monteith equation, in Pa K^-1
"""
calc_slope(tair) = (saturated_vapour_pressure(tair + 0.1oneunit(tair)) -
                  saturated_vapour_pressure(tair)) / 0.1oneunit(tair)

"""
Boundary layer conductance for heat - single sided, free convection
"""
function boundary_conductance_free(f::AbstractBoundaryConductance, v, p)
    gb = free_boundary_conductance(DHEAT, v.tleaf, v.tair, f.leafwidth)
    # Convert from m s-1 to mol m-2 s-1
    gb * cmolar(v.pressure, v.tair)
end

"""
Boundary layer conductance for heat - single sided, forced convection
"""
function boundary_conductance_forced(f::AbstractBoundaryConductance, v, p)
    gb = forced_boundary_conductance(v.windspeed, f.leafwidth)
    # Convert from m s-1 to mol m-2 s-1
    gb * cmolar(v.pressure, v.tair)
end

radiation_conductance(f::YingPingRadiationConductance, v, p) =
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
forced_boundary_conductance(wind, width) = 
    0.003m*s^-1 * sqrt(ustrip(wind / width))

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
grashof_number(Tl, Ta, w) = 1.6e8m^-3*°C^-1 * w^3.0 * abs(Tl - Ta)


"""
    latent_heat_water_vapour(tair)

Caculates the late hear of water vapour from the air temperature.
"""
latent_heat_water_vapour(Ta) = (H2OLV0 - 2.365e3J*kg^-1*°C^-1 * Ta) * H2OMW


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
saturated_vapour_pressure(Ta) = 613.75Pa * exp(17.502 * Ta / (240.97°C + Ta))


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

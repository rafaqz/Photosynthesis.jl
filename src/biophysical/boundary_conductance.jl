"""
Boundary conductance models.

Provide parameters for [`vapour_conductance!`](@ref), [`boundary_conductance_free`](@ref), 
and [`boundary_conductance_forced`](@ref).
"""
abstract type AbstractBoundaryConductance end

"""
    boundary_conductance_free(f::AbstractBoundaryConductance, v)

Boundary layer conductance for heat - single sided, free convection,
given variables `v`.

Retuns conductance in `u"mol*m-2*s-1"`.
"""
function boundary_conductance_free end

"""
    boundary_conductance_forced(f::AbstractBoundaryConductance, v)

Boundary layer conductance for heat - single sided, forced convection, 
given variables `v`.

Retuns conductance in `u"mol*m-2*s-1"`.
"""
function boundary_conductance_forced end

"""
    vapour_conductance!(v, f::AbstractBoundaryConductance)

Vapour conductance.

# TODO remove assingments to v
"""
function vapour_conductance! end

"""
    BoundaryConductance(leafwidth, gsc)
    BoundaryConductance(; leafwidth=0.05u"m", gsc=1.0u"mol*m^-2*s^-1")

Standard boundary conductance formulation parameters.
"""
@columns struct BoundaryConductance{M,MolMS} <: AbstractBoundaryConductance
    leafwidth::M | 0.05 | m             | (0.0, 1.0) | "Mean width of leaves"
    gsc::MolMS   | 1.0  | mol*m^-2*s^-1 | (0.0, 1.0) | "Stomatal conductance of the boundary layer to CO₂" # Check this
end

vapour_conductance!(v, f::AbstractBoundaryConductance) = begin
    v.gbhf = boundary_conductance_free(f, v)

    # Total boundary layer conductance for heat
    # See Leuning et al (1995) PCE 18:1183-1200 Eqn E5
    v.gbh = v.gbhu + v.gbhf
    # Total conductance for heat: two-sided
    v.gh = 2.0(v.gbh + v.gradn)

    # Total conductance for water vapour
    v.gbv = GBVGBH * v.gbh
    v.gsv = GSVGSC * f.gsc
    # gv = nsides * (gbv * gsv) / (gbv + gsv) # already one-sided value
    (v.gbv * v.gsv) / (v.gbv + v.gsv)
end

boundary_conductance_free(f::BoundaryConductance, v) = begin
    gb = free_boundary_conductance(v.tleaf, v.tair, f.leafwidth)
    # Convert from m s-1 to mol m-2 s-1
    gb * cmolar(v.pressure, v.tair)
end

boundary_conductance_forced(f::BoundaryConductance, v) = begin
    gb = forced_boundary_conductance(v.windspeed, f.leafwidth)
    # Convert from m s-1 to mol m-2 s-1
    gb * cmolar(v.pressure, v.tair)
end


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

See Leuning et al (1995) PCE 18:1183-1200 Eqn E3
"""
free_boundary_conductance(tleaf, tair, leafwidth) = 
    0.5DHEAT * (grashof_number(tleaf, tair, leafwidth)^(1/4)) / leafwidth


"""
    cmolar(ρ, Ta)

Convert from m*s^-1 to mol*m^-2*s^-1
"""
cmolar(pressure, tair) = pressure / (R * (tair |> K))


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

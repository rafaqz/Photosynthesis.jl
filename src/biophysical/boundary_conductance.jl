"""
Boundary conductance models.

Provide parameters for [`boundary_conductance_free`](@ref), 
and [`boundary_conductance_forced`](@ref).
"""
abstract type AbstractBoundaryConductance end

"""
    BoundaryConductance(leafwidth)
    BoundaryConductance(; leafwidth=0.05u"m")

Standard boundary conductance formulation parameters.

$(FIELDDOCTABLE)
"""
@columns struct BoundaryConductance{LW} <: AbstractBoundaryConductance
    leafwidth::LW | 0.05 | m | (0.0, 1.0) | "Mean width of leaves"
end

"""
    boundary_conductance_free(f::AbstractBoundaryConductance, v)

Boundary layer conductance for heat - single sided, free convection,
given variables `v`.

Retuns conductance in `u"mol*m-2*s-1"`.  """
boundary_conductance_free(f::BoundaryConductance, v) =
    cmolar(v) * 0.5DHEAT * (grashof_number(v.tleaf, v.tair, f.leafwidth)^(1/4)) / f.leafwidth

"""
    boundary_conductance_forced(f::AbstractBoundaryConductance, v)

Boundary layer conductance for heat - single sided, forced convection, 
given variables `v`.

Retuns conductance in `u"mol*m-2*s-1"`.
"""
boundary_conductance_forced(f::BoundaryConductance, v) =
    # TODO name this 0.003
    0.003m*s^-1 * sqrt(ustrip(u"s^-1", v.windspeed / f.leafwidth)) * cmolar(v)

"""
    cmolar(pressure, airtemp)

Convert from m.s-1 to mol.m-2.s-1
"""
cmolar(pressure, airtemp) = pressure / (R * K(airtemp))
cmolar(v) = cmolar(v.pressure, v.tair) 

"""
    grashof_number(tleaf, tair, leafwidth)

Calculates the Grashof number given leaf temperature, 
air tempereature and leaf width.

return: dimensionless

See Leuning et al (1995) PCE 18:1183-1200 Eqn E4
"""
grashof_number(tleaf, tair, leafwidth) = 
    1.6e8m^-3*K^-1 * leafwidth^3 * abs(tleaf - tair)

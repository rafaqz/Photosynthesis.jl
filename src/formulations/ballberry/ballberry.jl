"""
Stomatal conductance submodels.

Many stomatal conductance formulations are identical to Ball-Berry, 
with changes only to this submodel.
"""
abstract type AbstractBallBerryStomatalConductanceSubModel <: AbstractStomatalConductanceSubModel end


# mixin macro that adds stomatal conductance fields to a struct
@mix @columns struct MixinBallBerryStomatalConductanceSubModel{μMoMo,F}
    gamma::μMoMo    | 0.0    | μmol*mol^-1  | (0.0, 10.0) | "Gamma for all Ball-Berry type models"
    g1::F           | 7.0    | _            | (0.0, 10.0) | "Slope parameter"
end

gamma(m::AbstractStomatalConductanceSubModel) = m.gamma
g1(m::AbstractStomatalConductanceSubModel) = m.g1


"""
    BallBerryStomatalConductanceSubModel(gamma, g1)

Basic Ball-Berry stomatal conductance formulation. 
Parameters are `gamma` in μmol*mol^-1 and `g1` scalar, 
which are also used for all Ball-Berry type models.

(modelgs = 2 in maestra)

$(FIELDDOCTABLE)
"""
@MixinBallBerryStomatalConductanceSubModel struct BallBerryStomatalConductanceSubModel{} <: AbstractStomatalConductanceSubModel end

gs_div_a(m::BallBerryStomatalConductanceSubModel, v) = g1(m) * v.rh / (v.cs - gamma(m))


"""
Ball-Berry implementation of stomatal conductance.
Most other models inherit fields and behaviours from this.
"""
abstract type AbstractBallBerryStomatalConductance <: AbstractStomatalConductance end

# Mixin fields for objects inheriting AbstractBallBerryStomatalConductance
@mix @columns struct MixinBallBerryStomatalConductance{MoM2S,GS,SM}
    g0::MoM2S       | 0.03                                   | mol*m^-2*s^-1 | (0.0, 0.2) | "Stomatal leakiness (gs when photosynthesis is zero)"
    gs_submodel::GS | BallBerryStomatalConductanceSubModel() | _             | _          | _
    soil_model::SM  | PotentialSoilMethod()                  | _             | _          | _
end

g0(m::AbstractBallBerryStomatalConductance) = m.g0
gs_submodel(m::AbstractBallBerryStomatalConductance) = m.gs_submodel
soil_model(m::AbstractBallBerryStomatalConductance) = m.soil_model


"""
    BallBerryStomatalConductance(g0, gs_submodel, soil_model)

Simple Ball-Berry stomatal conductance model.

$(FIELDDOCTABLE)
"""
@MixinBallBerryStomatalConductance struct BallBerryStomatalConductance{} <: AbstractBallBerryStomatalConductance end

"""
    BallBerryVars(tair, windspeed, par, rnet, soilmoist, pressure, tleaf, swp, swpshade, vpd, ca, rh,
                  cs, vpdleaf, rhleaf, fheat, gbhu, gbhf, gh, gbh, gv, gradn, lhv, et, slope, decoup,
                  gs_div_a, km, ci, gammastar, gs, gsv, gbv, jmax, vcmax, rd, ac, aleaf, vj, aj, fsoil)

Variables for Ball-Berry stomatal conductance models.

$(FIELDDOCTABLE)
"""
@MixinFvCBVars mutable struct BallBerryVars{} end


"""
    stomatal_conductance!(v, m::AbstractBallBerryStomatalConductance)

Stomatal conductance calculations.

Returns a tuple of leaf assimilation and stomatal conductance 
`(aleaf, gs)` in `u"μmol*m^-2*s^-1" and `u"mol*m^-2*s^-1"`
"""
function stomatal_conductance!(v, m::AbstractBallBerryStomatalConductance)
    v.fsoil = soilmoisture_conductance(soil_model(m), v)
    v.gs_div_a = gs_div_a(gs_submodel(m), v) * v.fsoil

    ac = rubisco_limited_rate(m, v, v.gs_div_a)
    aj = transport_limited_rate(m, v, v.gs_div_a)

    aleaf = min(ac, aj) - v.rd
    gs = max(g0(m) + v.gs_div_a * aleaf, g0(m))

    return aleaf, gs
end

"""
    rubisco_limited_rate(m::AbstractBallBerryStomatalConductance, v, gs_div_a)

Solution for assimilation when Rubisco activity is limiting, for all Ball-Berry models

Returns assimilation rate in `u"μmol*m^-2*s^-1"`
"""
function rubisco_limited_rate(m::AbstractBallBerryStomatalConductance, v, gs_div_a)
    a = g0(m) + gs_div_a * (v.vcmax - v.rd)
    b = (1.0 - v.cs * gs_div_a) * (v.vcmax - v.rd) + g0(m) * (v.km - v.cs) -
        gs_div_a * (v.vcmax * v.gammastar + v.km * v.rd)
    c = -(1.0 - v.cs * gs_div_a) * (v.vcmax * v.gammastar + v.km * v.rd) -
        g0(m) * v.km * v.cs
    cic = quad(Upper(), a, b, c)

    ac = if (cic <= zero(cic)) || (cic > v.cs)
        zero(v.vcmax)
    else
        v.vcmax * (cic - v.gammastar) / (cic + v.km)
    end
    return ac
end


"""
    transport_limited_rate(m::AbstractBallBerryStomatalConductance, v, gs_div_a)

Solution for assimilation rate when electron transport rate is limiting, 
for all Ball-Berry type models.

Returns assimilation rate in `u"μmol*m^-2*s^-1"`
"""
function transport_limited_rate(m::AbstractBallBerryStomatalConductance, v, gs_div_a)
    a = g0(m) + gs_div_a * (v.vj - v.rd)
    b = (1 - v.cs * gs_div_a) * (v.vj - v.rd) + g0(m) * (2v.gammastar - v.cs) -
        gs_div_a * (v.vj * v.gammastar + 2v.gammastar * v.rd)
    c = -(1 - v.cs * gs_div_a) * v.gammastar * (v.vj + 2v.rd) - g0(m) * 2v.gammastar * v.cs
    cij = quad(Upper(), a, b, c)
    aj = v.vj * (cij - v.gammastar) / (cij + 2v.gammastar)

    if (aj - v.rd < 1e-6oneunit(v.rd)) # Below light compensation point. TODO: why the magic number?
        cij = v.cs
        aj = v.vj * (cij - v.gammastar) / (cij + 2v.gammastar)
    end
    return aj
end

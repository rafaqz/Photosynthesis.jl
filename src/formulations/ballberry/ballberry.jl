abstract type AbstractGSsubModel end

# Ball-Berry implementation of stomatal conductance

" @Gs mixin macro adds stomatal conductance fields to a struct"
@mix @columns struct MixinBallBerryGs{μMoMo,F}
    gamma::μMoMo    | 0.0    | μmol*mol^-1  | (0.0, 10.0) | "Gamma for all Ball-Berry type models"
    g1::F           | 7.0    | _            | (0.0, 10.0) | "Slope parameter"
end

"""
Ball-Berry stomatal conductance formulation parameters
`gamma` in μmol*mol^-1 and `g1`scalar, also used for all Ball-Berry type models.
(modelgs = 2 in maestra)
"""
@MixinBallBerryGs struct BallBerryGSsubModel{} <: AbstractGSsubModel end

"""
    gsdiva(::BallBerryGSsubModel, v) """
gsdiva(f::BallBerryGSsubModel, v) = f.g1 * v.rh / (v.cs - f.gamma)
# Ball-Berry implementation of photosynthesis.
# Most other models inherit fields and behaviours from this.



# Ball-Berry abstract/mixin implementation of photosynthesis

abstract type AbstractBallBerryStomatalConductance <: AbstractStomatalConductance end

@mix @columns struct MixinBallBerryStomCond{MoM2S,GS,SM}
    g0::MoM2S       | 0.03                  | mol*m^-2*s^-1 | (0.0, 0.2) | "Stomatal leakiness (gs when photosynthesis is zero)"
    gs_submodel::GS | BallBerryGSsubModel() | _             | _          | _
    soilmethod::SM  | PotentialSoilMethod() | _             | _          | _
end

@MixinBallBerryStomCond struct BallBerryStomatalConductance{} <: AbstractBallBerryStomatalConductance end

@MixinFvCBVars mutable struct BallBerryVars{} end


update_extremes!(f::AbstractBallBerryStomatalConductance, v) = begin
    v.aleaf = -v.rd
    v.gs = f.g0
end


"""
    stomatal_conductance!(f::AbstractBallBerryStomatalConductance, v)
Stomatal conductance calculations.
"""
function stomatal_conductance!(f::AbstractBallBerryStomatalConductance, v)
    v.fsoil = soilmoisture_conductance(f.soilmethod, v)
    v.gsdiva = gsdiva(f.gs_submodel, v) * v.fsoil

    ac = rubisco_limited_rate(f, v, v.gsdiva)
    aj = transport_limited_rate(f, v, v.gsdiva)

    aleaf = min(ac, aj) - v.rd
    gs = max(f.g0 + v.gsdiva * aleaf, f.g0)
    # println("assimilation: ", aleaf, "\nrubisco limit: ", ac, "\ntransport limit: ", aj, "\nrespiration: ", v.rd, "\n")

    aleaf, gs
end

"""
    rubisco_limited_rate(f::AbstractBallBerryStomatalConductance, v)
Solution when Rubisco activity is limiting for all Ball-Berry models
"""
function rubisco_limited_rate(f::AbstractBallBerryStomatalConductance, v, gsdiva)
    a = f.g0 + gsdiva * (v.vcmax - v.rd)
    b = (1.0 - v.cs * gsdiva) * (v.vcmax - v.rd) + f.g0 * (v.km - v.cs) -
        gsdiva * (v.vcmax * v.gammastar + v.km * v.rd)
    c = -(1.0 - v.cs * gsdiva) * (v.vcmax * v.gammastar + v.km * v.rd) -
        f.g0 * v.km * v.cs
    cic = quad(Upper(), a, b, c)

    if (cic <= zero(cic)) || (cic > v.cs)
        ac = zero(v.vcmax)
    else
        ac = v.vcmax * (cic - v.gammastar) / (cic + v.km)
    end
    ac
end


"""
    transport_limited_rate(f::AbstractBallBerryStomatalConductance, v)
Solution for when electron transport rate is limiting for all Ball-Berry type models
"""
function transport_limited_rate(f::AbstractBallBerryStomatalConductance, v, gsdiva)
    a = f.g0 + gsdiva * (v.vj - v.rd)
    b = (1 - v.cs * gsdiva) * (v.vj - v.rd) + f.g0 * (2v.gammastar - v.cs) -
        gsdiva * (v.vj * v.gammastar + 2v.gammastar * v.rd)
    c = -(1 - v.cs * gsdiva) * v.gammastar * (v.vj + 2v.rd) - f.g0 * 2v.gammastar * v.cs
    cij = quad(Upper(), a, b, c)
    aj = v.vj * (cij - v.gammastar) / (cij + 2v.gammastar)

    if (aj - v.rd < 1e-6oneunit(v.rd)) # Below light compensation point. TODO: why the magic number?
        cij = v.cs
        aj = v.vj * (cij - v.gammastar) / (cij + 2v.gammastar)
    end

    aj
end

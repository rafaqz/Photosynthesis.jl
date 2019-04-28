
# Ball-Berry implementation of stomatal conductance

" @Gs mixin macro adds stomatal conductance fields to a struct"
@mix @columns struct MixinBallBerryGs{μMoMo,F}
    gamma::μMoMo    | 0.0    | μmol*mol^-1    | Gamma(1, 2)     | (0.0, 10.0) | "Gamma for all Ball-Berry type models"
    g1::F           | 7.0    | _              | Gamma(10, 7/10) | (0.0, 10.0) | "Slope parameter"
end

"""
Ball-Berry stomatal conductance formulation parameters
`gamma` in μmol*mol^-1 and `g1`scalar, also used for all Ball-Berry type models.
(modelgs = 2 in maestra)
"""
@MixinBallBerryGs struct BallBerryStomatalConductance{} <: AbstractStomatalConductance end

"""
    gsdiva(::BallBerryStomatalConductance, v) """
gsdiva(f::BallBerryStomatalConductance, v) = f.g1 * v.rh / (v.cs - f.gamma) * v.fsoil
# Ball-Berry implementation of photosynthesis.
# Most other models inherit fields and behaviours from this.



# Ball-Berry abstract/mixin implementation of photosynthesis

abstract type AbstractBallBerryPhotosynthesis <: AbstractFvCBPhotosynthesis end

@MixinFvCBPhoto @mix struct MixinBallBerryPhoto{MoM2S}
    g0::MoM2S         | 0.03 | mol*m^-2*s^-1 | Gamma(10, 0.03/10) | (0.0, 0.2) | "Stomatal leakiness (gs when photosynthesis is zero)."
end

@MixinBallBerryPhoto struct BallBerryPhotosynthesis{} <: AbstractBallBerryPhotosynthesis end

@MixinFvCBVars mutable struct BallBerryVars{} end


photo_init!(::AbstractBallBerryPhotosynthesis, v) = v.rhleaf = v.rh

photo_update!(::AbstractBallBerryPhotosynthesis, v, tleaf) =
    v.rhleaf = 1.0 - v.vpdleaf / saturated_vapour_pressure(tleaf)

update_extremes!(f::AbstractBallBerryPhotosynthesis, v) = begin
    v.aleaf = -v.rd
    v.gs = f.g0
end


"""
    stomatal_conductance!(f::AbstractBallBerryPhotosynthesis, v)
Stomatal conductance calculations.
"""
function stomatal_conductance!(f::AbstractBallBerryPhotosynthesis, v)
    v.gsdiva = gsdiva(f.gsmodel, v)

    v.ac = rubisco_limited_rate(f, v)
    v.aj = transport_limited_rate(f, v)
    v.aleaf = min(v.ac, v.aj) - v.rd

    v.gs = max(f.g0 + v.gsdiva * v.aleaf, f.g0)
    nothing
end

"""
    rubisco_limited_rate(f::AbstractBallBerryPhotosynthesis, v)
Solution when Rubisco activity is limiting for all Ball-Berry models
"""
function rubisco_limited_rate(f::AbstractBallBerryPhotosynthesis, v)
    a = f.g0 + v.gsdiva * (v.vcmax - v.rd)
    b = (1.0 - v.cs * v.gsdiva) * (v.vcmax - v.rd) + f.g0 * (v.km - v.cs) -
        v.gsdiva * (v.vcmax * v.gammastar + v.km * v.rd)
    c = -(1.0 - v.cs * v.gsdiva) * (v.vcmax * v.gammastar + v.km * v.rd) -
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
    transport_limited_rate(f::AbstractBallBerryPhotosynthesis, v)
Solution for when electron transport rate is limiting for all Ball-Berry type models
"""
function transport_limited_rate(f::AbstractBallBerryPhotosynthesis, v)
    a = f.g0 + v.gsdiva * (v.vj - v.rd)
    b = (1.0 - v.cs * v.gsdiva) * (v.vj - v.rd) + f.g0 * (2.0v.gammastar - v.cs) -
        v.gsdiva * (v.vj * v.gammastar + 2.0v.gammastar * v.rd)
    c = -(1.0 - v.cs * v.gsdiva) * v.gammastar * (v.vj + 2.0v.rd) - f.g0 * 2.0v.gammastar * v.cs
    cij = quad(Upper(), a, b, c)
    aj = v.vj * (cij - v.gammastar) / (cij + 2.0v.gammastar)

    if (aj - v.rd < 1e-6oneunit(v.rd)) # Below light compensation point
        cij = v.cs
        aj = v.vj * (cij - v.gammastar) / (cij + 2.0v.gammastar)
    end

    aj
end

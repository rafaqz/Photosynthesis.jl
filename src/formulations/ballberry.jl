"""
Ball-Berry stomaratl conductance formulation parameters
`gamma` in μmol*mol^-1 and `g1`scalar, also used for all Ball-Berry type models.
(modelgs = 2 in maestra)
"""
@Gs struct BallBerryStomatalConductance{} <: AbstractStomatalConductance end


abstract type AbstractBallBerryModel{GS,SM} <: AbstractPhotoModel end

" @BB mixin macro adds base Balll-Berry fields to a struct"
@mix @columns struct BB{MoM2S}
    g0::MoM2S | 0.03 | mol*m^-2*s^-1 | Gamma(10, 0.03/10) | (0.0, 0.2) | "Stomatal leakiness (gs when photosynthesis is zero)."
end

"""
General Ball-Berry stomatal conductance model.

THis model has various submodels and soil water methods defined in gsmodel
and soilmethod.
"""
@BB struct BallBerryModel{GS,SM} <: AbstractBallBerryModel{GS,SM}
    gsmodel::GS    | BallBerryStomatalConductance() | _ | _ | _ | _
    soilmethod::SM | ConstantSoilMethod()           | _ | _ | _ | _
end

@Vars mutable struct BallBerryVars{} end


photo_init!(::AbstractBallBerryModel, v) = v.rhleaf = v.rh

photo_update!(::AbstractBallBerryModel, v, tleaf) = 
    v.rhleaf = 1.0 - v.vpdleaf / saturated_vapour_pressure(tleaf)


extremes!(f::AbstractBallBerryModel, v) = begin
    v.aleaf = -v.rd
    v.gs = f.g0
end


"""
    stomatal_conductance!(f::BallBerryModel, v, p)
Stomatal conductance calculations for the Jarvis model
"""
function stomatal_conductance!(f::BallBerryModel, v, p)
    v.gsdiva = gsdiva(f.gsmodel, v, p)
    assimilation!(v, p)
    v.gs = max(f.g0 + v.gsdiva * v.aleaf, f.g0)
    nothing
end


""" 
    rubisco_limited_rate(f::AbstractBallBerryModel, v, p)
Solution when Rubisco activity is limiting for all Ball-Berry models
"""
function rubisco_limited_rate(f::AbstractBallBerryModel, v, p)
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
    transport_limited_rate(f::AbstractBallBerryModel, v, p)
Solution for when electron transport rate is limiting for all Ball-Berry type models
"""
function transport_limited_rate(f::AbstractBallBerryModel, v, p)
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

"""
    gsdiva(::BallBerryStomatalConductance, v, p) """
gsdiva(f::BallBerryStomatalConductance, v, p) = 
    f.g1 * v.rh / (v.cs - f.gamma) * v.fsoil

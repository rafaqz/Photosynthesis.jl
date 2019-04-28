
# Tuzet mplementation of stomatal conductance

""" Tuzet stomatal conductance formulation parameters
(modelgs = 5 in maestra)
"""
@MixinBallBerryGs struct TuzetStomatalConductance{} <: AbstractStomatalConductance end

"""
    gsdiva(::TuzetStomatalConductance, v) """
gsdiva(f::TuzetStomatalConductance, v) =
    (f.g1 / (v.cs - f.gamma)) * fpsil(v.psilin, v.sf, v.psiv)

"""
    fpsil(psil, sf, psiv)
Tuzet et al. 2003 leaf water potential function
"""
fpsil(psil, sf, psiv) = (1 + exp(sf * psiv)) / (1 + exp(sf * (psiv - psil)))


# Tuzet mplementation of soil method

@MixinMaespaSoil struct TuzetSoilMethod{} <: AbstractSoilMethod{T} end

soilmoisture_conductance!(f::TuzetSoilMethod, v) = begin
    # What is this from? v.psif = (1 + exp(v.sf * v.psiv)) / (1 + exp(v.sf * (v.psiv - v.psil)))
    pd = non_stomatal_potential_dependence(f.non_stomatal, v.weightedswp)
    v.vcmax *= pd
    v.jmax *= pd
end



# Implementation of photosynthesis

"""
Tuzet photosynthesis model.

This model is limited to using `TuzetStomatalConductance` and `TuzetSoilMethods`.
"""
@MixinMaespaPhoto struct TuzetPhotosynthesis{} <: AbstractMaespaPhotosynthesis end

@redefault struct TuzetPhotosynthesis
    gsmodel    | TuzetStomatalConductance()
    soilmethod | TuzetSoilMethod()
end

@MixinMaespaVars mutable struct TuzetVars{kPa,pkPa}
    psilin::kPa      | -999.0      | kPa                   | _
    psiv::kPa        | -1.9        | kPa                   | _
    sf::pkPa         | 3.2         | kPa^-1                | _
end



# Tuzet implementation of energy balance
# Essentially a wrapper. Runs a standard FvCB energy with Tuzet 
# stomatal conductance and soil moisture routine in a root finder

@default struct TuzetEnergyBalance <: AbstractEnergyBalance 
    energy_balance | FvCBEnergyBalance(photosynthesis = TuzetPhotosynthesis())
end

"""
    run_enbal!(v, pm::TuzetEnergyBalance)
Runs energy balance inside a root finding algorithm to calculate leaf water potential.
"""
function enbal!(m::TuzetEnergyBalance, v, p)
    v.psilin = -0.1oneunit(v.psilin)
    v.psil = -0.1oneunit(v.psil)
    bracket = -100.0oneunit(v.psil), zero(v.psil)
    tolz = 1e-03oneunit(v.psil)

    v.psilin = find_zero(x -> leaf_water_potential_finder(x, v, p), bracket, FalsePosition(), atol=tolz)
    nothing
end

"""
    leaf_water_potential_finder(psilin, v, p)
PSIL finder. A wrapper function that
return the squared difference in PSILIN and PSIL
"""
function leaf_water_potential_finder(psilin, v, p)
    v.ci = zero(v.ci)
    enbal!(p.energy_balance, v)
    v.psilin - v.psil
end

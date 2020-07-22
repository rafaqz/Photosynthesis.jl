
# Tuzet implementation of stomatal conductance ---------------------------------------

""" 
    TuzetStomatalConductanceSubModel(gamma, g1)

Tuzet stomatal conductance formulation parameters. (modelgs = 5 in maestra)

$(FIELDDOCTABLE)
"""
@MixinBallBerryStomatalConductanceSubModel struct TuzetStomatalConductanceSubModel{} <: AbstractStomatalConductanceSubModel end

gs_div_a(f::TuzetStomatalConductanceSubModel, v) = (f.g1 / (v.cs - f.gamma)) * fpsil(v.psilin, v.sf, v.psiv)

"""
    fpsil(psil, sf, psiv)

Tuzet et al. 2003 leaf water potential function
"""
fpsil(psil, sf, psiv) = (1 + exp(sf * psiv)) / (1 + exp(sf * (psiv - psil)))


# Tuzet mplementation of soil method

@default_kw struct TuzetSoilMethod{T,NS} <: AbstractSoilMethod
    soilmethod::T    | ConstantSoilMethod()
    potentialdependence_model::NS | ZhouPotentialDependence()
end

# TODO: Fix this
# soilmoisture_conductance!(v, f::TuzetSoilMethod) = begin
    # What is this from? v.psif = (1 + exp(v.sf * v.psiv)) / (1 + exp(v.sf * (v.psiv - v.psil)))
    # non_stomatal_potential_dependence(f.potentialdependence_model, v.weightedswp)
# end

"""
    TuzetStomatalConductance(gs_submodel, soilmethod)

Tuzet stomatal conductance model.

This model is limited to using `TuzetStomatalConductance` and `TuzetSoilMethods`.

$(FIELDDOCTABLE)
"""
@MixinBallBerryStomatalConductance struct TuzetStomatalConductance{} <: AbstractBallBerryStomatalConductance end

# Change defaults for gs_submodel and soilmethod
@default TuzetStomatalConductance begin
    gs_submodel | TuzetStomatalConductanceSubModel()
    soilmethod  | TuzetSoilMethod()
end

"""
    TuzetVars()

Variables object for Tuzet models.

$(FIELDDOCTABLE)
"""
@MixinFvCBVars mutable struct TuzetVars{kPa,pkPa}
#   Field       | Default | Units  | 
    psilin::kPa | -999.0  | kPa    | _
    psiv::kPa   | -1.9    | kPa    | _
    sf::pkPa    | 3.2     | kPa^-1 | _
end


# Tuzet implementation of energy balance ---------------------------------------------

"""
    TuzetEnergyBalance(energy_balance)
    TuzetEnergyBalance(; energy_balance)

Tuzet implementation of energy balance.
Essentially just a wrapper that runs a standard FvCB energy balance, with Tuzet 
stomatal conductance formulation and soil moisture routines in a root finder.

$(FIELDDOCTABLE)
"""
@default struct TuzetEnergyBalance{EB} <: AbstractFvCBEnergyBalance 
    energy_balance::EB | FvCBEnergyBalance(photosynthesis=FvCBPhotosynthesis(stomatal_conductance=TuzetStomatalConductance))
end

"""
    enbal!(v, pm::TuzetEnergyBalance)

Runs energy balance inside a root finding algorithm to calculate leaf water potential.
"""
function enbal!(v, p::TuzetEnergyBalance)
    v.psilin = -0.1oneunit(v.psilin)
    v.psil = -0.1oneunit(v.psil)
    bracket = -100.0oneunit(v.psil), zero(v.psil)
    tolz = 1e-03oneunit(v.psil)

    v.psilin = findzero(x -> leaf_water_potential_finder(p.energy_balance, v, x), bracket; atol=tolz)
    nothing
end

"""
    leaf_water_potential_finder(p, v, psilin)

PSIL finder. A wrapper function for the energy balance that returns the squared 
difference in PSILIN and PSIL, for use in a root finder.
"""
function leaf_water_potential_finder(p, v, psilin)
    v.ci = zero(v.ci)
    enbal!(v, p.energy_balance)
    v.psilin - v.psil
end

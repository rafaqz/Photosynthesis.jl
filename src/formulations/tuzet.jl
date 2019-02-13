# Tuzet model.

@MaespaSoil struct TuzetSoilMethod{} <: AbstractSoilMethod{T} end

""" Tuzet stomatal conductance formulation parameters
(modelgs = 5 in maestra)
"""
@Gs struct TuzetStomatalConductance{} <: AbstractStomatalConductance end

"""
Tuzet photosynthesis model.

This model is limited to using `TuzetStomatalConductance` and `TuzetSoilMethods`.
"""
@Maespa struct TuzetModel{GS<:TuzetStomatalConductance,
                          SM<:TuzetSoilMethod{<:AbstractPotentialDependence}
                         } <: AbstractMaespaModel
    gsmodel::GS    | TuzetStomatalConductance() | _ | _ | _ | _
    soilmethod::SM | TuzetSoilMethod()          | _ | _ | _ | _
end

@MaespaVars mutable struct TuzetVars{kPa,pkPa}
    psilin::kPa      | -999.0      | kPa                   | _
    psiv::kPa        | -1.9        | kPa                   | _
    sf::pkPa         | 3.2         | kPa^-1                | _
end


soilmoisture_conductance!(f::TuzetSoilMethod, v) = begin
    # What is this from? v.psif = (1 + exp(v.sf * v.psiv)) / (1 + exp(v.sf * (v.psiv - v.psil)))
    vjmod = potential_dependence(f.soildata, v.weightedswp)
    v.vcmax *= vjmod
    v.jmax *= vjmod
    v.fsoil = oneunit(v.fsoil)
end

"""
    gsdiva(::TuzetStomatalConductance, v) """
gsdiva(f::TuzetStomatalConductance, v) =
    (f.g1 / (v.cs - f.gamma)) * fpsil(v.psilin, v.sf, v.psiv)

""" 
    fpsil(psil, sf, psiv)
Tuzet et al. 2003 leaf water potential function 
"""
fpsil(psil, sf, psiv) = (1 + exp(sf * psiv)) / (1 + exp(sf * (psiv - psil)))

"""
    run_phototrans!(v, pm::TuzetModel)
Potosynthesis runner for `TuzetModel` runs phototranspiration inside a root
finding algorithm to calculate leaf water potential
"""
function run_phototrans!(m::TuzetModel, v, p)
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
    phototranspiration!(v, p)
    v.psilin - v.psil
end

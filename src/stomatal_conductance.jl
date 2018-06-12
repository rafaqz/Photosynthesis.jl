"Deal with extremes where jmax or vcmax are less than zero"
function extremes! end

extremes!(v, p) = begin
    if v.jmax <= zero(v.jmax) || v.vcmax <= zero(v.vcmax)
        println("extreme")
        extremes!(v, p.model, p)
        return true
    end
    false
end
extremes!(v, f::JarvisModel, p) = begin
    v.aleaf = -v.rd
    v.gs = f.gsmin
end
extremes!(v, f::AbstractBallBerryModel, p) = begin
    v.aleaf = -v.rd
    v.gs = f.g0
end

"""
    stomatal_conductance!(v, ::AbstractPhotoModel, p)
Stomatal conductance and intercellular CO2 partial pressure calculations.

v.aleaf is NET leaf photosynthesis.
"""
function stomatal_conductance!(v, f, p) end

"""
    stomatal_conductance!(v, f::JarvisModel, p)
Stomatal conductance and for the Jarvis model
"""
function stomatal_conductance!(v, f::JarvisModel, p)
    v.gs = factor_conductance(f, v, p)
    assimilation_gs_known!(v, p)
    nothing
end

"""
    stomatal_conductance!(v, f::BallBerryModel, p)
Stomatal conductance calculations for the Jarvis model
"""
function stomatal_conductance!(v, f::BallBerryModel, p)
    assimilation_gs_unknown!(v, p)
    v.gs = max(f.g0 + v.gsdiva * v.aleaf, f.g0)
    nothing
end

"""
    stomatal_conductance!(v, f::TuzetModel, p)
Stomatal conductance calculations for the Jarvis model
"""
function stomatal_conductance!(v, f::AbstractMaespaModel, p)
    assimilation_gs_unknown!(v, p)

    v.gs = shape_gs(f.gsshape, v, f)
    psil_gs!(v, f, p)

    if v.gs > zero(v.gs) && v.aleaf > zero(v.aleaf)
        v.ci = v.cs - v.aleaf / v.gs
    else
        v.ci = v.cs
    end

    nothing
end

function shape_gs(f::HardMinimumGS, v, p) 
    gs = p.g0 + v.gsdiva * v.aleaf
    max(p.g0, gs)
end

""" 
Hyperbolic minimum; use for optimization to avoid discontinuity.
HMSHAPE determines how smooth the transition is.
"""
function shape_gs(f::HyperbolicMinimumGS, v, p) 
    aleafhypmin = (ac + aj - sqrt((v.ac + v.aj)^2 - 4f.hmshape * v.ac * v.aj)) /
              (2 * f.hmshape) - v.rd
    p.g0 + v.gsdiva * aleafhypmin
end

function psil_gs!(v, f, p) end

"""
Stomatal conductance calculations for the Emax model.
Uses the Ball-Berry model but contains transpiration rate calculations calculated 
in other ways in other models 
"""
function psil_gs!(v, f::EmaxModel, p)
    # Maximum transpiration rate
    v.emaxleaf = v.ktot * (v.weightedswp - v.minleafwp)
    # Leaf transpiration in mmol m-2 s-1 - ignoring boundary layer effects
    etest = (v.vpd / v.pressure) * v.gs * GSVGSC
    # Leaf water potential
    v.psil = v.weightedswp - etest / v.ktot

    if etest > v.emaxleaf
       v.fsoil = v.emaxleaf / etest # Just for output
       gsv = max(v.emaxleaf / (v.vpd / v.pressure), f.g0)
       # From yplantmc
       v.gs = gsv / GSVGSC
       # Minimum leaf water potential reached: recalculate psil
       v.psil = v.weightedswp - v.emaxleaf / v.ktot
       # Recalculate aleaf using jarvis method
       assimilation_gs_known!(v, f)
    end
end

""" Ball-Berry assimilation, gs not known"""
function assimilation_gs_unknown!(v, p)
    v.gsdiva = gsdiva(p.model.gsmodel, v, p)
    v.ac = rubisco_limited_rate(p.model, v, p)
    v.aj = transport_limited_rate(p.model, v, p)
    v.aleaf = min(v.ac, v.aj) - v.rd
end

""" Jarvis assimilation where we already know gs. Also used in the emax model """
function assimilation_gs_known!(v, p)
    v.ac = rubisco_limited_rate(v, p)
    v.aj = transport_limited_rate(v, p)
    v.aleaf = min(v.ac, v.aj)
end


""" 
    rubisco_limited_rate(v, p)
Solution when Rubisco activity is limiting for the Jarvis model """
function rubisco_limited_rate(v, p)
    a = 1.0 / v.gs
    b = (v.rd - v.vcmax) / v.gs - v.cs - v.km
    c = v.vcmax * (v.cs - v.gammastar) - v.rd * (v.cs + v.km)
    quad(Val{:lower}, a, b, c)
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
    cic = quad(Val{:upper}, a, b, c)

    if (cic <= zero(cic)) || (cic > v.cs)
        ac = zero(v.vcmax)
    else
        ac = v.vcmax * (cic - v.gammastar) / (cic + v.km)
    end
    ac
end

""" 
    transport_limited_rate(f::JarvisModel, v, p)
Solution when electron transport rate is limiting for the Jarvis model
"""
function transport_limited_rate(v, p)
    a = 1.0 / v.gs
    b = (v.rd - v.vj) / v.gs - v.cs - 2v.gammastar
    c = v.vj * (v.cs - v.gammastar) - v.rd * (v.cs + 2v.gammastar)
    quad(Val{:lower}, a, b, c)
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
    cij = quad(Val{:upper}, a, b, c)
    aj = v.vj * (cij - v.gammastar) / (cij + 2.0v.gammastar)


    if (aj - v.rd < 1e-6oneunit(v.rd)) # Below light compensation point
        cij = v.cs
        aj = v.vj * (cij - v.gammastar) / (cij + 2.0v.gammastar)
    end

    aj
end

"""
    gsdiva(::AbstractModelGS, v, p)
Formulation-specific component for the Ball-Berry family of stomatal conductance models.
"""
function gsdiva end

"""
    gsdiva(::BallBerryStomatalConductance, v, p) """
gsdiva(f::BallBerryStomatalConductance, v, p) = 
    f.g1 * v.rh / (v.cs - f.gamma) * v.fsoil
"""
    gsdiva(::LeuningStomatalConductance, v, p) 
From R. Leuning, A critical appraisal of a combined stomatal-photosynthesis
model for C3 plants. Plant, Celt and Environment (1995) 18, 339-355
"""
gsdiva(f::LeuningStomatalConductance, v, p) =
    f.g1 / ( v.cs - f.gamma) / (1.0 + v.vpdleaf / f.d0l) * v.fsoil

"""
    gsdiva(::MedlynStomatalConductance, v, p) """
gsdiva(f::MedlynStomatalConductance, v, p) =
    (1.0 + (f.g1 * v.fsoil) / sqrt(max(f.vpdmin, v.vpdleaf)u"Pa^-1")) / v.cs

"""
    gsdiva(::ThreeParStomatalConductance, v, p) """
gsdiva(f::ThreeParStomatalConductance, v, p) =
    f.g1 / (v.cs - f.gamma) / (max(f.vpdmin, v.vpdleaf)u"Pa^-1")^f.gk

"""
    gsdiva(::TuzetStomatalConductance, v, p) """
gsdiva(f::TuzetStomatalConductance, v, p) =
    (f.g1 / (v.cs - f.gamma)) * fpsil(v.psilin, v.sf, v.psiv)


""" 
    fpsil(psil, sf, psiv)
Tuzet et al. 2003 leaf water potential function """
fpsil(psil, sf, psiv) = (1 + exp(sf * psiv)) / (1 + exp(sf * (psiv - psil)))

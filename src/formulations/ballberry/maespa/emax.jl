
# Emax implementation of soil method

@MixinMaespaSoil struct EmaxSoilMethod{T} <: AbstractSoilMethod{T} end

soilmoisture_conductance!(f::EmaxSoilMethod, v) = begin
    pd = non_stomatal_potential_dependence(f.non_stomatal, v.weightedswp)
    v.vcmax *= pd
    v.jmax *= pd
    v.fsoil = oneunit(v.fsoil)
end

# Emax implementation of photosynthesis

"""
Emax photosynthesis inherited from Maespa.
The same options are available for specialised stomatal conducance,
but using emax soil water methods.
"""
@MixinMaespaPhoto struct EmaxPhotosynthesis{} <: AbstractMaespaPhotosynthesis end

struct JarvisMode <: AbstractJarvisPhotosynthesis end

# Set specific defaults for Emax stomatal conductance and soil method
@redefault struct EmaxPhotosynthesis
    gsmodel    | BallBerryStomatalConductance()
    soilmethod | EmaxSoilMethod()             
end

@MixinMaespaVars mutable struct EmaxVars{kPa,mMoM2S}
    minleafwp::kPa   | 0.1         | kPa                   | _
    emaxleaf::mMoM2S | 400.0       | mmol*m^-2*s^-1        | _
end


"""
    stomatal_conductance!(f::AbstractMaespaPhotosynthesis, v)
Stomatal conductance calculations for the Jarvis model
"""
function stomatal_conductance!(f::EmaxPhotosynthesis, v)
    v.gsdiva = gsdiva(f.gsmodel, v)

    # Maximum transpiration rate
    emaxleaf = v.ktot * (v.weightedswp - v.minleafwp)

    # Leaf transpiration: ignoring boundary layer effects
    etest = (v.vpd / v.pressure) * v.gs * GSVGSC

    # Leaf water potential
    v.psil = v.weightedswp - etest / v.ktot

    if etest > emaxleaf
        # Just for output
        v.fsoil = emaxleaf / etest

        gsv = emaxleaf / (v.vpd / v.pressure)
        v.gs = gsv / GSVGSC

        # Minimum leaf water potential reached, recalculate psil
        v.psil = v.weightedswp - emaxleaf / v.ktot

        v.ac = rubisco_limited_rate(JarvisMode(), v)
        # v.aj = transport_limited_rate(JarvisMode(), v)
        # v.aleaf = min(v.ac, v.aj) - v.rd

        # v.gs = shape_gs(f.gsshape, v, f)

    end

    nothing
end

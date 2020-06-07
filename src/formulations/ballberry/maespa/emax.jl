# Emax implementation of soil method
@udefault_kw struct EmaxSoilMethod{T,NS} <: AbstractSoilMethod 
    soilmethod::T    | ConstantSoilMethod()
    non_stomatal::NS | ZhouPotentialDependence()
end

soilmoisture_conductance!(v, m::EmaxSoilMethod) = begin
    pd = non_stomatal_potential_dependence(m.non_stomatal, v.swp)
    v.vcmax *= pd
    v.jmax *= pd
    oneunit(v.fsoil)
end

# Emax implementation of stomatal conductance
abstract type AbstractEmaxStomatalConductance <: AbstractBallBerryStomatalConductance end

"""
Emax photosynthesis inherited from Maespa.
The same options are available for specialised stomatal conducance,
but using emax soil water methods.
"""
@MixinBallBerryStomatalConductance struct EmaxStomatalConductance{SH} <: AbstractEmaxStomatalConductance 
    gsshape::SH    | HardMinimum() | _  | _ | _
end

gsshape(m::AbstractEmaxStomatalConductance) = m.gsshape

struct JarvisMode <: AbstractJarvisStomatalConductance end


@MixinFvCBVars mutable struct EmaxVars{M,EL,KT,P}
    minleafwp::M   | 0.1         | kPa                   | _
    emaxleaf::EL   | 400.0       | mmol*m^-2*s^-1        | _
    ktot::KT       | 2.0         | mmol*m^-2*s^-1*MPa^-1 | _
    psil::P        | -111.0      | kPa                   | _
end


"""
    stomatal_conductance!(v, f::AbstractEmaxStomatalConductance)
Stomatal conductance calculations for the Jarvis model
"""
function stomatal_conductance!(v, m::AbstractEmaxStomatalConductance)
    v.gsdiva = gsdiva(gs_submodel(m), v)

    # Maximum transpiration rate
    emaxleaf = v.ktot * (v.swp - v.minleafwp)

    # Leaf transpiration: ignoring boundary layer effects
    etest = (v.vpd / v.pressure) * v.gs * GSVGSC

    # Leaf water potential
    v.psil = v.swp - etest / v.ktot

    if etest > emaxleaf
        # Just for output
        v.fsoil = emaxleaf / etest

        gsv = emaxleaf / (v.vpd / v.pressure)
        v.gs = gsv / GSVGSC

        # Minimum leaf water potential reached, recalculate psil
        v.psil = v.swp - emaxleaf / v.ktot

        v.ac = rubisco_limited_rate(JarvisMode(), v)
        v.aj = transport_limited_rate(JarvisMode(), v)
        aleaf = min(v.ac, v.aj) - v.rd

        gs = shape_gs(gsshape(m), v, m)
    end

    aleaf, gs 
end


abstract type AbstractEmaxEnergyBalance <: AbstractFvCBEnergyBalance end

@columns struct EmaxEnergyBalance{EB,SR,PK} <: AbstractEmaxEnergyBalance 
    energy_balance_model::EB | FvCBEnergyBalance(photosynthesis=FvCBPhotosynthesis(
                                               stomatal_conductance=EmaxStomatalConductance(
                                                   soilmethod=EmaxSoilMethod()))) | _ | _ | _
    totsoilres::SR           | 0.5  | m^2*s^1*MPa^1*mmol^-1 | (0.0, 10.0) | _
    plantk::PK               | 3.0  | mmol*m^-2*s^-1*MPa^-1 | (0.0, 10.0) | _
end

energy_balance_model(m::AbstractEmaxEnergyBalance) = m.energy_balance_model
totsoilres(m::AbstractEmaxEnergyBalance) = m.totsoilres
plantk(m::AbstractEmaxEnergyBalance) = m.plantk

enbal!(v, m::AbstractEmaxEnergyBalance) = begin
    v.ktot = 10 / (totsoilres(m) + 1.0 / plantk(m))
    enbal!(v, energy_balance_model(m))
end

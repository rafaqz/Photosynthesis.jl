"""
    EmaxSoilMethod(soilmethod, non_stomatal)

Emax implementation of soil method.

$(FIELDDOCTABLE)
"""
@default_kw struct EmaxSoilMethod{T,NS} <: AbstractSoilMethod
    soilmethod::T    | ConstantSoilMethod()
    non_stomatal::NS | ZhouPotentialDependence()
end

soilmoisture_conductance!(v, m::EmaxSoilMethod) = begin
    pd = non_stomatal_potential_dependence(m.non_stomatal, v.swp)
    v.vcmax *= pd
    v.jmax *= pd
    oneunit(v.fsoil)
end

"""
    EmaxStomatalConductance(BallBerryStomatalConductance(gsshape, g0, gs_submodel, soil_model)

The same options are available for specialised stomatal conducance,
but using emax soil water methods.

WARNING: Currently not passing tests

$(FIELDDOCTABLE)
"""
@MixinBallBerryStomatalConductance struct EmaxStomatalConductance{SH} <: AbstractBallBerryStomatalConductance
    gsshape::SH    | HardMinimum() | _  | _ | _
end

gsshape(m::EmaxStomatalConductance) = m.gsshape

struct JarvisMode <: AbstractJarvisStomatalConductance end

function stomatal_conductance!(v, m::EmaxStomatalConductance)
    v.gs_div_a = gs_div_a(gs_submodel(m), v)

    # Maximum transpiration rate
    emaxleaf = v.ktot * (v.swp - v.minleafwp)

    # Leaf transpiration: ignoring boundary layer effects
    etest = (v.vpdleaf / v.pressure) * v.gs * GSVGSC

    # Leaf water potential
    v.psil = v.swp - etest / v.ktot

    if etest > emaxleaf
        # Just for output
        v.fsoil = emaxleaf / etest

        gsv = emaxleaf / (v.vpdleaf / v.pressure)
        v.gs = gsv / GSVGSC

        # Minimum leaf water potential reached, recalculate psil
        v.psil = v.swp - emaxleaf / v.ktot

        v.ac = rubisco_limited_rate(JarvisMode(), v)
        v.aj = transport_limited_rate(JarvisMode(), v)
        v.aleaf = min(v.ac, v.aj) - v.rd

        v.gs = shape_gs(gsshape(m), v, m)
    end

    v.aleaf, v.gs
end


"""
    EmaxVars()

Varbles for Emax models

$(FIELDDOCTABLE)
"""
@MixinEnviroVars @MixinMaespaVars @MixinFvCBVars mutable struct EmaxVars{MLWP,EML,KT,PSIL}
    minleafwp::MLWP | 1e-4        | MPa                   | _
    emaxleaf::EML   | 400.0       | mmol*m^-2*s^-1        | _
    ktot::KT        | 2.0         | mmol*m^-2*s^-1*MPa^-1 | _
    psil::PSIL      | -0.111      | MPa                   | _
end


"""
    EmaxEnergyBalance(energy_balance_model, totsoilres, plantk)

Wrapper to MaespaEnergyBalance model, adding `totsoilres` and `plantk` parameters.

$(FIELDDOCTABLE)
"""
@columns struct EmaxEnergyBalance{EB,SR,PK} <: AbstractMaespaEnergyBalance
    energy_balance_model::EB | MaespaEnergyBalance(photosynthesis=FvCBPhotosynthesis(
                                                   stomatal_conductance=EmaxStomatalConductance(
                                                     soilmethod=EmaxSoilMethod()))) | _ | _ | _
    totsoilres::SR           | 0.5  | m^2*s^1*MPa^1*mmol^-1 | (0.0, 10.0) | _
    plantk::PK               | 3.0  | mmol*m^-2*s^-1*MPa^-1 | (0.0, 10.0) | _
end

energy_balance_model(m::EmaxEnergyBalance) = m.energy_balance_model
totsoilres(m::EmaxEnergyBalance) = m.totsoilres
plantk(m::EmaxEnergyBalance) = m.plantk

enbal!(v, m::EmaxEnergyBalance) = begin
    v.ktot = 1 / (totsoilres(m) + 1.0 / plantk(m))
    enbal!(v, energy_balance_model(m))

    #= MAESPA comments:
    Return re-calculated leaf water potential (using ET without boundary layer conductance).
    We use etest otherwise psil < psilmin quite frequently when soil is dry.
    This is difficult to interpret, especially because PHOTOSYN does not account
    for boundary layer conductance.
    =#
    etest = (v.vpd / v.pressure) * v.gsv
    v.psil = v.swp - (etest / v.ktot)

    @show etest
    @show v.swp
    @show v.pressure
    @show v.vpd
    @show v.gsv
    @show v.psil
    return
end

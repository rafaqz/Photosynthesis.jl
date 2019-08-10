using Photosynthesis, Unitful, Test

@testset "it all actually builds" begin
    BadgerCollatzCompensation()
    BernacchiCompensation()

    JarvisLinearCO2()
    JarvisNonlinearCO2()
    JarvisHyperbolicVPD()
    JarvisLohammerVPD()
    JarvisFractionDeficitVPD()
    JarvisLinearDeclineVPD()

    BallBerryGSsubModel()
    LeuningGSsubModel()
    MedlynGSsubModel()
    TuzetGSsubModel()

    Jmax()
    NoOptimumVcmax()
    OptimumVcmax()
    Flux()
    DukeFlux()
    RubiscoRegen()
    Respiration()

    YingPingRadiationConductance()
    BoundaryConductance()
    McNaughtonJarvisDecoupling()
    NoDecoupling()
    DeficitSoilData()
    ContentSoilData()
    SimulatedSoilData()
    LinearPotentialDependence()
    ZhouPotentialDependence()
    NoPotentialDependence()
    VolumetricSoilMethod()
    ConstantSoilMethod()
    DeficitSoilMethod()
    PotentialSoilMethod()
    EmaxSoilMethod()
    # TuzetSoilMethod()
    BallBerryStomatalConductance()
    EmaxStomatalConductance()
    JarvisStomatalConductance()
    BallBerryVars()
    EmaxVars()
    TuzetVars()
    JarvisVars()
    FvCBPhotosynthesis()
    FvCBEnergyBalance()
end

@testset "it all actuall runs" begin
    p = FvCBEnergyBalance()
    v = EmaxVars()


    # Biophysical #################################################################
    
    latent_heat_water_vapour(v.tair)
    arrhenius(42.75u"J/mol", 37830.0u"J/mol", v.tleaf, 300.0u"K")

    # Radiation conductance
    radiation_conductance(YingPingRadiationConductance(), v)

    # Shape
    shape_gs(HardMinimumGS(), v, BallBerryStomatalConductance()) 
    shape_gs(HyperbolicMinimumGS(), v, BallBerryStomatalConductance()) 

    # Boundary conductance
    grashof_number(v.tleaf, v.tair, BoundaryConductance().leafwidth)
    cmolar(v.pressure, v.tair)
    boundary_conductance_free(BoundaryConductance(), v)
    boundary_conductance_forced(BoundaryConductance(), v)

    # Decoupling
    decoupling(McNaughtonJarvisDecoupling(), v)
    decoupling(NoDecoupling(), v)

    # Evapotranspiration
    penman_monteith(v.pressure, v.slope, v.lhv, v.rnet, v.vpd, 1.0u"mol*m^-2*s^-1", 1.0u"mol*m^-2*s^-1") # TODO this seems weird
    evapotranspiration(PenmanMonteithEvapotranspiration(), v)


    # Core ########################################################################

    # Compensation
    co2_compensation_point(BadgerCollatzCompensation(), v.tleaf)
    co2_compensation_point(BernacchiCompensation(), v.tleaf)
    rubisco_compensation_point(BadgerCollatzCompensation(), v.tleaf)
    rubisco_compensation_point(BernacchiCompensation(), v.tleaf)

    # Flux
    flux(Flux(), v)
    flux(DukeFlux(), v)
    flux(Flux(), v)
    flux(DukeFlux(), v)

    # Respiration
    respiration(Respiration(), v.tleaf)

    # Rubisco regeneration
    rubisco_regeneration(RubiscoRegen(), v)

    # Stomatal conductance submodels
    gsdiva(BallBerryGSsubModel(), BallBerryVars())
    gsdiva(LeuningGSsubModel(), BallBerryVars())
    gsdiva(MedlynGSsubModel(), BallBerryVars())
    gsdiva(TuzetGSsubModel(), TuzetVars())

    # Stomatal conductance models
    stomatal_conductance!(BallBerryStomatalConductance(gs_submodel=BallBerryGSsubModel()), v)
    stomatal_conductance!(BallBerryStomatalConductance(gs_submodel=LeuningGSsubModel()), v)
    stomatal_conductance!(BallBerryStomatalConductance(gs_submodel=MedlynGSsubModel()), v)
    stomatal_conductance!(EmaxStomatalConductance(gs_submodel=BallBerryGSsubModel()), v)
    stomatal_conductance!(EmaxStomatalConductance(gs_submodel=LeuningGSsubModel()), v)
    stomatal_conductance!(EmaxStomatalConductance(gs_submodel=MedlynGSsubModel()), v)
    stomatal_conductance!(JarvisStomatalConductance(), JarvisVars())
    # stomatal_conductance!(TuzetStomatalConductance(), TuzetVars())


    # Formulations ################################################################

    v = BallBerryVars()
    ph = FvCBPhotosynthesis(stomatal_conductance=BallBerryStomatalConductance())
    p = FvCBEnergyBalance(photosynthesis=ph)
    enbal!(p, v)
    photosynthesis!(p.photosynthesis, v)
    v.aleaf
    v.tleaf
    v.gs
    v.cs
    enbal!(p, v)
    
    v = JarvisVars()
    ph = FvCBPhotosynthesis(stomatal_conductance=JarvisStomatalConductance())
    p = FvCBEnergyBalance(photosynthesis=ph)
    factor_conductance(ph.stomatal_conductance, v)
    enbal!(p, v)
    photosynthesis!(p.photosynthesis, v)
    v.aleaf
    v.tleaf
    v.gs
    v.cs
    enbal!(p, v)

    # ph = FvCBPhotosynthesis(stomatal_conductance=TuzetStomatalConductance())
    # p = TuzetEnergyBalance(photosynthesis=ph)
    # v = TuzetVars()
    # enbal!(p, v)
    # photosynthesis!(p.photosynthesis, v)
    # v.aleaf
    # v.tleaf
    # v.gs
    # v.cs
    # enbal!(p, v)

    v = BallBerryVars()
    ph = FvCBPhotosynthesis(stomatal_conductance=EmaxStomatalConductance())
    v = EmaxVars()
    enbal!(p, v)
    photosynthesis!(p.photosynthesis, v)
    v.aleaf
    v.tleaf
    v.gs
    v.cs
    enbal!(p, v)

end

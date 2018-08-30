
@testset "it all actually builds" begin
    BadgerCollatzCompensation()
    BernacchiCompensation()
    JarvisLinearCO2()
    JarvisNonlinearCO2()
    JarvisHyperbolicVPD()
    JarvisLohammerVPD()
    JarvisFractionDeficitVPD()
    JarvisLinearDeclineVPD()
    BallBerryStomatalConductance()
    LeuningStomatalConductance()
    MedlynStomatalConductance()
    ThreeParStomatalConductance()
    TuzetStomatalConductance()
    Jmax()
    NoOptimumVcmax()
    OptimumVcmax()
    VcJmax()
    DukeVcJmax()
    RubiscoRegen()
    Respiration()
    YingPingRadiationConductance()
    BoundaryConductance()
    McNaughtonJarvisDecoupling()
    NoDecoupling()
    DeficitSoilData()
    ContentSoilData()
    SimulatedSoilData()
    PotentialSoilData()
    NoSoilData()
    LinearPotentialDependence()
    ZhouPotentialDependence()
    NoPotentialDependence()
    VolumetricSoilMethod()
    ConstantSoilMethod()
    DeficitSoilMethod()
    PotentialSoilMethod()
    EmaxSoilMethod()
    TuzetSoilMethod()
    BallBerryModel()
    TuzetModel()
    EmaxModel()
    JarvisModel()
    FvCBPhoto()
    PhotoVars()
    EmaxVars()
    TuzetVars()
    BallBerryVars()
    JarvisVars()
end

@testset "it all actuall runs" begin
    p = EnergyBalance()
    v = EmaxVars()

    penman_monteith(v.pressure, v.slope, v.lhv, v.rnet, v.vpd, 1.0u"mol*m^-2*s^-1", 1.0u"mol*m^-2*s^-1") # TODO this seems weird

    grashof_number(v.tleaf, v.tair, p.boundary_conductance.leafwidth)
    cmolar(v.pressure, v.tair)
    v.gradn = radiation_conductance(p.radiation_conductance, v, p)
    boundary_conductance_free(p.boundary_conductance, v, p)
    boundary_conductance_forced(p.boundary_conductance, v, p)
    latent_heat_water_vapour(v.tair)
    arrhenius(42.75u"J/mol", 37830.0u"J/mol", v.tleaf, 25.0u"°C")

    co2_compensation_point(BadgerCollatzCompensation(), v, p)
    co2_compensation_point(BernacchiCompensation(), v, p)
    rubisco_compensation_point(BadgerCollatzCompensation(), v, p)
    rubisco_compensation_point(BernacchiCompensation(), v, p)

    max_electron_transport_rate(VcJmax(), v, p)
    max_electron_transport_rate(DukeVcJmax(), v, p)
    max_rubisco_activity(VcJmax(),v, p)
    max_rubisco_activity(DukeVcJmax(),v, p)
    v.cs
    respiration(Respiration(), v, p)

    gsdiva(BallBerryStomatalConductance(), BallBerryVars(), p.photo)
    gsdiva(LeuningStomatalConductance(), BallBerryVars(), p.photo)
    gsdiva(MedlynStomatalConductance(), BallBerryVars(), p.photo)
    gsdiva(ThreeParStomatalConductance(), BallBerryVars(), p.photo)
    gsdiva(TuzetStomatalConductance(), TuzetVars(), FvCBPhoto(model=TuzetModel()))

    stomatal_conductance!(BallBerryVars(), BallBerryModel(gsmodel=BallBerryStomatalConductance()), FvCBPhoto())
    stomatal_conductance!(BallBerryVars(), BallBerryModel(gsmodel=LeuningStomatalConductance()), FvCBPhoto())
    stomatal_conductance!(BallBerryVars(), BallBerryModel(gsmodel=MedlynStomatalConductance()), FvCBPhoto())
    stomatal_conductance!(BallBerryVars(), BallBerryModel(gsmodel=ThreeParStomatalConductance()), FvCBPhoto())
    stomatal_conductance!(EmaxVars(), EmaxModel(gsmodel=BallBerryStomatalConductance()), FvCBPhoto())
    stomatal_conductance!(EmaxVars(), EmaxModel(gsmodel=LeuningStomatalConductance()), FvCBPhoto())
    stomatal_conductance!(EmaxVars(), EmaxModel(gsmodel=MedlynStomatalConductance()), FvCBPhoto())
    stomatal_conductance!(EmaxVars(), EmaxModel(gsmodel=ThreeParStomatalConductance()), FvCBPhoto())
    stomatal_conductance!(JarvisVars(), JarvisModel(), FvCBPhoto())
    stomatal_conductance!(TuzetVars(), TuzetModel(), FvCBPhoto())

    ballberryplant = EnergyBalance(photo=FvCBPhoto(model=BallBerryModel()))
    tuzetplant = EnergyBalance(photo=FvCBPhoto(model=TuzetModel()))
    emaxplant = EnergyBalance(photo=FvCBPhoto(model=EmaxModel()))
    jarvisplant = EnergyBalance(photo=FvCBPhoto(model=JarvisModel()))

    factor_conductance(JarvisModel(), v, jarvisplant)
    v = BallBerryVars()
    p = ballberryplant
    photosynthesis!(v, p.photo)
    v.aleaf
    v.tleaf
    v.gs
    v.cs
    phototranspiration!(v, p)

    p = tuzetplant
    v = TuzetVars()
    photosynthesis!(v, p.photo)
    v.aleaf
    v.tleaf
    v.gs
    v.cs
    phototranspiration!(v, p)

    p = emaxplant
    v = EmaxVars()
    photosynthesis!(v, p.photo)
    v.aleaf
    v.tleaf
    v.gs
    v.cs
    phototranspiration!(v, p)

    run_phototrans!(EmaxVars(), emaxplant)
    run_phototrans!(TuzetVars(), tuzetplant)
end

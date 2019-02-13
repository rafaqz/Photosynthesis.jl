using Photosynthesis,
      Unitful, 
      Test

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
    PhotoVars()
    EmaxVars()
    TuzetVars()
    BallBerryVars()
    JarvisVars()
    FvCBEnergyBalance()
end

@testset "it all actuall runs" begin
    global p = FvCBEnergyBalance()
    global v = EmaxVars()

    penman_monteith(v.pressure, v.slope, v.lhv, v.rnet, v.vpd, 1.0u"mol*m^-2*s^-1", 1.0u"mol*m^-2*s^-1") # TODO this seems weird

    grashof_number(v.tleaf, v.tair, p.boundary_conductance.leafwidth)
    cmolar(v.pressure, v.tair)
    v.gradn = radiation_conductance(p.radiation_conductance, v)
    boundary_conductance_free(p.boundary_conductance, v)
    boundary_conductance_forced(p.boundary_conductance, v)
    latent_heat_water_vapour(v.tair)
    arrhenius(42.75u"J/mol", 37830.0u"J/mol", v.tleaf, 300.0u"K")

    co2_compensation_point(BadgerCollatzCompensation(), v)
    co2_compensation_point(BernacchiCompensation(), v)
    rubisco_compensation_point(BadgerCollatzCompensation(), v)
    rubisco_compensation_point(BernacchiCompensation(), v)

    max_electron_transport_rate(VcJmax(), v)
    max_electron_transport_rate(DukeVcJmax(), v)
    max_rubisco_activity(VcJmax(),v)
    max_rubisco_activity(DukeVcJmax(),v)

    respiration(Respiration(), v)

    rubisco_regeneration(RubiscoRegen(), v)

    gsdiva(BallBerryStomatalConductance(), BallBerryVars())
    gsdiva(LeuningStomatalConductance(), BallBerryVars())
    gsdiva(MedlynStomatalConductance(), BallBerryVars())
    gsdiva(ThreeParStomatalConductance(), BallBerryVars())
    gsdiva(TuzetStomatalConductance(), TuzetVars())

    stomatal_conductance!(BallBerryModel(gsmodel=BallBerryStomatalConductance()), BallBerryVars())
    stomatal_conductance!(BallBerryModel(gsmodel=LeuningStomatalConductance()), BallBerryVars())
    stomatal_conductance!(BallBerryModel(gsmodel=MedlynStomatalConductance()), BallBerryVars())
    stomatal_conductance!(BallBerryModel(gsmodel=ThreeParStomatalConductance()), BallBerryVars())
    stomatal_conductance!(EmaxModel(gsmodel=BallBerryStomatalConductance()), EmaxVars())
    stomatal_conductance!(EmaxModel(gsmodel=LeuningStomatalConductance()), EmaxVars())
    stomatal_conductance!(EmaxModel(gsmodel=MedlynStomatalConductance()), EmaxVars())
    stomatal_conductance!(EmaxModel(gsmodel=ThreeParStomatalConductance()), EmaxVars())
    stomatal_conductance!(JarvisModel(), JarvisVars())
    stomatal_conductance!(TuzetModel(), TuzetVars())

    global v = JarvisVars()
    global p = FvCBEnergyBalance(photo=JarvisModel())
    photosynthesis!(p.photo, v)
    v.aleaf
    v.tleaf
    v.gs
    v.cs
    enbal!(p, v)
    run_enbal!(p, v)
    factor_conductance(JarvisModel(), v)

    global v = BallBerryVars()
    global p = FvCBEnergyBalance(photo=BallBerryModel())
    photosynthesis!(p.photo, v)
    v.aleaf
    v.tleaf
    v.gs
    v.cs
    enbal!(p, v)

    global p = FvCBEnergyBalance(photo=TuzetModel())
    global v = TuzetVars()
    photosynthesis!(p.photo, v)
    v.aleaf
    v.tleaf
    v.gs
    v.cs
    enbal!(p, v)
    run_enbal!(p, v)

    global p = FvCBEnergyBalance(photo=EmaxModel())
    global v = EmaxVars()
    photosynthesis!(p.photo, v)
    v.aleaf
    v.tleaf
    v.gs
    v.cs
    enbal!(p, v)
    run_enbal!(p, v)
end

nothing

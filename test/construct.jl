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
    BallBerryPhotosynthesis()
    TuzetPhotosynthesis()
    EmaxPhotosynthesis()
    JarvisPhotosynthesis()
    BallBerryVars()
    EmaxVars()
    TuzetVars()
    JarvisVars()
    FvCBEnergyBalance()
end

@testset "it all actuall runs" begin
    p = FvCBEnergyBalance()
    v = EmaxVars()
    enbal!(p, v)


    # Biophysical 
    latent_heat_water_vapour(v.tair)
    arrhenius(42.75u"J/mol", 37830.0u"J/mol", v.tleaf, 300.0u"K")

    # Radiation conductance
    radiation_conductance(YingPingRadiationConductance(), v)

    # Shape
    shape_gs(HardMinimumGS(), v, BallBerryPhotosynthesis()) 
    shape_gs(HyperbolicMinimumGS(), v, BallBerryPhotosynthesis()) 

    # Boundary conductance
    grashof_number(v.tleaf, v.tair, BoundaryConductance().leafwidth)
    cmolar(v.pressure, v.tair)
    boundary_conductance_free(BoundaryConductance(), v)
    boundary_conductance_forced(BoundaryConductance(), v)

    # Decoupling
    calc_decoupling(McNaughtonJarvisDecoupling(), v)
    calc_decoupling(NoDecoupling(), v)

    # Evapotranspiration
    penman_monteith(v.pressure, v.slope, v.lhv, v.rnet, v.vpd, 1.0u"mol*m^-2*s^-1", 1.0u"mol*m^-2*s^-1") # TODO this seems weird
    calc_evapotranspiration(PenmanMonteithEvapotranspiration(), v)


    # Core

    # Compensation
    co2_compensation_point(BadgerCollatzCompensation(), v.tleaf)
    co2_compensation_point(BernacchiCompensation(), v.tleaf)
    rubisco_compensation_point(BadgerCollatzCompensation(), v.tleaf)
    rubisco_compensation_point(BernacchiCompensation(), v.tleaf)

    # Flux
    max_electron_transport_rate(VcJmax(), v.tleaf)
    max_electron_transport_rate(DukeVcJmax(), v.tleaf)
    max_rubisco_activity(VcJmax(), v.tleaf)
    max_rubisco_activity(DukeVcJmax(), v.tleaf)

    # Respiration
    respiration(Respiration(), v.tleaf)

    # Rubisco regeneration
    rubisco_regeneration(RubiscoRegen(), v)

    # Stomatal conductance submodels
    gsdiva(BallBerryStomatalConductance(), BallBerryVars())
    gsdiva(LeuningStomatalConductance(), BallBerryVars())
    gsdiva(MedlynStomatalConductance(), BallBerryVars())
    gsdiva(ThreeParStomatalConductance(), BallBerryVars())
    gsdiva(TuzetStomatalConductance(), TuzetVars())

    # Stomatal conductance models
    stomatal_conductance!(BallBerryPhotosynthesis(gsmodel=BallBerryStomatalConductance()), v)
    stomatal_conductance!(BallBerryPhotosynthesis(gsmodel=LeuningStomatalConductance()), v)
    stomatal_conductance!(BallBerryPhotosynthesis(gsmodel=MedlynStomatalConductance()), v)
    stomatal_conductance!(BallBerryPhotosynthesis(gsmodel=ThreeParStomatalConductance()), v)
    stomatal_conductance!(EmaxPhotosynthesis(gsmodel=BallBerryStomatalConductance()), v)
    stomatal_conductance!(EmaxPhotosynthesis(gsmodel=LeuningStomatalConductance()), v)
    stomatal_conductance!(EmaxPhotosynthesis(gsmodel=MedlynStomatalConductance()), v)
    stomatal_conductance!(EmaxPhotosynthesis(gsmodel=ThreeParStomatalConductance()), v)
    stomatal_conductance!(JarvisPhotosynthesis(), JarvisVars())
    stomatal_conductance!(TuzetPhotosynthesis(), TuzetVars())


    # Formulations
    
    v = JarvisVars()
    p = FvCBEnergyBalance(photosynthesis=JarvisPhotosynthesis())
    factor_conductance(JarvisPhotosynthesis(), v)
    enbal!(p, v)
    photosynthesis!(p.photosynthesis, v)
    v.aleaf
    v.tleaf
    v.gs
    v.cs
    enbal!(p, v)

    v = BallBerryVars()
    p = FvCBEnergyBalance(photosynthesis=BallBerryPhotosynthesis())
    enbal!(p, v)
    photosynthesis!(p.photosynthesis, v)
    v.aleaf
    v.tleaf
    v.gs
    v.cs
    enbal!(p, v)

    p = FvCBEnergyBalance(photosynthesis=TuzetPhotosynthesis())
    v = TuzetVars()
    enbal!(p, v)
    photosynthesis!(p.photosynthesis, v)
    v.aleaf
    v.tleaf
    v.gs
    v.cs
    enbal!(p, v)

    p = FvCBEnergyBalance(photosynthesis=EmaxPhotosynthesis())
    v = EmaxVars()
    enbal!(p, v)
    photosynthesis!(p.photosynthesis, v)
    v.aleaf
    v.tleaf
    v.gs
    v.cs
    enbal!(p, v)

end

nothing

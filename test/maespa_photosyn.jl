using Unitful: °C, K
include("shared.jl")

function fortran_photosyn(p, v::V, ieco, ismaespabool, modelgs, wsoilmethod, soildata, vfun) where V
    ismaespa = unsigned(0) + ismaespabool
    # Use current vars or dummy if they arent needed
    tzvars = V <: TuzetVars ? v : TuzetVars()
    emaxvars = V <: EmaxVars ? v : EmaxVars()

    ph = p.photosynthesis_model
    sc = ph.stomatal_conductance_model
    gs = sc.gs_submodel
    vcj = Photosynthesis.fluxparams(ph.flux_model)
    vc = vcj.vcmaxformulation
    j = vcj.jmaxformulation
    gk = typeof(sc.gs_submodel) <: MedlynStomatalConductanceSubModel ? sc.gs_submodel.gk : 0.0
    d0l = typeof(sc.gs_submodel) <: LeuningStomatalConductanceSubModel ? sc.gs_submodel.d0l.val : 0.0
    par =      v.par.val
    tleaf =    Float32[ustrip(v.tleaf |> °C)]
    tmove =    AcclimatizedRespiration().tmove.val
    cs =       v.cs.val
    ca =       v.ca.val
    rh =       v.rh
    rnet =     v.rnet.val
    tair =     ustrip(v.tair |> °C)
    wind =     v.windspeed.val
    vpd =      v.vpd.val
    press =    v.pressure.val
    vmfd =     JarvisStomatalConductance().vmfd.val
    jmax25 =   vcj.jmaxformulation.jmax25.val
    eavj =     j.eavj.val
    edvj =     j.edvj.val
    delsj =    j.delsj.val
    vcmax25 =  vc.vcmax25.val
    eavc =     vc.eavc.val
    if typeof(vc) <: OptimumVcmax
        edvc = vc.edvc.val
        delsc = vc.delsc.val
    else
        edvc = 0.0
        delsc = 0.0
    end
    tvjup =    ustrip(DukeFlux().tvjup |> °C)
    tvjdn =    ustrip(DukeFlux().tvjdn |> °C)
    theta =    ph.rubisco_regen_model.theta
    ajq =      ph.rubisco_regen_model.ajq
    rd0 =      ph.respiration_model.rd0.val
    q10f =     ph.respiration_model.q10f.val
    k10f =     AcclimatizedRespiration().k10f.val
    tref =     ustrip(ph.respiration_model.tref |> °C)
    rtemp =    ustrip(ph.respiration_model.tref |> °C)
    dayresp =  ph.respiration_model.dayresp
    tbelow =   ustrip(ph.respiration_model.tbelow |> °C)
    gsref =    JarvisStomatalConductance().gsref.val
    gsmin =    JarvisStomatalConductance().gsmin.val
    i0 =       JarvisLight().i0.val
    d0 =       JarvisLinearDeclineVPD().d0.val
    vk1 =      JarvisHyperbolicVPD().vk1
    vk2 =      JarvisHyperbolicVPD().vk2
    vpd1 =     JarvisLohammerVPD().vpd1.val
    vpd2 =     JarvisLohammerVPD().vpd2.val
    vmfd0 =    JarvisFractionDeficitVPD().vmfd0.val
    gsja =     JarvisLinearCO2().gsja.val
    gsjb =     JarvisNonlinearCO2().gsjb.val
    t0 =       ustrip(JarvisTemp1().t0 |> °C)
    tref =     ustrip(JarvisTemp1().tref |> °C)
    tmax =     ustrip(JarvisTemp1().tmax |> °C)
    soilmoisture = typeof(sc.soil_model) <: PotentialSoilMethod ? v.swp.val : v.soilmoist  
    emaxleaf = emaxvars.emaxleaf.val
    plantk =   EmaxEnergyBalance().plantk.val
    totsoilres = EmaxEnergyBalance().totsoilres.val

    smd1 =     DeficitSoilData().smd1
    smd2 =     DeficitSoilData().smd2
    wc1 =      VolumetricSoilMethod().wc1
    wc2 =      VolumetricSoilMethod().wc2
    swpexp =   PotentialSoilMethod().swpexp.val
    fsoil =    Float32[v.fsoil]
    g0 =       sc.g0.val
    gamma =    sc.gs_submodel.gamma.val
    vpdmin =   MedlynStomatalConductanceSubModel().vpdmin.val
    g1 =       sc.gs_submodel.g1
    gk =       MedlynStomatalConductanceSubModel().gk
    gs =       Float32[v.gs.val]
    aleaf =    Float32[v.aleaf.val]
    rd =       Float32[v.rd.val]
    minleafwp = emaxvars.minleafwp.val
    ktot =     emaxvars.ktot.val
    weightedswp = emaxvars.swp.val
    vpara =    LinearPotentialDependence().vpara.val
    vparb =    LinearPotentialDependence().vparb.val
    vparc =    0.0 # unused
    fheat =    0.0 # unused
    etest =    0.0 # unused
    gbh =      v.gbh.val
    sf =       tzvars.sf.val
    psiv =     tzvars.psiv.val
    hmshape =  HyperbolicMinimum().hmshape
    psilin =   tzvars.psilin.val
    psil =     Float32[emaxvars.psil.val]
    ci =       Float32[v.ci.val]

    photo = Libdl.dlsym(photosynlib, :photosyn_)
    ccall(photo, Nothing, (
    Ref{Float32},
    Ref{Float32},
    Ref{Float32},
    Ref{Float32},
    Ref{Float32}, Ref{Float32},
    Ref{Float32},
    Ref{Float32},
    Ref{Int32},
    Ref{Float32},
    Ref{Float32},
    Ref{Float32},
    Ref{Float32},
    Ref{Float32},
    Ref{Float32},
    Ref{Float32},
    Ref{Float32},
    Ref{Float32},
    Ref{Float32},
    Ref{Float32},
    Ref{Float32},
    Ref{Float32},
    Ref{Float32},
    Ref{Float32},
    Ref{Float32},
    Ref{Float32},
    Ref{Int32},
    Ref{Float32},
    Ref{Float32},
    Ref{Float32},
    Ref{Float32},
    Ref{Float32},
    Ref{Float32},
    Ref{Float32},
    Ref{Float32},
    Ref{Float32},
    Ref{Float32},
    Ref{Float32},
    Ref{Float32},
    Ref{Float32},
    Ref{Float32},
    Ref{Int32},
    Ref{Float32},
    Ref{Float32},
    Ref{Float32},
    Ref{Float32},
    Ref{Float32},
    Ref{Float32},
    Ref{Int32},
    Ref{Float32},
    Ref{Float32},
    Ref{Float32},
    Ref{Float32},
    Ref{Float32},
    Ref{Float32},
    Ref{Float32},
    Ref{Float32},
    Ref{Float32},
    Ref{Float32},
    Ref{Float32},
    Ref{Float32},
    Ref{Float32},
    Ref{Float32},
    Ref{Float32},
    Ref{Float32},
    Ref{Float32},
    Ref{Float32},
    Ref{Float32},
    Ref{Float32},
    Ref{Float32},
    Ref{Float32},
    Ref{Float32},
    Ref{Float32},
    Ref{Int32}),
    par,              
    tleaf,            
    tmove,            
    cs,               
    rh,               
    vpd,              
    vmfd,             
    jmax25,           
    ieco,             
    eavj,             
    edvj,             
    delsj,            
    vcmax25,          
    eavc,             
    edvc,             
    delsc,            
    tvjup,            
    tvjdn,            
    theta,            
    ajq,              
    rd0,              
    q10f,             
    k10f,             
    rtemp,            
    dayresp,          
    tbelow,           
    modelgs,          
    gsref,            
    gsmin,            
    i0,               
    d0,               
    vk1,              
    vk2,              
    vpd1,             
    vpd2,             
    vmfd0,            
    gsja,             
    gsjb,             
    t0,               
    tref,             
    tmax,             
    wsoilmethod,      
    soilmoisture,     
    emaxleaf,         
    smd1,             
    smd2,             
    wc1,              
    wc2,              
    soildata,         
    swpexp,           
    fsoil,            
    g0,               
    d0l,              
    gamma,            
    vpdmin,           
    g1,               
    gk,               
    gs,               
    aleaf,            
    rd,               
    minleafwp,        
    ktot,             
    weightedswp,      
    vpara,            
    vparb,            
    vparc,            
    vfun,             
    sf,               
    psiv,             
    hmshape,          
    psilin,           
    psil,             
    ci,               
    ismaespa)         

    (rd[1], fsoil[1], aleaf[1], gs[1], ci[1], ismaespa)
end

@testset "Ball Berry" begin
    p = FvCBEnergyBalance(photosynthesis_model=FvCBPhotosynthesis(
                               stomatal_conductance_model=BallBerryStomatalConductance(
                                    gs_submodel=BallBerryStomatalConductanceSubModel(),
                                    soil_model=NoSoilMethod()),
                               flux_model=DukeFlux(),
                               compensation=BadgerCollatzCompensation()))
    v = BallBerryVars()
    ieco = BADGERCOLLATZ
    modelgs = BALLBERRY_GS
    wsoilmethod = 0
    soildata = SOILDATA_NONE
    ismaespa = false
    vfun = 1
    # Prime variables
    enbal!(v, p)

    # Run fortran
    (rd, fsoil, aleaf, gs, ci) =
        fortran_photosyn(p, v, ieco, ismaespa, modelgs, wsoilmethod, soildata, vfun)
    # Run julia
    photosynthesis!(v, p.photosynthesis_model)

    v.rd, v.fsoil, v.aleaf, v.gs, v.ci
    v.tleaf |> °C

    @test rd ≈ v.rd.val
    @test fsoil ≈ v.fsoil rtol=1e-7
    @test gs ≈ v.gs.val rtol=1e-6
    @test aleaf ≈ v.aleaf.val# rtol=1e-6
end

# using DataFrames, Unitful, Flatten
# DataFrame(name = [FieldMetadata.fieldnames(BallBerryVars)...], 
#           val = ulflatten(Vector, v),
#           default = [FieldMetadata.default(BallBerryVars)...], 
#           units = [FieldMetadata.units(BallBerryVars)...])


@testset "Ball Berry water potential" begin
    p = FvCBEnergyBalance(
         photosynthesis_model=FvCBPhotosynthesis(
             stomatal_conductance_model=BallBerryStomatalConductance(
                 gs_submodel=BallBerryStomatalConductanceSubModel(),
                 soil_model=PotentialSoilMethod(),
             ),
             flux_model=DukeFlux(),
             compensation_model=BadgerCollatzCompensation(),
         )
    )
    v = BallBerryVars()
    ieco = BADGERCOLLATZ
    modelgs = BALLBERRY_GS
    wsoilmethod = 2
    soildata = SOILDATA_POTENTIAL
    ismaespa = false
    vfun = 1
    # Prime variables
    enbal!(v, p)

    # Run fortran
    (rd, fsoil, aleaf, gs, ci) =
        fortran_photosyn(p, v, ieco, ismaespa, modelgs, wsoilmethod, soildata, vfun)
    # Run julia
    photosynthesis!(v, p.photosynthesis_model)

    v.rd, v.fsoil, v.aleaf, v.gs, v.ci
    v.tleaf |> °C

    @test fsoil ≈ v.fsoil
    @test rd ≈ v.rd.val
    @test aleaf ≈ v.aleaf.val# rtol=1e-6
    @test gs ≈ v.gs.val rtol=1e-6
end

@testset "Leuning" begin
    p = FvCBEnergyBalance(
         photosynthesis_model=FvCBPhotosynthesis(
             stomatal_conductance_model=BallBerryStomatalConductance(
                 gs_submodel=LeuningStomatalConductanceSubModel(),
                 soil_model=NoSoilMethod(),
             ),
             flux_model=DukeFlux(),
             compensation_model=BadgerCollatzCompensation()
        )
    )
    v = BallBerryVars()
    ieco = BADGERCOLLATZ
    modelgs = BALLBERRY_GS
    wsoilmethod = 0
    soildata = SOILDATA_NONE
    ismaespa = false
    vfun = 1
    # Prime variables
    enbal!(v, p)

    # Run fortran
    (rd, fsoil, aleaf, gs, ci) =
        fortran_photosyn(p, v, ieco, ismaespa, modelgs, wsoilmethod, soildata, vfun)
    # Run julia
    photosynthesis!(v, p.photosynthesis_model)

    v.rd, v.fsoil, v.aleaf, v.gs, v.ci
    v.tleaf |> °C

    @test rd ≈ v.rd.val
    @test fsoil ≈ v.fsoil rtol=1e-7
    @test_broken gs ≈ v.gs.val rtol=1e-6
    @test_broken aleaf ≈ v.aleaf.val# rtol=1e-6
end

@testset "Medlyn" begin
    p = FvCBEnergyBalance(
        photosynthesis_model=FvCBPhotosynthesis(
            stomatal_conductance=BallBerryStomatalConductance(
                gs_submodel=MedlynStomatalConductanceSubModel(),
                soil_model=NoSoilMethod()),
            flux_model=DukeFlux(),
            compensation_model=BadgerCollatzCompensation()
        )
    )
    v = BallBerryVars()
    ieco = BADGERCOLLATZ
    modelgs = MEDLYN_GS
    wsoilmethod = 0
    soildata = SOILDATA_NONE
    ismaespa = false
    vfun = 1
    # Prime variables
    enbal!(v, p)

    # Run fortran
    (rd, fsoil, aleaf, gs, ci) =
        fortran_photosyn(p, v, ieco, ismaespa, modelgs, wsoilmethod, soildata, vfun)
    # Run julia
    photosynthesis!(v, p.photosynthesis_model)

    v.rd, v.fsoil, v.aleaf, v.gs, v.ci
    v.tleaf |> °C

    @test rd ≈ v.rd.val
    @test fsoil ≈ v.fsoil rtol=1e-7
    @test_broken gs ≈ v.gs.val rtol=1e-6
    @test_broken aleaf ≈ v.aleaf.val# rtol=1e-6
end

@testset "Emax" begin
    p1 = FvCBEnergyBalance(
        photosynthesis_model=FvCBPhotosynthesis(
            stomatal_conductance_model=EmaxStomatalConductance(
                gs_submodel=BallBerryStomatalConductanceSubModel(),
                soil_model=EmaxSoilMethod()
            ),
            flux_model=DukeFlux(),
        compensation_model=BadgerCollatzCompensation())
    )
    p = EmaxEnergyBalance(energy_balance_model=p1)
    v = EmaxVars()
    ieco = BADGERCOLLATZ
    modelgs = BALLBERRY_GS
    wsoilmethod = 1
    soildata = SOILDATA_POTENTIAL
    ismaespa = true
    vfun = 1

    # Prime variables
    enbal!(v, p)
    # Run fortran
    (rd, fsoil, aleaf, gs, ci) =
        fortran_photosyn(p.energy_balance_model, v, ieco, ismaespa, modelgs, wsoilmethod, soildata, vfun)
    # Run julia
    photosynthesis!(v, p.energy_balance_model.photosynthesis_model)
    v.rd, v.fsoil, v.aleaf, v.gs, v.ci

    @test fsoil ≈ v.fsoil rtol=1e-7
    @test rd ≈ v.rd.val
    @test aleaf ≈ v.aleaf.val# rtol=1e-6
    @test gs ≈ v.gs.val rtol=1e-7
end

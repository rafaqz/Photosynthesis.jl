using Photosynthesis

include(joinpath(dirname(pathof(Photosynthesis)), "../test/shared.jl"))

function fortran_photosyn(p, v, vfun=1)

    # Use current vars or dummy if they arent needed
    tzvars = v isa TuzetVars ? v : TuzetVars()
    emaxvars = v isa EmaxVars ? v : EmaxVars()

    ph = p.photosynthesis_model
    sc = ph.stomatal_conductance_model
    gs = sc.gs_submodel
    vcj = flux_model(ph.flux_model)
    vc = vcj.vcmaxformulation
    j = vcj.jmaxformulation
    sm = sc.soil_model

    ieco = ph.compensation_model isa BadgerCollatzCompensation ? BADGERCOLLATZ : BERNACCI
    if gs isa BallBerryStomatalConductanceSubModel 
        vpdmin = 0.0
        gk = 0.0
        modelgs = BALLBERRY_GS
    elseif gs isa LeuningStomatalConductanceSubModel 
        vpdmin = 0.0
        gk = 0.0
        modelgs = LEUNING_GS
    elseif gs isa MedlynStomatalConductanceSubModel 
        vpdmin =   ustrip(u"kPa", MedlynStomatalConductanceSubModel().vpdmin)
        gk =       MedlynStomatalConductanceSubModel().gk
        modelgs = MEDLYN_GS
    end
    if sm isa NoSoilMethod
        wc1 = 0.0
        wc2 = 0.0
        swpexp = 0.0
        soilmoisture = v.soilmoist  
        wsoilmethod = NOSOILMETHOD
        soildata = SOILDATA_NONE
    elseif sm isa PotentialSoilMethod
        wc1 = 0.0
        wc2 = 0.0
        soilmoisture = ustrip(u"Pa", v.swp)
        swpexp = ustrip(u"Pa^-1", sm.swpexp)
        wsoilmethod = POTENTIALSOILMETHOD
        soildata = SOILDATA_POTENTIAL
    elseif sm isa VolumetricSoilMethod
        wc1 = sm.wc1
        wc2 = sm.wc2
        soilmoisture = v.soilmoist  
        swpexp = 0.0
        wsoilmethod = VOLUMETRICSOILMETHOD
        soildata = SOILDATA_CONTENT
    elseif sm isa EmaxSoilMethod
        soilmoisture = v.soilmoist  
        swpexp = 0.0
        wc1 = 0.0
        wc2 = 0.0
        wsoilmethod = EMAXSOILMETHOD
        soildata = SOILDATA_NONE
    end
    ismaespabool = sc isa EmaxStomatalConductance ? true : false
    ismaespa = unsigned(0) + ismaespabool

    gk = typeof(sc.gs_submodel) <: MedlynStomatalConductanceSubModel ? sc.gs_submodel.gk : 0.0
    D0 = typeof(sc.gs_submodel) <: LeuningStomatalConductanceSubModel ? ustrip(u"Pa", sc.gs_submodel.D0) : 0.0
    par =      ustrip(u"μmol*m^-2*s^-1", v.par)
    tleaf =    Float32[ustrip(°C, v.tleaf)]
    cs =       ustrip(u"μmol/mol", v.cs)
    ca =       ustrip(u"μmol/mol", v.ca)
    rh =       v.rh
    rnet =     ustrip(u"J*m^-2*s^-1", v.rnet)
    tair =     ustrip(°C, v.tair)
    wind =     ustrip(u"m*s^-1", v.windspeed)
    vpd =      ustrip(u"Pa", v.vpd)
    press =    ustrip(u"Pa", v.pressure)
    vmfd =     JarvisStomatalConductance().vmfd.val
    jmax25 =   ustrip(u"μmol*m^-2*s^-1", vcj.jmaxformulation.jmax25)
    eavj =     ustrip(u"J*mol^-1", j.eavj)
    edvj =     ustrip(u"J*mol^-1", j.edvj)
    delsj =    ustrip(u"J*K^-1*mol^-1", j.delsj)
    vcmax25 =  ustrip(u"μmol*m^-2*s^-1", vc.vcmax25)
    eavc =     ustrip(u"J*mol^-1", vc.eavc)
    if typeof(vc) <: OptimumVcmax
        edvc = ustrip(u"J*mol^-1", vc.edvc)
        delsc = ustrip(u"J*K^-1*mol^-1", vc.delsc)
    else
        edvc = 0.0
        delsc = 0.0
    end
    tvjup =    ustrip(°C, DukeFlux().tvjup)
    tvjdn =    ustrip(°C, DukeFlux().tvjdn)
    theta =    ph.rubisco_regen_model.theta
    ajq =      ph.rubisco_regen_model.ajq
    rd0 =      ustrip(u"μmol*m^-2*s^-1", ph.respiration_model.rd0)
    q10f =     ustrip(u"K^-1", ph.respiration_model.q10f)
    if ph.respiration_model isa AcclimatizedRespiration
        k10f = ustrip(u"K^-1", AcclimatizedRespiration().k10f)
        tmove = ustrip(K, AcclimatizedRespiration().tmove)
    else
        k10f = 0.0
        tmove = 0.0
    end
    tref =     ustrip(°C, ph.respiration_model.tref)
    rtemp =    ustrip(°C, ph.respiration_model.tref)
    dayresp =  ph.respiration_model.dayresp
    tbelow =   ustrip(°C, ph.respiration_model.tbelow)
    gsref =    JarvisStomatalConductance().gsref.val
    gsmin =    JarvisStomatalConductance().g0.val
    i0 =       JarvisLight().i0.val
    d0 =       JarvisLinearDeclineVPD().d0.val
    vk1 =      JarvisHyperbolicVPD().vk1
    vk2 =      JarvisHyperbolicVPD().vk2
    vpd1 =     JarvisLohammerVPD().vpd1.val
    vpd2 =     JarvisLohammerVPD().vpd2.val
    vmfd0 =    JarvisFractionDeficitVPD().vmfd0.val
    gsja =     JarvisLinearCO2().gsja.val
    gsjb =     JarvisNonlinearCO2().gsjb.val
    t0 =       ustrip(°C, JarvisTemp1().t0)
    tmax =     ustrip(°C, JarvisTemp1().tmax)
    emaxleaf = emaxvars.emaxleaf.val
    plantk =   EmaxEnergyBalance().plantk.val
    totsoilres = EmaxEnergyBalance().totsoilres.val

    smd1 =     DeficitSoilData().smd1
    smd2 =     DeficitSoilData().smd2
    fsoil =    Float32[v.fsoil]
    g0 =       ustrip(u"mol*m^-2*s^-1", sc.g0)
    gamma =    sc.gs_submodel.gamma
    g1 =       sc.gs_submodel.g1
    gs =       Float32[ustrip(u"mol*m^-2*s^-1", v.gs)]
    aleaf =    Float32[ustrip(u"mol*m^-2*s^-1", v.aleaf)]
    rd =       Float32[ustrip(u"mol*m^-2*s^-1", v.rd)]
    minleafwp = emaxvars.minleafwp.val
    ktot =     emaxvars.ktot.val
    weightedswp = emaxvars.swp.val
    vpara =    ustrip(u"Pa", LinearPotentialDependence().vpara)
    vparb =    ustrip(u"Pa", LinearPotentialDependence().vparb)
    vparc =    0.0 # unused
    fheat =    0.0 # unused
    etest =    0.0 # unused
    gbh =      ustrip(u"mol*m^-2*s^-1", v.gbh)
    sf =       tzvars.sf.val
    psiv =     tzvars.psiv.val
    hmshape =  HyperbolicMinimum().hmshape
    psilin =   tzvars.psilin.val
    psil =     Float32[emaxvars.psil.val]
    ci =       Float32[ustrip(u"μmol*mol^-1", v.ci)]

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
    D0,              
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

    (rd[1], fsoil[1], aleaf[1], gs[1], ci[1], tleaf[1])
end

@testset "Ball Berry" begin
    p = FvCBEnergyBalance(
         photosynthesis_model=FvCBPhotosynthesis(
             stomatal_conductance_model=BallBerryStomatalConductance(
                 gs_submodel=BallBerryStomatalConductanceSubModel(),
                 soil_model=NoSoilMethod(),
             ),
             flux_model=DukeFlux(),
             compensation_model=BadgerCollatzCompensation(),
             respiration_model=Respiration(),
         )
    )

    v = BallBerryVars()
    # Prime variables
    enbal!(v, p)

    # Run fortran
    (rd, fsoil, aleaf, gs, ci, tleaf) =
        fortran_photosyn(p, v)
    # Run julia
    photosynthesis!(v, p.photosynthesis_model)

    @test rd ≈ ustrip(u"μmol*m^-2*s^-1", v.rd)
    @test fsoil ≈ v.fsoil
    @test ci ≈ ustrip(u"μmol*mol^-1", v.ci)
    @test gs ≈ ustrip(u"mol*m^-2*s^-1", v.gs)
    @test aleaf ≈ ustrip(u"μmol*m^-2*s^-1", v.aleaf)
end

@testset "Ball Berry water potential" begin
    p = FvCBEnergyBalance(
         photosynthesis_model=FvCBPhotosynthesis(
             stomatal_conductance_model=BallBerryStomatalConductance(
                 gs_submodel=BallBerryStomatalConductanceSubModel(),
                 soil_model=PotentialSoilMethod(),
             ),
             flux_model=DukeFlux(),
             compensation_model=BadgerCollatzCompensation(),
             respiration_model=Respiration(),
         )
    )
    v = BallBerryVars()
    # Prime variables
    enbal!(v, p)

    # Run fortran
    rd, fsoil, aleaf, gs, ci, tleaf = fortran_photosyn(p, v)
    # Run julia
    photosynthesis!(v, p.photosynthesis_model)

    @test tleaf ≈ ustrip(°C, v.tleaf)
    @test rd ≈ ustrip(u"μmol*m^-2*s^-1", v.rd)
    @test fsoil ≈ v.fsoil
    @test ci ≈ ustrip(u"μmol*mol^-1", v.ci)
    @test gs ≈ ustrip(u"mol*m^-2*s^-1", v.gs)
    @test aleaf ≈ ustrip(u"μmol*m^-2*s^-1", v.aleaf)
end

@testset "Leuning" begin
    p = FvCBEnergyBalance(
         photosynthesis_model=FvCBPhotosynthesis(
             stomatal_conductance_model=BallBerryStomatalConductance(
                 gs_submodel=LeuningStomatalConductanceSubModel(),
                 soil_model=NoSoilMethod(),
             ),
             flux_model=DukeFlux(),
             compensation_model=BadgerCollatzCompensation(),
             respiration_model=Respiration(),
         )
    )
    v = BallBerryVars()
    # Prime variables
    enbal!(v, p)

    # Run fortran
    rd, fsoil, aleaf, gs, ci, tleaf = fortran_photosyn(p, v)
    # Run julia
    photosynthesis!(v, p.photosynthesis_model)

    @test tleaf ≈ ustrip(°C, v.tleaf)
    @test rd ≈ ustrip(u"μmol*m^-2*s^-1", v.rd)
    @test fsoil ≈ v.fsoil
    @test ci ≈ ustrip(u"μmol*mol^-1", v.ci)
    @test gs ≈ ustrip(u"mol*m^-2*s^-1", v.gs)
    @test aleaf ≈ ustrip(u"μmol*m^-2*s^-1", v.aleaf)
end

@testset "Medlyn" begin
    p = FvCBEnergyBalance(
        photosynthesis_model=FvCBPhotosynthesis(
            stomatal_conductance_model=BallBerryStomatalConductance(
                gs_submodel=MedlynStomatalConductanceSubModel(),
                soil_model=NoSoilMethod()),
            flux_model=DukeFlux(),
            compensation_model=BernacchiCompensation(),
            respiration_model=Respiration(),
        )
    )
    v = BallBerryVars()
    # Prime variables
    enbal!(v, p)

    # Run fortran
    rd, fsoil, aleaf, gs, ci, tleaf = fortran_photosyn(p, v)
    # Run julia
    photosynthesis!(v, p.photosynthesis_model)

    @test tleaf ≈ ustrip(°C, v.tleaf)
    @test rd ≈ ustrip(u"μmol*m^-2*s^-1", v.rd)
    @test fsoil ≈ v.fsoil
    @test ci ≈ ustrip(u"μmol*mol^-1", v.ci)
    @test gs ≈ ustrip(u"mol*m^-2*s^-1", v.gs)
    @test aleaf ≈ ustrip(u"μmol*m^-2*s^-1", v.aleaf)
end

@testset "Emax" begin
    p1 = FvCBEnergyBalance(
        photosynthesis_model=FvCBPhotosynthesis(
            stomatal_conductance_model=EmaxStomatalConductance(
                gs_submodel=BallBerryStomatalConductanceSubModel(),
                soil_model=EmaxSoilMethod()
            ),
            flux_model=DukeFlux(),
            respiration_model=Respiration(),
            compensation_model=BadgerCollatzCompensation(),
        )
    )
    p = EmaxEnergyBalance(energy_balance_model=p1)
    v = EmaxVars()
    # Prime variables
    enbal!(v, p)

    # Run fortran
    rd, fsoil, aleaf, gs, ci, tleaf = fortran_photosyn(p1, v)
    # Run julia
    photosynthesis!(v, p.energy_balance_model.photosynthesis_model)

    @test_broken fsoil ≈ v.fsoil rtol=1e-7
    @test rd ≈ v.rd.val
    @test_broken aleaf ≈ v.aleaf.val# rtol=1e-6
    @test_broken gs ≈ v.gs.val rtol=1e-7
end

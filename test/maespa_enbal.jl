using Photosynthesis, Flatten, FieldMetadata
using Photosynthesis: flux_model

using FieldMetadata: units

include(joinpath(dirname(pathof(Photosynthesis)), "../test/shared.jl"))

function run_fortran_enbal(p, v, vfun=1)

    emaxenbal = if p isa EmaxEnergyBalance 
        p
    else
        EmaxEnergyBalance()
    end

    tuzetvars = TuzetVars()
    emaxvars = EmaxVars()
    psil = Float32[ustrip(units(EmaxVars, :psil), emaxvars.psil)]
    # Use current vars or dummy if they arent needed
    if v isa TuzetVars
        tuzetvars = v
        psil = Float32[ustrip(units(TuzetVars, :psil), tuzetvars.psil)]
    elseif v isa EmaxVars
        emaxvars = v
    end

    ph = p.photosynthesis_model
    sc = ph.stomatal_conductance_model
    gs = sc.gs_submodel
    vcj = flux_model(ph.flux_model)
    vc = vcj.vcmaxformulation
    j = vcj.jmaxformulation
    sm = sc.soil_model
    # This is not even used in Maespa, its a Maestra relic but the parameters
    # are still passed in...
    jsc = JarvisStomatalConductance()

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
        vpdmin = ustrip(units(MedlynStomatalConductanceSubModel, :vpdmin), MedlynStomatalConductanceSubModel().vpdmin)
        gk = MedlynStomatalConductanceSubModel().gk
        modelgs = MEDLYN_GS
    end
    if sm isa NoSoilMethod
        wc1 = 0.0
        wc2 = 0.0
        swpexp = 0.0
        soilmoisture = v.soilmoist  
        wsoilmethod = SOILMETHOD_NONE
        soildata = SOILDATA_NONE
    elseif sm isa PotentialSoilMethod
        wc1 = 0.0
        wc2 = 0.0
        soilmoisture = ustrip(units(typeof(v), :swp), v.swp)
        swpexp = ustrip(units(typeof(sm), :swpexp), sm.swpexp)
        wsoilmethod = SOILMETHOD_POTENTIAL
        soildata = SOILDATA_POTENTIAL
    elseif sm isa VolumetricSoilMethod
        wc1 = sm.wc1
        wc2 = sm.wc2
        soilmoisture = v.soilmoist  
        swpexp = 0.0
        wsoilmethod = SOILMETHOD_VOLUMETRIC
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
    D0 = if typeof(sc.gs_submodel) <: LeuningStomatalConductanceSubModel 
        ustrip(units(LeuningStomatalConductanceSubModel, :D0), sc.gs_submodel.D0) 
    else
        0.0
    end
    par =      ustrip(units(typeof(v), :par), v.par)
    tleaf =    Float32[ustrip(°C, v.tleaf)]
    cs =       ustrip(units(typeof(v), :cs), v.cs)
    ca =       ustrip(units(typeof(v), :ca), v.ca)
    rh =       v.rh
    rnet =     ustrip(units(typeof(v), :rnet), v.rnet)
    tair =     ustrip(°C, v.tair)
    wind =     ustrip(units(typeof(v), :windspeed), v.windspeed)
    vpd =      ustrip(units(typeof(v), :vpd), v.vpd)
    press =    ustrip(units(typeof(v), :pressure), v.pressure)
    vmfd =     JarvisStomatalConductance().vmfd.val
    jmax25 =   ustrip(units(typeof(vcj.jmaxformulation), :jmax25), vcj.jmaxformulation.jmax25)
    eavj =     ustrip(units(typeof(j), :eavj), j.eavj)
    edvj =     ustrip(units(typeof(j), :edvj), j.edvj)
    delsj =    ustrip(units(typeof(j), :delsj), j.delsj)
    vcmax25 =  ustrip(units(typeof(vc), :vcmax25), vc.vcmax25)
    eavc =     ustrip(units(typeof(vc), :eavc), vc.eavc)
    if typeof(vc) <: OptimumVcmax
        edvc = ustrip(units(typeof(vc), :edvc), vc.edvc)
        delsc = ustrip(units(typeof(vc), :delsc), vc.delsc)
    else
        edvc = 0.0
        delsc = 0.0
    end
    if ph.flux_model isa DukeFlux
        tvjup = ustrip(°C, ph.flux_model.tvjup)
        tvjdn = ustrip(°C, ph.flux_model.tvjdn)
    else
        tvjup = -100.0
        tvjdn = -100.0
    end
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
    gsref =    jsc.gsref.val
    gsmin =    jsc.g0.val
    i0 =       jsc.lightmethod.i0.val
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
    plantk =   emaxenbal.plantk.val
    totsoilres = emaxenbal.totsoilres.val

    smd1 =     DeficitSoilData().smd1
    smd2 =     DeficitSoilData().smd2
    fsoil =    Float32[v.fsoil]
    g0 =       ustrip(units(typeof(sc), :g0), sc.g0)
    gamma =    sc.gs_submodel.gamma
    g1 =       sc.gs_submodel.g1
    gs =       Float32[ustrip(units(typeof(v), :gs), v.gs)]
    aleaf =    Float32[ustrip(units(typeof(v), :aleaf), v.aleaf)]
    rd =       Float32[ustrip(units(typeof(v), :rd), v.rd)]
    minleafwp = ustrip(units(EmaxVars, :minleafwp), emaxvars.minleafwp)
    ktot =     ustrip(units(EmaxVars, :ktot), emaxvars.ktot)
    weightedswp = ustrip(units(typeof(emaxvars), :swp), emaxvars.swp)
    vpara =    ustrip(units(LinearPotentialDependence, :vparb), LinearPotentialDependence().vpara)
    vparb =    ustrip(units(LinearPotentialDependence, :vparb), LinearPotentialDependence().vparb)
    vparc =    0.0 # unused
    fheat =    0.0 # unused
    etest =    0.0 # unused
    gbh =      ustrip(units(typeof(v), :gbh), v.gbh)
    sf =       ustrip(units(TuzetVars, :sf), tuzetvars.sf)
    psiv =     ustrip(units(TuzetVars, :psiv), tuzetvars.psiv)
    hmshape =  HyperbolicMinimum().hmshape
    psilin =   ustrip(units(TuzetVars, :psilin), tuzetvars.psilin)
    ci =       Float32[ustrip(units(typeof(v), :ci), v.ci)]

    iday = 1
    ihour = 1
    nsides = 1
    itermax = p.max_iterations
    yp = WangRadiationConductance()
    rdfipt =   yp.rdfipt
    tuipt =    yp.tuipt
    tdipt =    yp.tdipt
    wleaf =    ustrip(u"m", p.boundary_conductance_model.leafwidth)
    gsc =      Float32[0.0]
    et =       Float32[0.0]

    pstranspif = Libdl.dlsym(maespa_photosynlib, :pstranspif_)
    ccall(pstranspif, Nothing, (
    Ref{Int32},
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
    Ref{Int32},
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
    Ref{Bool}),
    iday,
    ihour,
    rdfipt,
    tuipt,
    tdipt,
    rnet,
    wind,
    par,
    tair,
    tmove,
    ca,
    rh,
    vpd,
    vmfd,
    press,
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
    wsoilmethod,
    emaxleaf,
    soilmoisture,
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
    wleaf,
    nsides,
    vpara,
    vparb,
    vparc,
    vfun,
    sf,
    psiv,
    itermax,
    gsc,
    aleaf,
    rd,
    et,
    fheat,
    tleaf,
    gbh,
    plantk,
    totsoilres,
    minleafwp,
    weightedswp,
    ktot,
    hmshape,
    psil,
    etest,
    ci,
    ismaespa)
    (tleaf[1], rd[1], emaxleaf[1], psil[1], fsoil[1], aleaf[1], gs[1], ci[1], et[1], gsc[1])
end

basetypeof(x) = typeof(x).name.wrapper

@testset "Tuzet inner" begin
    p = MaespaEnergyBalance(
        photosynthesis=FvCBPhotosynthesis(
            stomatal_conductance=EmaxStomatalConductance(
                gs_submodel=BallBerryStomatalConductanceSubModel(),
                soilmethod=TuzetSoilMethod(),
            ),
            flux=DukeFlux(),
            compensation=BadgerCollatzCompensation(),
            respiration_model=Respiration(),
        ),
        atol=0.005K,
    )
    for varmult in (0.95, 1.0, 1.04)
        # p = TuzetEnergyBalance(energy_balance_model=p1)
        v = TuzetVars()
        v = Flatten.modify(x -> varmult * x, v)
        println("var_multiplier: ", varmult)
        enbal!(v, p)
        tleaf, rd, emaxleaf, psil, fsoil, aleaf, gs, ci, et = run_fortran_enbal(p, v)
        @test tleaf ≈ ustrip(v.tleaf |> °C)
        @test rd ≈ v.rd.val
        @test_broken psil ≈ ustrip(u"MPa", v.psil)
        @test fsoil ≈ v.fsoil
        @test aleaf ≈ ustrip(u"μmol*m^-2*s^-1", v.aleaf)
        @test gs ≈ ustrip(u"mol*m^-2*s^-1", v.gs)
        @test ci ≈ ustrip(u"μmol*mol^-1", v.ci)
        @test et ≈ ustrip(u"μmol*m^-2*s^-1", v.et)
    end
end

@testset "Emax" begin
    p1 = MaespaEnergyBalance(
        photosynthesis=FvCBPhotosynthesis(
            stomatal_conductance=EmaxStomatalConductance(
                gs_submodel=BallBerryStomatalConductanceSubModel(),
                soilmethod=EmaxSoilMethod(),
            ),
            flux=DukeFlux(),
            compensation=BadgerCollatzCompensation(),
            respiration_model=Respiration(),
        ),
        atol=0.005K,
    )
    p = EmaxEnergyBalance(energy_balance_model=p1)
    for varmult in (0.95, 1.0, 1.001),
        parmult in (0.9999999, 1.0, 1.0000001)
        v = EmaxVars()
        v = Flatten.modify(x -> varmult * x, v)
        p = Flatten.modify(x -> x * parmult, p)
        println("var_multiplier: ", varmult)
        enbal!(v, p)
        tleaf, rd, emaxleaf, psil, fsoil, aleaf, gs, ci, et = run_fortran_enbal(p1, v);
        @test tleaf ≈ ustrip(v.tleaf |> °C)
        @test rd ≈ v.rd.val
        @test emaxleaf ≈ ustrip(u"mmol*m^-2*s^-1", v.emaxleaf)
        @test psil ≈ ustrip(u"MPa", v.psil)
        @test fsoil ≈ v.fsoil
        @test aleaf ≈ ustrip(u"μmol*m^-2*s^-1", v.aleaf)
        @test gs ≈ ustrip(u"mol*m^-2*s^-1", v.gs)
        @test ci ≈ ustrip(u"μmol*mol^-1", v.ci)
        @test et ≈ ustrip(u"μmol*m^-2*s^-1", v.et)
    end
end

function test_components(submodel, compensation, soilmethod, resp, flux, v)
    println("Testing: ", basetypeof.((submodel, compensation, soilmethod, resp, flux)))
    p = MaespaEnergyBalance(
        photosynthesis_model=FvCBPhotosynthesis(
            stomatal_conductance_model=BallBerryStomatalConductance(
                gs_submodel=submodel,
                soil_model=soilmethod,
            ),
            compensation_model=compensation,
            respiration_model=resp,
            flux_model=flux,
        ),
        atol=0.005K,
    )
    enbal!(v, p)
    tleaf, rd, emaxleaf, psil, fsoil, aleaf, gs, ci, et = run_fortran_enbal(p, v)
    @test tleaf ≈ ustrip(°C, v.tleaf)
    @test rd ≈ ustrip(u"μmol*m^-2*s^-1", v.rd)
    @test fsoil ≈ v.fsoil
    @test aleaf ≈ ustrip(u"μmol*m^-2*s^-1", v.aleaf)
    @test gs ≈ ustrip(u"mol*m^-2*s^-1", v.gs)
    @test ci ≈ ustrip(u"μmol*mol^-1", v.ci)
    @test et ≈ ustrip(u"μmol*m^-2*s^-1", v.et)
end

@testset "Test combinatorics of MaespaEnergyBalance/BallBerryStomatalConductance models" begin
    gs_submodels = (BallBerryStomatalConductanceSubModel(),
                    LeuningStomatalConductanceSubModel(),
                    MedlynStomatalConductanceSubModel(),)
    compensation = BernacchiCompensation(), BadgerCollatzCompensation()
    soilmethods = (NoSoilMethod(), PotentialSoilMethod())
    respiration = (Respiration(),)# AcclimatizedRespiration())
    flux = (DukeFlux(), Flux())
    # Multiply the variable defaults to make sure it handles more values
    # Muliplying v.tair by 1.05 (setting it to 313k) breaks things.
    # So this range is the current limit.
    var_multipliers = (0.96, 1.0, 1.03)
    par_multipliers = (0.9, 1.0, 1.04)
    EB = MaespaEnergyBalance
    for gs in gs_submodels, comp in compensation, sm in soilmethods, 
        resp in respiration, fl in flux, varmult in var_multipliers,
        parmult in par_multipliers
        v = Flatten.modify(x -> x * varmult, BallBerryVars())
        # resp fails lower
        # gs fails both
        comp, sm, fl = map((comp, sm, fl)) do p
            Flatten.modify(x -> x * parmult, p)
        end
        println("parmult: ", parmult, "varmult: ", varmult)
        println(map(x -> nameof(typeof(x)), (gs, comp, sm, resp, fl)))
        test_components(gs, comp, sm, resp, fl, v)
    end
end

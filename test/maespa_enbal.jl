using Unitful: °C, K
using Photosynthesis: fluxparams
include("shared.jl")

function run_fortran_enbal(p, v, ieco, ismaespa, modelgs, wsoilmethod, soildata, vfun)
    tzvars = TuzetVars()
    emaxvars = EmaxVars()
    ph = p.photosynthesis_model
    sc = ph.stomatal_conductance_model
    vcj = fluxparams(ph.flux_model)
    vc = vcj.vcmaxformulation
    j = vcj.jmaxformulation
    gk = 0.0
    v.rhleaf = v.rh
    d0l = 0.0
    typeof(sc.gs_submodel) <: LeuningStomatalConductanceSubModel && (d0l = sc.gs_submodel.d0l.val)
    iday = 1
    ihour = 1
    nsides = 1
    itermax = 100
    yp = YingPingRadiationConductance()
    rdfipt =   yp.rdfipt
    tuipt =    yp.tuipt
    tdipt =    yp.tdipt
    par =      v.par.val
    tleaf =    Float32[ustrip(v.tleaf |> °C)]
    tmove =    AcclimatizedRespiration().tmove.val
    cs =       v.cs.val
    ca =       v.ca.val
    rh =       v.rh
    rnet =     v.rnet.val
    tair =     ustrip(v.tair |> °C)
    wind =     0.0
    wleaf =    p.boundary_conductance_model.leafwidth.val
    gsc =      p.boundary_conductance_model.gsc.val
    vpd =      v.vpd.val
    press =    v.pressure.val
    vmfd =     JarvisStomatalConductance().vmfd.val
    jmax25 =   vcj.jmaxformulation.jmax25.val
    eavj =     j.eavj.val
    edvj =     j.edvj.val
    delsj =    j.delsj.val
    vcmax25 =  vc.vcmax25.val
    eavc =     vc.eavc.val
    edvc = 0.0
    delsc = 0.0
    # edvc =     vc.edvc.val
    # delsc =    vc.delsc.val
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
    soilmoisture = v.soilmoist
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
    d0l =      LeuningStomatalConductanceSubModel().d0l.val
    gamma =    sc.gs_submodel.gamma.val
    vpdmin =   MedlynStomatalConductanceSubModel().vpdmin.val
    g1 =       sc.gs_submodel.g1
    gk =       0.0 #ThreeParStomatalConductance().gk
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
    sf =       TuzetVars().sf.val
    psiv =     TuzetVars().psiv.val
    hmshape =  HyperbolicMinimum().hmshape
    psilin =   TuzetVars().psilin.val
    psil =     Float32[emaxvars.psil.val]
    ci =       Float32[v.ci.val]
    pstranspif = Libdl.dlsym(photosynlib, :pstranspif_)
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
    d0l,
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
    rnet,
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
    (tleaf[1], rd[1], emaxleaf[1], psil[1], fsoil[1], aleaf[1], gs[1], ci[1])
end


# Emax
sc = EmaxStomatalConductance(gs_submodel=BallBerryStomatalConductanceSubModel())
photo = FvCBPhotosynthesis(stomatal_conductance=sc,
                           flux=DukeFlux(),
                           compensation=BadgerCollatzCompensation())
p = FvCBEnergyBalance(photosynthesis=photo)
v = EmaxVars()
ieco = 1
modelgs = 2
wsoilmethod = 1
soildata = 1
ismaespa = true
vfun = 1

(p, v, ieco, ismaespa, modelgs, wsoilmethod, soildata, vfun)

(tleaf, rd, emaxleaf, psil, fsoil, aleaf, gs, ci) = 
    run_fortran_enbal(p, v, ieco, ismaespa, modelgs, wsoilmethod, soildata, vfun)
enbal!(v, p)

println(v.vcmax)

enbal!(v, p);

# println("tleaf: ", tleaf)
@test tleaf ≈ ustrip(v.tleaf |> °C)
# println("rd: ", rd)
@test rd ≈ v.rd.val
# println("emaxleaf: ", emaxleaf)
@test emaxleaf ≈ v.emaxleaf.val
# println("psil: ", psil)
@test psil ≈ v.psil.val
# println("fsoil: ", fsoil)
@test fsoil ≈ v.fsoil rtol=1e-7
# println("aleaf: ", aleaf)
@test aleaf ≈ v.aleaf.val# rtol=1e-6
# println("gs: ", gs)
@test gs ≈ v.gs.val rtol=1e-7
@test ci ≈ v.ci.val rtol=1e-6


# Ball Berry
p = FvCBEnergyBalance(
    photosynthesis_model=FvCBPhotosynthesis(
        stomatal_conductance=BallBerryStomatalConductance(
            gs_submodel=BallBerryStomatalConductance(),
            soil_model=NoSoilMethod()
        ),
        flux_model=DukeFlux(),
        compensation_model=BadgerCollatzCompensation()
    )
)
v = BallBerryVars()
ieco = 0
modelgs = 2
wsoilmethod = 1
soildata = 1
ismaespa = false

(tleaf, rd, emaxleaf, psil, fsoil, aleaf, gs, ci) = 
    run_fortran_enbal(p, v, ieco, ismaespa, modelgs, wsoilmethod, soildata, vfun)
enbal!(v, p)

@test tleaf ≈ ustrip(v.tleaf |> °C)
@test rd ≈ v.rd.val
@test fsoil ≈ v.fsoil rtol=1e-7
@test_broken aleaf ≈ v.aleaf.val# rtol=1e-6
@test_broken gs ≈ v.gs.val rtol=1e-7
@test ci ≈ v.ci.val rtol=1e-6


# Leuning

p = FvCBEnergyBalance(
    photosynthesis_model=FvCBPhotosynthesis(
        stomatal_conductance=BallBerryStomatalConductance(
            gs_submodel=LeuningStomatalConductanceSubModel(),
            soil_model=NoSoilMethod()
        ),
        flux_model=DukeFlux(),
        compensation_model=BadgerCollatzCompensation()
    )
)
v = BallBerryVars()
modelgs = 3

(tleaf, rd, emaxleaf, psil, fsoil, aleaf, gs, ci) = 
    run_fortran_enbal(p, v, ieco, ismaespa, modelgs, wsoilmethod, soildata, vfun)
enbal!(v, p)

@test tleaf ≈ ustrip(v.tleaf |> °C)
# println("tleaf: ", tleaf)
@test rd ≈ v.rd.val
# println("fsoil: ", fsoil)
@test fsoil ≈ v.fsoil rtol=1e-7
# println("aleaf: ", aleaf)
@test_broken aleaf ≈ v.aleaf.val# rtol=1e-6
# println("gs: ", gs)
@test_broken gs ≈ v.gs.val rtol=1e-7
# println("ci: ", ci)
@test ci ≈ v.ci.val rtol=1e-6

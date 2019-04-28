using Photosynthesis, Unitful, Test, Libdl
using Unitful: °C, K

photosynlib = dlopen(joinpath(ENV["MAESPA"], "physiol"))

function run_fortran_enbal(p, v, ieco, ismaespa, modelgs, wsoilmethod, soildata, vfun)
    tzvars = TuzetVars()
    emaxvars = EmaxVars()
    ph = p.photo
    vcj = ph.vcjmax
    vc = vcj.vcmaxformulation
    j = vcj.jmaxformulation
    gk = 0.0
    v.rhleaf = v.rh
    d0l = 0.0
    typeof(ph.gsmodel) <: LeuningStomatalConductance && (d0l = ph.gsmodel.d0l.val)
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
    wleaf =    p.boundary_conductance.leafwidth.val
    gsc =      p.boundary_conductance.gsc.val
    vpd =      v.vpd.val
    press =    v.pressure.val
    vmfd =     JarvisModel().vmfd.val
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
    tvjup =    ustrip(vcj.tvjup |> °C)
    tvjdn =    ustrip(vcj.tvjdn |> °C)
    theta =    ph.rubisco_regen.theta
    ajq =      ph.rubisco_regen.ajq
    rd0 =      ph.respiration.rd0.val
    q10f =     ph.respiration.q10f.val
    k10f =     AcclimatizedRespiration().k10f.val
    tref =     ustrip(ph.respiration.tref |> °C)
    rtemp =    ustrip(ph.respiration.tref |> °C)
    dayresp =  ph.respiration.dayresp
    tbelow =   ustrip(ph.respiration.tbelow |> °C)
    gsref =    JarvisModel().gsref.val
    gsmin =    JarvisModel().gsmin.val
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
    plantk =   EmaxModel().plantk.val
    totsoilres = EmaxModel().totsoilres.val
    smd1 =     DeficitSoilData().smd1
    smd2 =     DeficitSoilData().smd2
    wc1 =      VolumetricSoilMethod().wc1
    wc2 =      VolumetricSoilMethod().wc2
    swpexp =   PotentialSoilData().swpexp
    fsoil =    Float32[v.fsoil]
    g0 =       ph.g0.val
    d0l =      LeuningStomatalConductance().d0l.val
    gamma =    ph.gsmodel.gamma.val
    vpdmin =   MedlynStomatalConductance().vpdmin.val
    g1 =       ph.gsmodel.g1
    gk =       ThreeParStomatalConductance().gk
    gs =       Float32[v.gs.val]
    aleaf =    Float32[v.aleaf.val]
    rd =       Float32[v.rd.val]
    minleafwp = emaxvars.minleafwp.val
    ktot =     emaxvars.ktot.val
    weightedswp = emaxvars.weightedswp.val
    vpara =    LinearPotentialDependence().vpara.val
    vparb =    LinearPotentialDependence().vparb.val
    vparc =    0.0 # unused
    fheat =    0.0 # unused
    etest =    0.0 # unused
    gbh =      v.gbh.val
    sf =       TuzetVars().sf.val
    psiv =     TuzetVars().psiv.val
    hmshape =  HyperbolicMinimumGS().hmshape
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

p = FvCBEnergyBalance(photo=EmaxModel(gsmodel=BallBerryStomatalConductance(),
                      vcjmax=DukeVcJmax(),
                      compensation=BadgerCollatzCompensation()))
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
run_enbal!(p, v)

println(v.vcmax)

enbal!(p, v);

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

p = FvCBEnergyBalance(photo=BallBerryModel(gsmodel=BallBerryStomatalConductance(),
                      vcjmax=DukeVcJmax(),
                      compensation=BadgerCollatzCompensation()))
v = BallBerryVars()
ieco = 0
modelgs = 2
wsoilmethod = 1
soildata = 1
ismaespa = false

(tleaf, rd, emaxleaf, psil, fsoil, aleaf, gs, ci) = 
    run_fortran_enbal(p, v, ieco, ismaespa, modelgs, wsoilmethod, soildata, vfun)
run_enbal!(p, v)

@test tleaf ≈ ustrip(v.tleaf |> °C)
@test rd ≈ v.rd.val
@test fsoil ≈ v.fsoil rtol=1e-7
@test_broken aleaf ≈ v.aleaf.val# rtol=1e-6
@test_broken gs ≈ v.gs.val rtol=1e-7
@test ci ≈ v.ci.val rtol=1e-6


# Leuning

p = FvCBEnergyBalance(photo=BallBerryModel(gsmodel=LeuningStomatalConductance(),
                      vcjmax=DukeVcJmax(),
                      compensation=BadgerCollatzCompensation()))
v = BallBerryVars()
modelgs = 3

(tleaf, rd, emaxleaf, psil, fsoil, aleaf, gs, ci) = 
    run_fortran_enbal(p, v, ieco, ismaespa, modelgs, wsoilmethod, soildata, vfun)
run_enbal!(p, v)

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

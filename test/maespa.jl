using Libdl
using Photosynthesis: quad

# test against the original FORTRAN routines
global emaxplant = EnergyBalance(photo=FvCBPhoto(model=EmaxModel()))
global p = emaxplant.photo
global v = EmaxVars()
global v.tleaf = 15u"°C"
global photosynlib = dlopen(joinpath(ENV["MAESPA"], "physiol.so"))

@testset "funcs" begin

    global resp = Libdl.dlsym(photosynlib, :resp_)
    global f = p.respiration
    global rd0 = ustrip(f.rd0)
    global rdacc = 1.0
    global tleaf = v.tleaf.val
    global q10f = f.q10f.val
    global tref = f.tref.val
    global dayresp = f.dayresp
    global tbelow = f.tbelow.val
    global resp_ref = ccall(resp, Float32, (Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}
                             ), rd0, rdacc, tleaf, q10f, tref, dayresp, tbelow)
    v.rd = respiration(f, v, p)
    # @test v.rd.val == resp_ref

    global vcmaxtfn = Libdl.dlsym(photosynlib, :vcmaxtfn_)
    global v = EmaxVars()
    global p = emaxplant.photo
    v.tleaf = 15.0u"°C"
    global tleaf = v.tleaf.val
    global f = VcJmax(vcmaxformulation=NoOptimumVcmax())
    global vc = f.vcmaxformulation
    global vcmax25 = vc.vcmax25.val
    global eavc = vc.eavc.val
    global edvc = 0.0
    global delsc = 0.0
    global tvjup = -100.0
    global tvjdn = -100.0
    global vcmax_ref = ccall(vcmaxtfn, Float32, (Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}
                     ), vcmax25,tleaf,eavc,edvc,delsc,tvjup,tvjdn)
    @test ustrip( max_rubisco_activity(f, v, p)) ≈ vcmax_ref rtol=1e-4


    global v = EmaxVars()
    global p = emaxplant.photo
    v.tleaf = 15.0u"°C"
    global tleaf = v.tleaf.val
    global f = VcJmax(vcmaxformulation=OptimumVcmax())
    global vc = f.vcmaxformulation
    global eavc = vc.eavc.val
    global edvc = vc.edvc.val
    global delsc = vc.delsc.val
    global tvjup = -100.0
    global tvjdn = -100.0
    global vcmax_ref = ccall(vcmaxtfn, Float32, (Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}
                             ), vcmax25,tleaf,eavc,edvc,delsc,tvjup,tvjdn)
    @test ustrip(max_rubisco_activity(f, v, p)) ≈ vcmax_ref rtol=1e-4

    global v = EmaxVars()
    global p = emaxplant.photo
    v.tleaf = 15.0u"°C"
    global tleaf = v.tleaf.val
    global f = DukeVcJmax(vcmaxformulation=OptimumVcmax())
    global vc = f.vcmaxformulation
    global vcmax25 = vc.vcmax25.val
    global eavc = vc.eavc.val
    global edvc = vc.edvc.val
    global delsc = vc.delsc.val
    global tvjup = f.tvjup.val
    global tvjdn = f.tvjdn.val
    global vcmax_ref = ccall(vcmaxtfn, Float32, (Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}
                             ), vcmax25,tleaf,eavc,edvc,delsc,tvjup,tvjdn)
    @test ustrip(max_rubisco_activity(f, v, p)) ≈ vcmax_ref rtol=1e-4


    global jmaxtfn = Libdl.dlsym(photosynlib, :jmaxtfn_)
    global f = Jmax()
    v.tleaf = 15.0u"°C"
    global tleaf = v.tleaf.val
    global jmax25 = f.jmax25.val
    global eavj = f.eavj.val
    global edvj = f.edvj.val
    global delsj = f.delsj.val
    global tvjup = -100.0
    global tvjdn = -100.0
    global jmax_ref = ccall(jmaxtfn, Float32, (Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}
                             ), jmax25,tleaf,eavj,edvj,delsj,tvjup,tvjdn)
    @test ustrip( max_electron_transport_rate(f, v, p)) ≈ jmax_ref rtol=1e-4

end

global v = EmaxVars()
global p = EnergyBalance(photo=FvCBPhoto(model=EmaxModel(gsmodel=BallBerryStomatalConductance()),
                                  vcjmax=DukeVcJmax(),
                                  compensation = BadgerCollatzCompensation(),
                                 ))
global ph = p.photo
global mod = p.photo.model
global vcj = ph.vcjmax
global vc = vcj.vcmaxformulation
global j = vcj.jmaxformulation
global mgs = 2
global wsm = 1
global iec = 1
global soild = 1
global gk = 0.0
global ismaespa = true
v.rhleaf = v.rh
global tzv = TuzetVars()
global d0l = 0.0
typeof(ph.model.gsmodel) <: LeuningStomatalConductance && (d0l = ph.model.gsmodel.d0l.val)
phototranspiration!(v, p)
photosynthesis!(v, ph)

global par =      v.par.val
global tleaf =    v.tleaf.val
global tmove =    AcclimatizedRespiration().tmove.val
global cs =       v.cs.val
global rh =       v.rh
global vpd =      v.vpd.val
global vmfd =     JarvisModel().vmfd.val
global jmax25 =   vcj.jmaxformulation.jmax25.val
global ieco =     iec 
global eavj =     j.eavj.val
global edvj =     j.edvj.val
global delsj =    j.delsj.val
global vcmax25 =  vc.vcmax25.val
global eavc =     vc.eavc.val
global edvc = 0.0
global delsc = 0.0
# edvc =     vc.edvc.val
# delsc =    vc.delsc.val
global tvjup =    vcj.tvjup.val
global tvjdn =    vcj.tvjdn.val
global theta =    ph.rubisco_regen.theta
global ajq =      ph.rubisco_regen.ajq
global rd0 =      ph.respiration.rd0.val
global q10f =     ph.respiration.q10f.val
global k10f =     AcclimatizedRespiration().k10f.val
global tref =    ph.respiration.tref.val
global dayresp =  ph.respiration.dayresp
global tbelow =   ph.respiration.tbelow.val
global modelgs =  mgs
global gsref =    JarvisModel().gsref.val
global gsmin =    JarvisModel().gsmin.val
global i0 =       JarvisLight().i0.val
global d0 =       JarvisLinearDeclineVPD().d0.val
global vk1 =      JarvisHyperbolicVPD().vk1
global vk2 =      JarvisHyperbolicVPD().vk2
global vpd1 =     JarvisLohammerVPD().vpd1.val
global vpd2 =     JarvisLohammerVPD().vpd2.val
global vmfd0 =    JarvisFractionDeficitVPD().vmfd0.val
global gsja =     JarvisLinearCO2().gsja.val
global gsjb =     JarvisNonlinearCO2().gsjb.val
global t0 =       JarvisTemp1().t0.val
global tref =     JarvisTemp1().tref.val
global tmax =     JarvisTemp1().tmax.val
global wsoilmethod = wsm
global soilmoisture = v.soilmoist
global emaxleaf = v.emaxleaf.val
global smd1 =     DeficitSoilData().smd1
global smd2 =     DeficitSoilData().smd2
global wc1 =      VolumetricSoilMethod().wc1
global wc2 =      VolumetricSoilMethod().wc2
global soildata = soild
global swpexp =   PotentialSoilData().swpexp
global fsoil =    Float32[v.fsoil]
global g0 =       mod.g0.val
global d0l =      LeuningStomatalConductance().d0l.val
global γ =        mod.gsmodel.gamma.val
global vpdmin =   MedlynStomatalConductance().vpdmin.val
global g1 =       mod.gsmodel.g1
global gk =       ThreeParStomatalConductance().gk
global gs =       Float32[v.gs.val]
global aleaf =    Float32[v.aleaf.val]
global rd =       v.rd.val
global minleafwp = v.minleafwp.val
global ktot =     v.ktot.val
global weightedswp = v.weightedswp.val
global vpara =    LinearPotentialDependence().vpara.val
global vparb =    LinearPotentialDependence().vparb.val
global vparc =    0.0 # unused
global vfun =     1
global sf =       TuzetVars().sf.val
global psiv =     TuzetVars().psiv.val
global hmshape =  HyperbolicMinimumGS().hmshape
global psilin =   TuzetVars().psilin.val
global psil =     Float32[v.psil.val]
global ci =       Float32[v.ci.val]
global ismaespa = true

@testset "photosynthesis" begin

    global photosyn = Libdl.dlsym(photosynlib, :photosyn_)
    ccall(photosyn,
        Nothing,
        (
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
    Ref{Bool}),
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
    tref,
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
    γ,
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

    gs[1]
    @test tleaf[1] == v.tleaf.val
    @test rd[1] == v.rd.val
    # @test km[1] ≈ v.km.val rtol=1e-4
    # @test gammastar[1] ≈ v.gammastar.val rtol=1e-4
    # @test v.jmax.val ≈ jmax[1] rtol=1e-4
    # @test v.vcmax.val ≈ vcmax[1] rtol=1e-4
    # @test v.vj.val == vj[1]
    # @test gsdv[1] == v.gsdiva.val
    # @test aj[1] ≈ v.aj.val rtol=1e-7
    # @test ac[1] ≈ v.ac.val rtol=1e-6
    @test emaxleaf[1] == v.emaxleaf.val
    @test psil[1] ≈ v.psil.val
    @test fsoil[1] ≈ v.fsoil rtol=1e-7
    @test_broken aleaf[1] == v.aleaf.val# rtol=1e-6
    @test_broken gs[1] ≈ v.gs.val rtol=1e-7
    @test ci[1] ≈ v.ci.val rtol=1e-6

end


@testset "other funcs" begin

    global kmfn = Libdl.dlsym(photosynlib, :kmfn_)
    global ieco = 0 # Bernacci
    global tleaf = v.tleaf.val
    global km_ref = ccall(kmfn, Float32, (Ref{Float32}, Ref{Int32}), tleaf, ieco)
    global km = rubisco_compensation_point(BernacchiCompensation(), v, p) # Michaelis-Menten for Rubisco, umol mol-1
    @test km.val ≈ km_ref rtol=1e-4
    global ieco = 1 # Badger-Collatz
    global tleaf = v.tleaf.val
    global km_ref = ccall(kmfn, Float32, (Ref{Float32}, Ref{Int32}), tleaf,ieco)
    global km = rubisco_compensation_point(BadgerCollatzCompensation(), v, p) # Michaelis-Menten for Rubisco, umol mol-1
    @test km.val ≈ km_ref rtol=1e-4

    global gammafn = Libdl.dlsym(photosynlib, :gammafn_)
    global gammastar_ref = ccall(gammafn, Float32, (Ref{Float32}, Ref{Int32}), tleaf, 0)
    global gammastar = co2_compensation_point(BernacchiCompensation(), v, p) # Michaelis-Menten for Rubisco, umol mol-1
    @test gammastar.val ≈ gammastar_ref rtol=1e-4
    global gammastar_ref = ccall(gammafn, Float32, (Ref{Float32}, Ref{Int32}), tleaf, 1)
    global gammastar = co2_compensation_point(BadgerCollatzCompensation(), v, p) # Michaelis-Menten for Rubisco, umol mol-1
    @test gammastar.val ≈ gammastar_ref rtol=1e-7

    global arrhfn = Libdl.dlsym(photosynlib, :arrh_)
    global arrh_ref = ccall(arrhfn, Float32, (Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}), 42.75, 37830.0, 30.0, 25.0)
    global arrh = arrhenius(42.75u"μmol*mol^-1", 37830.0u"J*mol^-1", 30.0u"°C" |> u"K", 25.0u"°C" |> u"K")
    @test arrh.val ≈ arrh_ref
    global vcmaxtfn = Libdl.dlsym(photosynlib, :vcmaxtfn_)
    global f = DukeVcJmax(vcmaxformulation=OptimumVcmax())
    global vcmax_ref = ccall(vcmaxtfn, Float32, (Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}
                             ), vcmax25,tleaf,eavc,edvc,delsc,tvjup,tvjdn)
    @test ustrip(max_rubisco_activity(f, v, p)) ≈ vcmax_ref

    global quadm_fort = Libdl.dlsym(photosynlib, :quadm_)
    global f = p.photo.rubisco_regen
    global a = theta[1]
    global b = -(ajq[1] * par[1] + v.jmax.val)
    global c = ajq[1] * par[1] * v.jmax.val
    global quad_ref = ccall(quadm_fort, Float32, (Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Int32}), a, b, c, 1)/4.0
    global quad_test = quad(Photosynthesis.Lower(), a,b,c)/4
    @test quad_ref ≈ quad_test atol=1e-5

end

nothing

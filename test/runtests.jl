using Photosynthesis
using Unitful

@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end


# TODO improve these tests. This was inherited without any testing, 
# it will take time to develope suitable tests for all the processes.

# Make sure everything constructs with defaults
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
p = EnergyBalance()
v = EmaxVars()


# Make sure everything actually runs
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


# test against the original FORTRAN routines
p = emaxplant
v = EmaxVars()
v.tleaf = 15u"°C"
p = emaxplant.photo
photosynlib = Libdl.dlopen("test/physiol.so")

resp = Libdl.dlsym(photosynlib, :resp_)
f = p.respiration
rd0 = ustrip(f.rd0)
rdacc = 1.0
tleaf = v.tleaf.val
q10f = f.q10f.val
rtemp = f.rtemp.val
dayresp = f.dayresp
tbelow = f.tbelow.val
resp_ref = ccall(resp, Float32, (Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}
                         ), rd0, rdacc, tleaf, q10f, rtemp, dayresp, tbelow)
v.rd = respiration(f, v, p)
# @test v.rd.val == resp_ref

vcmaxtfn = Libdl.dlsym(photosynlib, :vcmaxtfn_)
v = EmaxVars()
p = emaxplant.photo
v.tleaf = 15.0u"°C"
tleaf = v.tleaf.val
f = VcJmax(vcmaxformulation=NoOptimumVcmax())
vc = f.vcmaxformulation
vcmax25 = vc.vcmax25.val
eavc = vc.eavc.val
edvc = 0.0
delsc = 0.0
tvjup = -100.0
tvjdn = -100.0
vcmax_ref = ccall(vcmaxtfn, Float32, (Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}
                 ), vcmax25,tleaf,eavc,edvc,delsc,tvjup,tvjdn)
@test ustrip( max_rubisco_activity(f, v, p)) ≈ vcmax_ref rtol=1e-4


v = EmaxVars()
p = emaxplant.photo
v.tleaf = 15.0u"°C"
tleaf = v.tleaf.val
f = VcJmax(vcmaxformulation=OptimumVcmax())
vc = f.vcmaxformulation
eavc = vc.eavc.val
edvc = vc.edvc.val
delsc = vc.delsc.val
tvjup = -100.0
tvjdn = -100.0
vcmax_ref = ccall(vcmaxtfn, Float32, (Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}
                         ), vcmax25,tleaf,eavc,edvc,delsc,tvjup,tvjdn)
@test ustrip(max_rubisco_activity(f, v, p)) ≈ vcmax_ref rtol=1e-4

v = EmaxVars()
p = emaxplant.photo
v.tleaf = 15.0u"°C"
tleaf = v.tleaf.val
f = DukeVcJmax(vcmaxformulation=OptimumVcmax())
vc = f.vcmaxformulation
vcmax25 = vc.vcmax25.val
eavc = vc.eavc.val
edvc = vc.edvc.val
delsc = vc.delsc.val
tvjup = f.tvjup.val
tvjdn = f.tvjdn.val
vcmax_ref = ccall(vcmaxtfn, Float32, (Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}
                         ), vcmax25,tleaf,eavc,edvc,delsc,tvjup,tvjdn)
@test ustrip(max_rubisco_activity(f, v, p)) ≈ vcmax_ref rtol=1e-4


jmaxtfn = Libdl.dlsym(photosynlib, :jmaxtfn_)
f = Jmax()
v.tleaf = 15.0u"°C"
tleaf = v.tleaf.val
jmax25 = f.jmax25.val
eavj = f.eavj.val
edvj = f.edvj.val
delsj = f.delsj.val
tvjup = -100.0
tvjdn = -100.0
jmax_ref = ccall(jmaxtfn, Float32, (Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}
                         ), jmax25,tleaf,eavj,edvj,delsj,tvjup,tvjdn)
@test ustrip( max_electron_transport_rate(f, v, p)) ≈ jmax_ref rtol=1e-4


photosyn = Libdl.dlsym(photosynlib, :photosyn_)
v = EmaxVars()
p = EnergyBalance(photo=FvCBPhoto(model=EmaxModel(gsmodel=BallBerryStomatalConductance()),
                                  vcjmax=DukeVcJmax(),
                                  compensation = BadgerCollatzCompensation(),
                                 ))
ph = p.photo
mod = p.photo.model
vcj = ph.vcjmax
vc = vcj.vcmaxformulation
j = vcj.jmaxformulation
mgs = 2
wsm = 1
iec = 1
soild = 1
gk = 0.0
ismaespa = true
v.rhleaf = v.rh
tzv = TuzetVars()
d0l = 0.0
typeof(ph.model.gsmodel) <: LeuningStomatalConductance && (d0l = ph.model.gsmodel.d0l.val)
phototranspiration!(v, p)
photosynthesis!(v, ph)

par =      v.par.val
tleaf =    v.tleaf.val
tmove =    AcclimatizedRespiration().tmove.val
cs =       v.cs.val
rh =       v.rh
vpd =      v.vpd.val
vmfd =     JarvisModel().vmfd.val
jmax25 =   vcj.jmaxformulation.jmax25.val
ieco =     iec 
eavj =     j.eavj.val
edvj =     j.edvj.val
delsj =    j.delsj.val
vcmax25 =  vc.vcmax25.val
eavc =     vc.eavc.val
edvc =     vc.edvc.val
delsc =    vc.delsc.val
tvjup =    vcj.tvjup.val
tvjdn =    vcj.tvjdn.val
theta =    ph.rubisco_regen.theta
ajq =      ph.rubisco_regen.ajq
rd0 =      ph.respiration.rd0.val
q10f =     ph.respiration.q10f.val
k10f =     AcclimatizedRespiration().k10f.val
rtemp =    ph.respiration.rtemp.val
dayresp =  ph.respiration.dayresp
tbelow =   ph.respiration.tbelow.val
modelgs =  mgs
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
t0 =       JarvisTemp1().t0.val
tref =     JarvisTemp1().tref.val
tmax =     JarvisTemp1().tmax.val
wsoilmethod = wsm
soilmoisture = v.soilmoist
emaxleaf = v.emaxleaf.val
smd1 =     DeficitSoilData().smd1
smd2 =     DeficitSoilData().smd2
wc1 =      VolumetricSoilMethod().wc1
wc2 =      VolumetricSoilMethod().wc2
soildata = soild
swpexp =   PotentialSoilData().swpexp
fsoil =    Float32[v.fsoil]
g0 =       mod.g0.val
d0l =      LeuningStomatalConductance().d0l.val
gamma =    mod.gsmodel.gamma.val
vpdmin =   MedlynStomatalConductance().vpdmin.val
g1 =       mod.gsmodel.g1
gk =       ThreeParStomatalConductance().gk
gs =       Float32[v.gs.val]
aleaf =    Float32[v.aleaf.val]
rd =       v.rd.val
minleafwp = v.minleafwp.val
ktot =     v.ktot.val
weightedswp = v.weightedswp.val
vpara =    LinearPotentialDependence().vpara.val
vparb =    LinearPotentialDependence().vparb.val
vparc =    0.0 # unused
vfun =     1
sf =       TuzetVars().sf.val
psiv =     TuzetVars().psiv.val
hmshape =  HyperbolicMinimumGS().hmshape
psilin =   TuzetVars().psilin.val
psil =     Float32[v.psil.val]
ci =       Float32[v.ci.val]
ismaespa = true

ccall(photosyn,
    Void,
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
@test aleaf[1] == v.aleaf.val# rtol=1e-6
@test gs[1] ≈ v.gs.val rtol=1e-7
@test ci[1] ≈ v.ci.val rtol=1e-6

kmfn = Libdl.dlsym(photosynlib, :kmfn_)
ieco = 0 # Bernacci
tleaf = v.tleaf.val
km_ref = ccall(kmfn, Float32, (Ref{Float32}, Ref{Int32}), tleaf, ieco)
km = rubisco_compensation_point(BernacchiCompensation(), v, p) # Michaelis-Menten for Rubisco, umol mol-1
@test km.val ≈ km_ref rtol=1e-4

ieco = 1 # Badger-Collatz
tleaf = v.tleaf.val
km_ref = ccall(kmfn, Float32, (Ref{Float32}, Ref{Int32}), tleaf,ieco)
km = rubisco_compensation_point(BadgerCollatzCompensation(), v, p) # Michaelis-Menten for Rubisco, umol mol-1
@test km.val ≈ km_ref rtol=1e-4

gammafn = Libdl.dlsym(photosynlib, :gammafn_)
gammastar_ref = ccall(gammafn, Float32, (Ref{Float32}, Ref{Int32}), tleaf, 0)
gammastar = co2_compensation_point(BernacchiCompensation(), v, p) # Michaelis-Menten for Rubisco, umol mol-1
@test gammastar.val ≈ gammastar_ref rtol=1e-4

gammastar_ref = ccall(gammafn, Float32, (Ref{Float32}, Ref{Int32}), tleaf, 1)
gammastar = co2_compensation_point(BadgerCollatzCompensation(), v, p) # Michaelis-Menten for Rubisco, umol mol-1
@test gammastar.val ≈ gammastar_ref rtol=1e-7

arrhfn = Libdl.dlsym(photosynlib, :arrh_)
arrh_ref = ccall(arrhfn, Float32, (Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}), 42.75, 37830.0, 30.0, 25.0)
arrh = arrhenius(42.75u"μmol*mol^-1", 37830.0u"J*mol^-1", 30.0u"°C" |> u"K", 25.0u"°C" |> u"K")
@test arrh.val ≈ arrh_ref

f = DukeVcJmax(vcmaxformulation=OptimumVcmax())
vcmax_ref = ccall(vcmaxtfn, Float32, (Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}
                         ), vcmax25,tleaf,eavc,edvc,delsc,tvjup,tvjdn)
@test ustrip(max_rubisco_activity(f, v, p)) ≈ vcmax_ref rtol=1e-4


quadm_fort = Libdl.dlsym(photosynlib, :quadm_)
f = p.photo.rubisco_regen
a = theta[1]
b = -(ajq[1] * par[1] + v.jmax.val)
c = ajq[1] * par[1] * v.jmax.val
vj1 = ccall(quadm_fort, Float32, (Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Int32}), a, b, c, 1)/4.0
SimpleRoots.quadm(a,b,c)/4
vj2 = rubisco_regeneration(p.photo.rubisco_regen, v, p)


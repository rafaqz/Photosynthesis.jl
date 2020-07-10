using Unitful: °C, K

# include("shared.jl")

# Setup
emax = FvCBEnergyBalance(
    photosynthesis_model=FvCBPhotosynthesis(
        stomatal_conductance_model=EmaxStomatalConductance()
    )
)
ph = emax.photosynthesis_model
v = EmaxVars()
v.tleaf = 15°C
photosynlib = dlopen(joinpath(ENV["MAESPA"], "physiol"))


# Rubisco compensation
# kmfn: rubisco_compensation_point
kmfn_fortran = Libdl.dlsym(photosynlib, :kmfn_)
tleaf = ustrip(v.tleaf |> °C)

ieco = 0 # Bernacci
km_ref = ccall(kmfn_fortran, Float32, (Ref{Float32}, Ref{Int32}), tleaf, ieco)
km = rubisco_compensation_point(BernacchiCompensation(), v.tleaf) # Michaelis-Menten for Rubisco, umol mol-1
@test km.val ≈ km_ref rtol=1e-4

ieco = 1 # Badger-Collatz
km_ref = ccall(kmfn_fortran, Float32, (Ref{Float32}, Ref{Int32}), tleaf,ieco)
km = rubisco_compensation_point(BadgerCollatzCompensation(), v.tleaf) # Michaelis-Menten for Rubisco, umol mol-1
@test km.val ≈ km_ref rtol=1e-4


# CO2 compensation
# gammafn: co2_compensation_point
gammafn_fortran = Libdl.dlsym(photosynlib, :gammafn_)
gammastar_ref = ccall(gammafn_fortran, Float32, (Ref{Float32}, Ref{Int32}), tleaf, BERNACCI)
gammastar = co2_compensation_point(BernacchiCompensation(), v.tleaf) # Michaelis-Menten for Rubisco, umol mol-1
@test gammastar.val ≈ gammastar_ref rtol=1e-4 

gammastar_ref = ccall(gammafn_fortran, Float32, (Ref{Float32}, Ref{Int32}), tleaf, BADGERCOLLATZ)
gammastar = co2_compensation_point(BadgerCollatzCompensation(), v.tleaf) # Michaelis-Menten for Rubisco, umol mol-1
@test gammastar.val ≈ gammastar_ref rtol=1e-4


# Vcmax
# vcmaxtfn: max_rubisco_activity
vcmaxtfn_fortran = Libdl.dlsym(photosynlib, :vcmaxtfn_)
v = EmaxVars()
v.tleaf = 15.0°C |> K
tleaf = ustrip(v.tleaf |> °C)
vc = NoOptimumVcmax()
vcmax25 = vc.vcmax25.val
eavc = vc.eavc.val
edvc = 0.0
delsc = 0.0
tvjup = -100.0
tvjdn = -100.0
vcmax_ref = ccall(vcmaxtfn_fortran, Float32, (Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}
                 ), vcmax25, tleaf, eavc, edvc, delsc, tvjup, tvjdn)
@test ustrip(max_rubisco_activity(vc, v.tleaf)) ≈ vcmax_ref rtol=1e-4

vcmax_ref = ccall(vcmaxtfn_fortran, Float32, (Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}
                         ), vcmax25, tleaf, eavc, edvc, delsc, tvjup, tvjdn)
@test ustrip(max_rubisco_activity(vc, v.tleaf)) ≈ vcmax_ref

v = EmaxVars()
v.tleaf = 15.0°C
tleaf = ustrip(v.tleaf |> °C)
vc = OptimumVcmax()
eavc = vc.eavc.val
edvc = vc.edvc.val
delsc = vc.delsc.val
tvjup = -100.0
tvjdn = -100.0
vcmax_ref = ccall(vcmaxtfn_fortran, Float32,
                  (Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}
                  ), vcmax25, tleaf, eavc, edvc, delsc, tvjup, tvjdn)
@test ustrip(max_rubisco_activity(vc, v.tleaf)) ≈ vcmax_ref rtol=1e-4


# Maximum electron transport
# jmaxtfn: max_electron_transport_rate
jmaxtfn_fortran = Libdl.dlsym(photosynlib, :jmaxtfn_)
f = Jmax()
v.tleaf = 15.0°C
tleaf = ustrip(v.tleaf |> °C)
jmax25 = f.jmax25.val
eavj = f.eavj.val
edvj = f.edvj.val
delsj = f.delsj.val
tvjup = -100.0
tvjdn = -100.0
jmax_ref = ccall(jmaxtfn_fortran, Float32, 
                 (Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}), 
                 jmax25, tleaf, eavj, edvj, delsj, tvjup, tvjdn)
@test ustrip(max_electron_transport_rate(f, v.tleaf)) ≈ jmax_ref rtol=1e-4


# Maximum rubisco activity
# vcmaxtfn: max_rubisco_activity
v = EmaxVars()
v.tleaf = 15.0°C
tleaf = ustrip(v.tleaf |> °C)
vc = OptimumVcmax()
f = DukeFlux()
vcmax25 = vc.vcmax25.val
eavc = vc.eavc.val
edvc = vc.edvc.val
delsc = vc.delsc.val
tvjup = ustrip(f.tvjup |> °C)
tvjdn = ustrip(f.tvjdn |> °C)
vcmax_ref = ccall(vcmaxtfn_fortran, Float32, 
                  (Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}), 
                  vcmax25, tleaf, eavc, edvc, delsc, tvjup, tvjdn)
@test ustrip(max_rubisco_activity(vc, v.tleaf)) ≈ vcmax_ref rtol=1e-4


# Respriation
# resp: respiration
resp_fortran = Libdl.dlsym(photosynlib, :resp_)
f = Respiration()
rd0 = ustrip(f.rd0)
rdacc = 1.0
tleaf = ustrip(v.tleaf |> °C)
q10f = f.q10f.val
tref = ustrip(f.tref |> °C)
dayresp = f.dayresp
tbelow = ustrip(f.tbelow |> °C)
k10f = 0.0 # No acclimation 
tmove = 0.0 # No acclimation 
resp_ref = ccall(resp_fortran, Float32,
                 (Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}),
                 rd0, rdacc, tleaf, tmove, q10f, k10f, tref, dayresp, tbelow)
v.rd = respiration(f, v.tleaf)
@test v.rd.val ≈ resp_ref


# Acclimatised Respriation
# resp: respiration
resp_fortran = Libdl.dlsym(photosynlib, :resp_)
f = AcclimatizedRespiration()
rd0 = ustrip(f.rd0)
rdacc = 1.0 # this isn't actually a parameter
tleaf = ustrip(v.tleaf |> °C)
q10f = f.q10f.val
tref = ustrip(f.tref |> °C)
dayresp = f.dayresp
tbelow = ustrip(f.tbelow |> °C)
k10f = f.k10f.val
tmove = f.tmove.val
resp_ref = ccall(resp_fortran, Float32,
                 (Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}),
                 rd0, rdacc, tleaf, tmove, q10f, k10f, tref, dayresp, tbelow)
v.rd = respiration(f, v.tleaf)

# This is actually commented out in the maespa FORTRAN
# @test v.rd.val ≈ resp_ref

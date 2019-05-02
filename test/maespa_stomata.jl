using Photosynthesis, Unitful, Test, Libdl
using Unitful: °C, K

# Setup
emax = FvCBEnergyBalance(photosynthesis=EmaxPhotosynthesis())
ph = emax.photosynthesis
v = EmaxVars()
v.tleaf = 15°C |> K


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


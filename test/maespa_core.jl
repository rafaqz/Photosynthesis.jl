using Photosynthesis

include(joinpath(dirname(pathof(Photosynthesis)), "../test/shared.jl"))

# Setup
emax = FvCBEnergyBalance(
    photosynthesis_model=FvCBPhotosynthesis(
        stomatal_conductance_model=EmaxStomatalConductance()
    )
)
ph = emax.photosynthesis_model
v = EmaxVars()
v.tleaf = 15°C
tleaf = ustrip(°C, v.tleaf)

@testset "rubisco_compensation_point/kmfn" begin
    kmfn_fortran = Libdl.dlsym(photosynlib, :kmfn_)
    ieco = BERNACCI
    km_ref = ccall(kmfn_fortran, Float32, (Ref{Float32}, Ref{Int32}), tleaf, ieco)
    km = rubisco_compensation_point(BernacchiCompensation(), v.tleaf) # Michaelis-Menten for Rubisco, umol mol-1
    @test ustrip(u"μmol/mol", km) ≈ km_ref

    ieco = BADGERCOLLATZ
    km_ref = ccall(kmfn_fortran, Float32, (Ref{Float32}, Ref{Int32}), tleaf,ieco)
    km = rubisco_compensation_point(BadgerCollatzCompensation(), v.tleaf) # Michaelis-Menten for Rubisco, umol mol-1
    @test ustrip(u"μmol/mol", km) ≈ km_ref
end


@testset "co2_compensation_point/GAMMAFN" begin
    gammafn_fortran = Libdl.dlsym(photosynlib, :gammafn_)
    gammastar_ref = ccall(gammafn_fortran, Float32, (Ref{Float32}, Ref{Int32}), tleaf, BERNACCI)
    gammastar = co2_compensation_point(BernacchiCompensation(), v.tleaf) # Michaelis-Menten for Rubisco, umol mol-1
    @test ustrip(u"μmol/mol", gammastar) ≈ gammastar_ref

    gammastar_ref = ccall(gammafn_fortran, Float32, (Ref{Float32}, Ref{Int32}), tleaf, BADGERCOLLATZ)
    gammastar = co2_compensation_point(BadgerCollatzCompensation(), v.tleaf) # Michaelis-Menten for Rubisco, umol mol-1
    @test ustrip(u"μmol/mol", gammastar) ≈ gammastar_ref
end


@testset "max_rubisco_activity/VCMAXTFN" begin
    vcmaxtfn_fortran = Libdl.dlsym(photosynlib, :vcmaxtfn_)
    vc = NoOptimumVcmax()
    vcmax25 = ustrip(u"μmol*m^-2*s^-1", vc.vcmax25)
    eavc = ustrip(u"J*mol^-1", vc.eavc)
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

    vc = OptimumVcmax()
    eavc = ustrip(u"J*mol^-1", vc.eavc)
    edvc = ustrip(u"J*mol^-1", vc.edvc)
    delsc = ustrip(u"J*K^-1*mol^-1", vc.delsc)
    tvjup = -100.0
    tvjdn = -100.0
    vcmax_ref = ccall(vcmaxtfn_fortran, Float32,
                      (Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}
                      ), vcmax25, tleaf, eavc, edvc, delsc, tvjup, tvjdn)
    @test ustrip(max_rubisco_activity(vc, v.tleaf)) ≈ vcmax_ref

    vc = OptimumVcmax()
    f = DukeFlux()
    vcmax25 = ustrip(u"μmol*m^-2*s^-1", vc.vcmax25)
    eavc = ustrip(u"J*mol^-1", vc.eavc)
    edvc = ustrip(u"J*mol^-1", vc.edvc)
    delsc = ustrip(u"J*K^-1*mol^-1", vc.delsc)
    tvjup = ustrip(°C, f.tvjup)
    tvjdn = ustrip(°C, f.tvjdn)
    vcmax_ref = ccall(vcmaxtfn_fortran, Float32, 
                      (Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}), 
                      vcmax25, tleaf, eavc, edvc, delsc, tvjup, tvjdn)
    @test ustrip(max_rubisco_activity(vc, v.tleaf)) ≈ vcmax_ref
end


@testset "max_electron_transport_rate/JMAXTFN" begin
    jmaxtfn_fortran = Libdl.dlsym(photosynlib, :jmaxtfn_)
    f = Jmax()
    jmax25 = ustrip(u"μmol*m^-2*s^-1", f.jmax25)
    eavj = ustrip(u"J*mol^-1", f.eavj)
    edvj = ustrip(u"J*mol^-1", f.edvj)
    delsj = ustrip(u"J*K^-1*mol^-1", f.delsj)
    tvjup = -100.0
    tvjdn = -100.0
    jmax_ref = ccall(jmaxtfn_fortran, Float32, 
                     (Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}), 
                     jmax25, tleaf, eavj, edvj, delsj, tvjup, tvjdn)
    @test ustrip(max_electron_transport_rate(f, v.tleaf)) ≈ jmax_ref
end



@testset "respiration/RESP" begin
    resp_fortran = Libdl.dlsym(photosynlib, :resp_)
    f = Respiration()
    rd0 = ustrip(u"μmol*m^-2*s^-1", f.rd0)
    rdacc = 1.0
    q10f = ustrip(u"K^-1", f.q10f)
    tref = ustrip(°C, f.tref)
    dayresp = f.dayresp
    tbelow = ustrip(°C, f.tbelow)
    k10f = 0.0 # No acclimation 
    tmove = 0.0 # No acclimation 
    resp_ref = ccall(resp_fortran, Float32,
                     (Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}),
                     rd0, rdacc, tleaf, tmove, q10f, k10f, tref, dayresp, tbelow)
    v.rd = respiration(f, v.tleaf)
    @test ustrip(u"μmol*m^-2*s^-1", v.rd) ≈ resp_ref

    f = AcclimatizedRespiration()
    rd0 = ustrip(u"μmol*m^-2*s^-1", f.rd0)
    rdacc = 1.0 # this isn't actually a parameter
    q10f = ustrip(u"K^-1", f.q10f)
    tref = ustrip(°C, f.tref)
    dayresp = f.dayresp
    tbelow = ustrip(°C, f.tbelow)
    k10f = ustrip(u"K^-1", f.k10f)
    tmove = ustrip(u"K", f.tmove)
    resp_ref = ccall(resp_fortran, Float32,
                     (Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}),
                     rd0, rdacc, tleaf, tmove, q10f, k10f, tref, dayresp, tbelow)
    v.rd = respiration(f, v.tleaf)

    # This is actually commented out in the maespa FORTRAN
    # @test v.rd.val ≈ resp_ref
end

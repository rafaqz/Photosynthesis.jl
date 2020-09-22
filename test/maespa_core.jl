using Photosynthesis

include(joinpath(dirname(pathof(Photosynthesis)), "../test/shared.jl"))

# Setup
emax = MaespaEnergyBalance(
    photosynthesis_model=FvCBPhotosynthesis(
        stomatal_conductance_model=EmaxStomatalConductance()
    )
)
ph = emax.photosynthesis_model

@testset "rubisco_compensation_point/kmfn" begin
    for tleaf in (-50, -20.0, 0.0, 15.0, 25.0, 45.0, 60.0, 90.0)
        v = EmaxVars()
        v.tleaf = tleaf * °C
        kmfn_fortran = Libdl.dlsym(maespa_photosynlib, :kmfn_)
        ieco = BERNACCI
        km_ref = ccall(kmfn_fortran, Float32, (Ref{Float32}, Ref{Int32}), tleaf, ieco)
        km = rubisco_compensation_point(BernacchiCompensation(), v.tleaf) # Michaelis-Menten for Rubisco, umol mol-1
        @test ustrip(u"μmol/mol", km) ≈ km_ref
        println("Bernacchi: ", (v.tleaf, km))
        ieco = BADGERCOLLATZ
        km_ref = ccall(kmfn_fortran, Float32, (Ref{Float32}, Ref{Int32}), tleaf, ieco)
        km = rubisco_compensation_point(BadgerCollatzCompensation(), v.tleaf) # Michaelis-Menten for Rubisco, umol mol-1
        @test ustrip(u"μmol/mol", km) ≈ km_ref
        println("Badger-Collatz: ", (v.tleaf, km))
    end
end


@testset "co2_compensation_point/GAMMAFN" begin
    for tleaf in (-60, -20.0, 0.0, 15.0, 25.0, 45.0, 60.0, 100.0)
        v = EmaxVars()
        v.tleaf = tleaf * °C
        gammafn_fortran = Libdl.dlsym(maespa_photosynlib, :gammafn_)
        gammastar_ref = ccall(gammafn_fortran, Float32, (Ref{Float32}, Ref{Int32}), tleaf, BERNACCI)
        gammastar = co2_compensation_point(BernacchiCompensation(), v.tleaf) # Michaelis-Menten for Rubisco, umol mol-1
        @test ustrip(u"μmol/mol", gammastar) ≈ gammastar_ref
        println("Bernacchi: ", (v.tleaf, gammastar))
        gammastar_ref = ccall(gammafn_fortran, Float32, (Ref{Float32}, Ref{Int32}), tleaf, BADGERCOLLATZ)
        gammastar = co2_compensation_point(BadgerCollatzCompensation(), v.tleaf) # Michaelis-Menten for Rubisco, umol mol-1
        @test ustrip(u"μmol/mol", gammastar) ≈ gammastar_ref
        println("Badger-Collatz: ", (v.tleaf, gammastar))
    end
end


@testset "max_rubisco_activity/VCMAXTFN" begin
    vcmaxtfn_fortran = Libdl.dlsym(maespa_photosynlib, :vcmaxtfn_)

    for parmult in (0.001, 0.01, 0.5, 1.0, 1.03)
        # Breaks below 0.0, but only by small amounts
        for tleaf in (0.0, 10.0, 15.0, 25.0, 50.0)
            v = EmaxVars()
            v.tleaf = tleaf * °C
            vc = Flatten.modify(x -> x * rand(min(parmult, 1.0):0.000001:max(parmult)), NoOptimumVcmax())
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
        end

        for tleaf in (-40.0, 0.0, 10.0, 15.0, 25.0, 55.0)
            v = EmaxVars()
            v.tleaf = tleaf * °C

            # The FORTRAN breaks with parmult > 1.1
            vc = Flatten.modify(x -> x * parmult, OptimumVcmax())
            vcmax25 = ustrip(u"μmol*m^-2*s^-1", vc.vcmax25)
            eavc = ustrip(u"J*mol^-1", vc.eavc)
            edvc = ustrip(u"J*mol^-1", vc.edvc)
            delsc = ustrip(u"J*K^-1*mol^-1", vc.delsc)
            tvjup = -100.0
            tvjdn = -100.0
            vcmax_ref = ccall(vcmaxtfn_fortran, Float32,
                              (Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}
                              ), vcmax25, tleaf, eavc, edvc, delsc, tvjup, tvjdn)
            @test ustrip(max_rubisco_activity(vc, v.tleaf)) ≈ vcmax_ref

            # vc = OptimumVcmax()
            # f = DukeFlux()
            # vcmax25 = ustrip(u"μmol*m^-2*s^-1", vc.vcmax25)
            # eavc = ustrip(u"J*mol^-1", vc.eavc)
            # edvc = ustrip(u"J*mol^-1", vc.edvc)
            # delsc = ustrip(u"J*K^-1*mol^-1", vc.delsc)
            # tvjup = ustrip(°C, f.tvjup)
            # tvjdn = ustrip(°C, f.tvjdn)
            # vcmax_ref = ccall(vcmaxtfn_fortran, Float32,
            #                   (Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}),
            #                   vcmax25, tleaf, eavc, edvc, delsc, tvjup, tvjdn)
            # @test ustrip(max_rubisco_activity(vc, v.tleaf)) ≈ vcmax_ref

        end
    end
end


@testset "max_electron_transport_rate/JMAXTFN" begin
    params = Flatten.flatten(Jmax(), Number)
    fns = Flatten.fieldnameflatten(Jmax(), Number)
    for parmult in permutations((0.01, 0.1, 1.0, 1.04))
        for tleaf in (-50.0, -10, 0.0, 15.0, 25.0, 50.0)
            v = EmaxVars()
            v.tleaf = tleaf * °C
            ps = NamedTuple{fns}(params .* parmult)
            @show ps
            f = Flatten.reconstruct(Jmax(), ps, Number)
            jmaxtfn_fortran = Libdl.dlsym(maespa_photosynlib, :jmaxtfn_)
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
    end
end


@testset "respiration/RESP" begin
    resp_fortran = Libdl.dlsym(maespa_photosynlib, :resp_)

    params = Flatten.flatten(Respiration(), Number)
    fns = Flatten.fieldnameflatten(Respiration(), Number)
    for parmult in permutations((0.7, 0.9, 1.0, 1.05, 1.15))
        for tleaf in (-20.0, 0.0, 15.0, 25.0, 50.0)
            v = EmaxVars()
            v.tleaf = tleaf * °C
            ps = NamedTuple{fns}(params .* parmult)
            @show ps
            f = Flatten.reconstruct(Respiration(), ps, Number)
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
            println(v.rd)
            @test ustrip(u"μmol*m^-2*s^-1", v.rd) ≈ resp_ref
        end
    end

    # params = Flatten.flatten(AcclimatizedRespiration(), Number)
    # fns = Flatten.fieldnameflatten(AcclimatizedRespiration(), Number)
    # for parmult in permutations((0.20, 0.5, 0.8, 1.0, 1.1, 1.5, 2.0))
    #     ps = NamedTuple{fns}(params .* parmult)
    #     @show ps
    #     f = Flatten.reconstruct(AcclimatizedRespiration(), ps, Number)
    #     rd0 = ustrip(u"μmol*m^-2*s^-1", f.rd0)
    #     rdacc = 1.0 # this isn't actually a parameter
    #     q10f = ustrip(u"K^-1", f.q10f)
    #     tref = ustrip(°C, f.tref)
    #     dayresp = f.dayresp
    #     tbelow = ustrip(°C, f.tbelow)
    #     k10f = ustrip(u"K^-1", f.k10f)
    #     tmove = ustrip(u"K", f.tmove)
    #     resp_ref = ccall(resp_fortran, Float32,
    #                      (Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}),
    #                      rd0, rdacc, tleaf, tmove, q10f, k10f, tref, dayresp, tbelow)
    #     v.rd = respiration(f, v.tleaf)
    #     # This is actually commented out in the maespa FORTRAN
    #     # @test v.rd.val ≈ resp_ref
    # end

end

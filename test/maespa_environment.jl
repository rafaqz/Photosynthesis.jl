using Photosynthesis

include(joinpath(dirname(pathof(Photosynthesis)), "../test/shared.jl"))

@testset "boundary_conductance_forced/GBHFORCED" begin
    gbhforced_fortran = Libdl.dlsym(maespa_photosynlib, :gbhforced_)
    v = BallBerryVars()
    for tair in (-50.0, 0.0, 15.0, 28.0, 40.0, 60.0, 100.0),
        press in (80000.0, 101250.0, 150000.0),
        wind in (0.1, 1.0, 10.0),
        wleaf in (0.01, 0.05, 0.1, 1.0)
        v.windspeed = wind * u"m*s^-1"
        v.tair = tair * u"°C"
        v.pressure = press * u"Pa"
        bc_model = BoundaryConductance(wleaf * u"m")
        bcforced_ref = ccall(gbhforced_fortran, Float32, 
                             (Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}), 
                             tair, press, wind, wleaf)
        bcforced = boundary_conductance_forced(bc_model, v)
        upreferred(bcforced)
        @test ustrip(u"mol*m^-2*s^-1", bcforced) ≈ bcforced_ref
    end
end

@testset "boundary_conductance_free/GBHFREE" begin
    v = BallBerryVars()
    gbhfree_fortran = Libdl.dlsym(maespa_photosynlib, :gbhfree_)
    for tleaf in (-50.0, 0.0, 15.0, 28.0, 40.0, 60.0, 100.0),
        tair in (-50.0, 0.0, 15.0, 28.0, 40.0, 60.0, 100.0),
        press in (80000.0, 101250.0, 150000.0),
        wleaf in (0.01, 0.05, 0.1, 1.0)
        v.tleaf = tleaf * °C |> K
        v.tair = tair * °C |> K
        v.pressure = press * u"Pa"
        bc_model = BoundaryConductance(wleaf * u"m")
        wleaf = ustrip(u"m", bc_model.leafwidth)
        bcfree_ref = ccall(gbhfree_fortran, Float32, 
                           (Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}), 
                           tair, tleaf, press, wleaf)
        bcfree = boundary_conductance_free(bc_model, v)
        @test ustrip(u"mol*m^-2*s^-1", bcfree) ≈ bcfree_ref
    end
end

@testset "radiation_conductance/GRADIATION" begin
    v = BallBerryVars()
    for tair in (-50.0, 0.0, 15.0, 28.0, 40.0, 60.0, 100.0)
        rc_model = WangRadiationConductance()
        gradiation_fortran = Libdl.dlsym(maespa_photosynlib, :gradiation_)
        tair = ustrip(°C, v.tair)
        rdfipt = rc_model.rdfipt
        tuipt = rc_model.tuipt
        tdipt = rc_model.tdipt
        rc_ref = ccall(gradiation_fortran, Float32, 
                       (Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}), 
                       tair, rdfipt, tuipt, tdipt)
        rc = radiation_conductance(rc_model, v)
        @test ustrip(u"mol*m^-2*s^-1", rc) ≈ rc_ref
    end
end

@testset "saturated_vapour_pressure/SATUR" begin
    for tair in (-50.0, 0.0, 15.0, 28.0, 40.0, 100.0)
        satur_fortran = Libdl.dlsym(maespa_photosynlib, :satur_)
        svp_ref = ccall(satur_fortran, Float32, (Ref{Float32},), tair)
        svp = saturated_vapour_pressure(tair * °C)
        @test ustrip(u"Pa", svp) ≈ svp_ref
    end
end

@testset "latent_heat_water_vapour" begin
    @test ustrip(u"kJ/mol", latent_heat_water_vapour(40°C)) ≈ 43.35 rtol=1e-3
    @test ustrip(u"kJ/mol", latent_heat_water_vapour(25°C)) ≈ 43.99 rtol=1e-3
    @test ustrip(u"kJ/mol", latent_heat_water_vapour(0°C)) ≈ 45.06 rtol=1e-3
end

@testset "PenmanMonteithEvapotranspiration/PENMON" begin
    for press in (80000.0, 101250.0, 150000.0),
        vpd in (10.0, 100.0, 1000.0, 10000.0),
        lhv in (1000.0, 10000.0, 50000.0),
        gv in (0.0, 1.0, 10),
        gh in (0.0, 1.0, 10),
        rnet in (0.0, 10.0, 1000.0),
        slope in (10.0, 100.0, 1.000)
        penmon_fortran = Libdl.dlsym(maespa_photosynlib, :penmon_)
        p = MaespaEnergyBalance(
            photosynthesis_model=FvCBPhotosynthesis(
                stomatal_conductance_model=BallBerryStomatalConductance(
                    gs_submodel=MedlynStomatalConductanceSubModel(),
                    soil_model=NoSoilMethod()),
                flux_model=DukeFlux(),
                compensation_model=BernacchiCompensation(),
                respiration_model=Respiration(),
            ),
            evapotranspiration_model = PenmanMonteithEvapotranspiration(),
            max_iterations=0
        )
        # Prime the variables to reasonable values
        v = BallBerryVars()
        et_model = p.evapotranspiration_model
        enbal!(v, p)
        v.pressure = press * u"Pa"
        v.rnet = rnet * u"J*m^-2*s^-1" 
        v.vpd = vpd * u"Pa"
        v.gv = gv * u"mol*m^-2*s^-1"
        v.gh = gh * u"mol*m^-2*s^-1"
        v.lhv = lhv * u"J*mol^-1"
        v.slope = slope * u"Pa*K^-1"
        et_ref = ccall(
             penmon_fortran, Float32, 
             (Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}), 
             press, slope, lhv, rnet, vpd, gh, gv
        )
        et = evapotranspiration(et_model, v)
        @test ustrip(u"kg*mol*J^-1*s^-3", et) ≈ et_ref
    end
end


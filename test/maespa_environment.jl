using Photosynthesis

include(joinpath(dirname(pathof(Photosynthesis)), "../test/shared.jl"))

@testset "boundary_conductance_forced/GBHFORCED" begin
    v = BallBerryVars()

    gbhforced_fortran = Libdl.dlsym(photosynlib, :gbhforced_)

    bc_model = BoundaryConductance()

    tair = ustrip(u"°C", v.tair)
    press = ustrip(u"Pa", v.pressure)
    wind = ustrip(u"m*s^-1", v.windspeed)
    wleaf = ustrip(u"m", bc_model.leafwidth)

    bcforced_ref = ccall(gbhforced_fortran, Float32, 
                         (Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}), 
                         tair, press, wind, wleaf)
    bcforced = boundary_conductance_forced(bc_model, v)
    upreferred(bcforced)
    
    @test ustrip(u"mol*m^-2*s^-1", bcforced) ≈ bcforced_ref
end

@testset "boundary_conductance_free/GBHFREE" begin
    v = BallBerryVars()
    v.tleaf = 28.0°C
    bc_model = BoundaryConductance()
    gbhfree_fortran = Libdl.dlsym(photosynlib, :gbhfree_)

    tair = ustrip(u"°C", v.tair)
    tleaf = ustrip(u"°C", v.tleaf)
    press = ustrip(u"Pa", v.pressure)
    wleaf = ustrip(u"m", bc_model.leafwidth)

    bcfree_ref = ccall(gbhfree_fortran, Float32, 
                       (Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}), 
                       tair, tleaf, press, wleaf)
    bcfree = boundary_conductance_free(bc_model, v)
    
    @test ustrip(u"mol*m^-2*s^-1", bcfree) ≈ bcfree_ref
end

@testset "radiation_conductance/GRADIATION" begin
    v = BallBerryVars()
    rc_model = WangRadiationConductance()
    gradiation_fortran = Libdl.dlsym(photosynlib, :gradiation_)

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

@testset "saturated_vapour_pressure/SATUR" begin
    v = BallBerryVars()
    satur_fortran = Libdl.dlsym(photosynlib, :satur_)
    svp_ref = ccall(satur_fortran, Float32, (Ref{Float32},), 15.0)
    svp = saturated_vapour_pressure(15°C)
    @test ustrip(u"Pa", svp) ≈ svp_ref
end


@testset "latent_heat_water_vapour" begin
    @test ustrip(u"kJ/mol", latent_heat_water_vapour(40°C)) ≈ 43.35 rtol=1e-3
    @test ustrip(u"kJ/mol", latent_heat_water_vapour(25°C)) ≈ 43.99 rtol=1e-3
    @test ustrip(u"kJ/mol", latent_heat_water_vapour(0°C)) ≈ 45.06 rtol=1e-3
end

@testset "PenmanMonteithEvapotranspiration/PENMON" begin
    v = BallBerryVars()
    penmon_fortran = Libdl.dlsym(photosynlib, :penmon_)
    p = FvCBEnergyBalance(
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
    et_model = p.evapotranspiration_model
    # Prime the variables to reasonable values
    enbal!(v, p)

    press = ustrip(u"Pa", v.pressure)
    slope_ = ustrip(u"Pa*K^-1", v.slope)
    rnet = ustrip(u"J*m^-2*s^-1", v.rnet)
    vpd = ustrip(u"Pa", v.vpd)
    lhv =  ustrip(u"J*mol^-1", v.lhv)
    gh = ustrip(u"mol*m^-2*s^-1", v.gh)
    gv = ustrip(u"mol*m^-2*s^-1", v.gv)

    et_ref = ccall(
         penmon_fortran, Float32, 
         (Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}), 
         press, slope_, lhv, rnet, vpd, gh, gv
    )
    et = evapotranspiration(et_model, v)

    @test ustrip(u"kg*mol*J^-1*s^-3", et) ≈ et_ref
end

using SafeTestsets

@time @safetestset "test all constructors work" begin include("construct.jl") end
@time @safetestset "test against maespa core" begin include("maespa_core.jl") end
@time @safetestset "test stomatal conductance against maespa fortran" begin include("maespa_stomata.jl") end
@time @safetestset "test against maespa fortran" begin include("maespa_photosyn.jl") end
@time @safetestset "test against maespa fortran" begin include("maespa_enbal.jl") end
@time @safetestset "test against maespa fortran" begin include("maespa_utils.jl") end

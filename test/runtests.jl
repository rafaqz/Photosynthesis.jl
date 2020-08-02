using SafeTestsets

@time @safetestset "test all constructors work" begin include("construct.jl") end
@time @safetestset "test core functions against maespa core" begin include("maespa_core.jl") end
@time @safetestset "test environmental functions against maespa core" begin include("maespa_environment.jl") end
@time @safetestset "test stomatal conductance against maespa fortran" begin include("maespa_stomata.jl") end
@time @safetestset "test photosynthesis against maespa fortran" begin include("maespa_photosyn.jl") end
@time @safetestset "test energy balance against maespa fortran" begin include("maespa_enbal.jl") end
@time @safetestset "test utils against maespa fortran" begin include("maespa_utils.jl") end

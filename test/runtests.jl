using SafeTestsets

@time @safetestset "test all constructors work" begin include("construct.jl") end
@time @safetestset "test against maespa fortran" begin include("maespa_fortran.jl") end

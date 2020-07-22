
"""
Canopy-atmospheric decoupling models, 
calculated in [`decoupling`](@ref) method.
"""
abstract type AbstractDecoupling end

"""
    decoupling(f::AbstractDecoupling, v)

Calculate decoupling, returning a float between 0.0 and 1.0.
"""
function decoupling end

"""
    McNaughtonJarvisDecoupling()

Calculate decoupling coefficient following McNaughton and Jarvis 1986
"""
struct McNaughtonJarvisDecoupling <: AbstractDecoupling end

decoupling(f::McNaughtonJarvisDecoupling, v) = begin
    γc = CPAIR * AIRMA * v.pressure / v.lhv
    epsilon = ustrip(v.slope / γc) # TODO why is ustrip needed here?
    (1.0 + epsilon) / (1.0 + epsilon + v.gbv / v.gsv)
end

"""
    NoDecoupling()

Don't calculate decoupling.
"""
struct NoDecoupling <: AbstractDecoupling end

decoupling(f::NoDecoupling, v) = 0.0

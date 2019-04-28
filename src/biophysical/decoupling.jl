abstract type AbstractDecoupling end

struct McNaughtonJarvisDecoupling <: AbstractDecoupling end
struct NoDecoupling <: AbstractDecoupling end

"""
    calc_decoupling(f::McNaughtonJarvisDecoupling, v)
Calculate decoupling coefficient following McNaughton and Jarvis 1986
"""
calc_decoupling(f::McNaughtonJarvisDecoupling, v) = begin
    γc = CPAIR * AIRMA * v.pressure / v.lhv
    epsilon = ustrip(v.slope / γc) # TODO why is ustrip needed here?
    (1.0 + epsilon) / (1.0 + epsilon + v.gbv / v.gsv)
end

"""
    calc_decoupling(f::NoDecoupling, v)
Don't calculate decoupling.
"""
calc_decoupling(f::NoDecoupling, v) = 0.0

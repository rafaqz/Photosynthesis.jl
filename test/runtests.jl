using Photosynthesis
using Unitful

@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

include("construct.jl")
include("maespa.jl")

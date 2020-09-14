
"""
Energy balance models calculate leaf temperature, usually also running 
photosynthesis and all other model components along with environmental models
like radiation and boundary layer conductance.

They are run in [`enbal!`](@ref) methods.
"""
abstract type AbstractEnergyBalance end

"""
    enbal!(v, m::AbstractEnergyBalance)

Calculates leaf photosynthesis and transpiration for an 
[`AbstractEnergyBalance`](@ref) model `m` and variables `v`.

Results are written to `v`.
"""
function enbal! end

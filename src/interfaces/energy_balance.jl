
"""
Energy balance models
"""
abstract type AbstractEnergyBalance end

"""
    model_init!(v, f)

Runs any model initialisation that needs to happen at the start of energy balance
"""
function enbal_init! end

"""
    enbal!(v, m)

Calculates leaf photosynthesis and transpiration for an `AbstractEnergyBalance`
model `m` and variables `v`.

Results are written to v.
"""
function enbal! end

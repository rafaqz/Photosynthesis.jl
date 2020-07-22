# Interface for stomatal conductance

"""
Stomatal conductance models
"""
abstract type AbstractStomatalConductance end

"""
Stomatal conductance submodels
"""
abstract type AbstractStomatalConductanceSubModel end

"""
    stomatal_conductance!(v, p)
Stomatal conductance and intercellular CO2 partial pressure calculations.

v.aleaf is NET leaf photosynthesis.
"""
function stomatal_conductance! end

"""
    photo_init!(vars, params::AbstractStomatalConductance)

Initialise variables based on the AbstractStomatalConductance model. 
"""
function gs_init! end

"""
    photo_update!(vars, params::AbstractStomatalConductance)

Update variables based on the model. 
"""
function gs_update! end 

"""
    gs_div_a(m::AbstractStomatalConductanceSubModel, v)

Returns the value of stomatal conductance `gs` divided by 
assimilation `a` for sub-model `m` given variables `v`.
"""
function gs_div_a end

"""
    update_extremes!(v, m::AbstractStomatalConductance)

Update variables in extreme conditions.
"""
function update_extremes! end

"""
    transport_limited_rate(m::AbstractStomatalConductance, v, gs_div_a)

Transport limited rate of assimilation.
"""
function transport_limited_rate end

"""
    rubisco_limited_rate(m::AbstractStomatalConductance, v, gs_div_a)

Rubisco limited rate of assimilation.
"""
function rubisco_limited_rate end

# Interface for stomatal conductance

abstract type AbstractStomatalConductance end

abstract type AbstractStomatalConductanceSubModel end

"""
    stomatal_conductance!(p, v)
Stomatal conductance and intercellular CO2 partial pressure calculations.

v.aleaf is NET leaf photosynthesis.
"""
function stomatal_conductance! end

function rubisco_limited_rate end

function transport_limited_rate end

"""
    gsdiva(p, v)
Formulation-specific component for the Ball-Berry family of stomatal conductance models.
"""
function gsdiva end

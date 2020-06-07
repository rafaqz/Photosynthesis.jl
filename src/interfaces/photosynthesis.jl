"""
Abstract supertype for all photosynthesis models
"""
abstract type AbstractPhotosynthesis end

"""
    photosynthesis!(params::AbstractPhotosynthesis, vars)

Run a photosynthesis model.
"""
function photosynthesis! end

"""
    photo_init!(vars, params::AbstractStomatalConductance)

Initialise variables based on AbstractStomatalConductance model. 
"""
function photo_init! end

function photo_update! end 

function check_extremes! end

function update_extremes! end

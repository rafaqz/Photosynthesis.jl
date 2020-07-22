"""
Abstract supertype for all photosynthesis models
"""
abstract type AbstractPhotosynthesis end

"""
    photosynthesis!(vars, params::AbstractPhotosynthesis)

Run a photosynthesis model, writing results to `vars`.
"""
function photosynthesis! end

"""
    check_extremes!(v, p::AbstractFvCBPhotosynthesis)

Check extreme values are in tolerable ranges
"""
function check_extremes! end

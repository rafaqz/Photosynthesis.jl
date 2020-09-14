"""
Abstract supertype for all photosynthesis models.

These expect to run inside [`enbal!`](@ref), and otherwise
nead to have the `tleaf` and `vpdleaf` variables manually set.
"""
abstract type AbstractPhotosynthesis end

"""
    photosynthesis!(vars, params::AbstractPhotosynthesis)

Run a photosynthesis model, writing results to `vars`.
"""
function photosynthesis! end

"""
    check_extremes!(v, p::AbstractFvCBPhotosynthesis)

Check extreme values are in tolerable ranges.
"""
function check_extremes! end

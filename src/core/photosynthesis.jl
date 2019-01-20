

abstract type AbstractPhotoModel end

"""
Fixed parameters for photosynthesis.

THe majority of these are "composed" from submodels that both hold
the parameters and use their own specific methods during the photosynthesis
routines.

Calling `PhotoParams()` will give the default values for all of these submodels.
Any parameters and submodels can be overridden with keyword arguments:

`PhotoParams(model=TuzetModel, ca= 450 | μmol*mol^-1)`
"""
@default_kw struct FvCBPhoto{F<:AbstractPhotoModel,
                             V<:AbstractVcJmax,
                             KM<:AbstractCompensation,
                             Ru<:AbstractRubiscoRegen,
                             Re<:Union{Nothing,AbstractRespiration}}
    model::F          | BallBerryModel()       
    vcjmax::V         | VcJmax()               
    compensation::KM  | BernacchiCompensation()
    rubisco_regen::Ru | RubiscoRegen()         
    respiration::Re   | nothing          
end

"""
    photosynthesis!(v, p)
Calculates photosynthesis according to the ECOCRAFT
agreed formulation of the Farquharvon Caemmerer (1982) equations.

Farquhar, G.D., S. Caemmerer and J.A. Berry. 1980. 
A biochemical model of photosynthetic CO2 assimilation in leaves of C3 species. 
Planta. 149:78-90. 
"""
function photosynthesis!(v, p)
    v.gammastar = co2_compensation_point(p.compensation, v, p) # CO2 compensation point, umol mol-1
    v.km = rubisco_compensation_point(p.compensation, v, p) # Michaelis-Menten for Rubisco, umol mol-1
    v.jmax = max_electron_transport_rate(p.vcjmax, v, p) # Potential electron transport rate, umol m-2 s-1
    v.vcmax = max_rubisco_activity(p.vcjmax, v, p) # Maximum Rubisco activity, umol m-2 s-1

    v.rd = respiration(p.respiration, v, p, v.rd) # Day leaf respiration, umol m-2 s-1

    soilmoisture_conductance!(p.model.soilmethod, v, p)

    v.vj = rubisco_regeneration(p.rubisco_regen, v, p)

    extremes!(v, p) && return nothing
    stomatal_conductance!(p.model, v, p)

    nothing
end


"Deal with extremes where jmax or vcmax are less than zero"
function extremes! end

extremes!(v, p) = begin
    if v.jmax <= zero(v.jmax) || v.vcmax <= zero(v.vcmax)
        println("extreme")
        extremes!(p.model, v)
        return true
    end
    false
end

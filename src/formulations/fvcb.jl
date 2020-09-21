
"""
Abstract supertype for all Farquhar/von Caemmerer/Berry derived photosynthesis models
"""
abstract type AbstractFvCBPhotosynthesis <: AbstractPhotosynthesis end

flux_model(p::AbstractFvCBPhotosynthesis) = p.flux_model
compensation_model(p::AbstractFvCBPhotosynthesis) = p.compensation_model
rubisco_regen_model(p::AbstractFvCBPhotosynthesis) = p.rubisco_regen_model
respiration_model(p::AbstractFvCBPhotosynthesis) = p.respiration_model
stomatal_conductance_model(p::AbstractFvCBPhotosynthesis) = p.stomatal_conductance_model

"""
    FvCBPhotosynthesis(flux, compensation, rubisco_regen, respiration, stomatal_conductance)

General Farquhar von Caemmerer Berry model of photosynthesis. Organised as a modular
components of electron flux, CO2 and rubusco compensation, rubisco regeneration, 
respiration and stomatal conductance.

Calculates photosynthesis according to the ECOCRAFT
agreed formulation of the Farquharvon Caemmerer (1982) equations.

Farquhar, G.D., S. Caemmerer and J.A. Berry. 1980.
A biochemical model of photosynthetic CO2 assimilation in leaves of C3 species.
Planta. 149:78-90.

$(FIELDDOCTABLE)
"""
@default_kw struct FvCBPhotosynthesis{F,KM,Ru,Re,GS} <: AbstractFvCBPhotosynthesis
    flux_model::F                  | Flux()
    compensation_model::KM         | BernacchiCompensation()
    rubisco_regen_model::Ru        | RubiscoRegen()
    respiration_model::Re          | Respiration()
    stomatal_conductance_model::GS | BallBerryStomatalConductance()
end

check_extremes!(v, p::AbstractFvCBPhotosynthesis) =
    if v.jmax <= zero(v.jmax) || v.vcmax <= zero(v.vcmax)
        update_extremes!(v, stomatal_conductance_model(p))
        true
    else
        false
    end

update_extremes!(v, m::AbstractStomatalConductance) = begin
    v.aleaf = -v.rd
    v.gs = g0(m)
end

function photosynthesis!(v, m::AbstractFvCBPhotosynthesis)
    v.gammastar = co2_compensation_point(compensation_model(m), v.tleaf) # CO2 compensation point, umol mol-1
    v.km = rubisco_compensation_point(compensation_model(m), v.tleaf) # Michaelis-Menten for Rubisco, umol mol-1
    v.jmax, v.vcmax = flux(flux_model(m), v) # Potential electron transport rate and maximum Rubisco activity
    v.rd = respiration(respiration_model(m), v.tleaf) # Day leaf respiration, umol m-2 s-1
    v.vj = rubisco_regeneration(rubisco_regen_model(m), v)

    # Zero values and exit in extreme cases.
    check_extremes!(v, m) && return

    v.aleaf, v.gs = stomatal_conductance!(v, stomatal_conductance_model(m))

    v.ci = if v.gs > zero(v.gs) && v.aleaf > zero(v.aleaf) 
        v.cs - v.aleaf / v.gs
    else
        v.cs
    end
    return
end

"""
    @MixinEnviroVars

Mixin variables for [`AbstractMaespaEnergyBalance`](@ref) variables objects.
"""
@mix @vars struct MixinFvCBVars{GSDA,KM,CI,GAM,GS,JMX,VCMX,RD,AC,AL,VJ,AJ,FS}
    # photosynthesis
    gs_div_a::GSDA   | 0.0    | mol*μmol^-1    | _
    km::KM           | 0.0    | μmol*mol^-1    | _
    ci::CI           | 0.0    | μmol*mol^-1    | _
    gammastar::GAM   | 0.0    | μmol*mol^-1    | _
    gs::GS           | 0.0    | mol*m^-2*s^-1  | _
    jmax::JMX        | 0.0    | μmol*m^-2*s^-1 | _
    vcmax::VCMX      | 0.0    | μmol*m^-2*s^-1 | _
    rd::RD           | 0.0    | μmol*m^-2*s^-1 | _
    ac::AC           | 0.0    | μmol*m^-2*s^-1 | _
    aleaf::AL        | 0.0    | μmol*m^-2*s^-1 | _
    vj::VJ           | 0.0    | μmol*m^-2*s^-1 | _
    aj::AJ           | 0.0    | μmol*m^-2*s^-1 | _
    # soil
    fsoil::FS        | 0.0    | _              | _
end



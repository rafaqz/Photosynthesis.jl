@MaespaSoil struct EmaxSoilMethod{} <: AbstractSoilMethod{T} end

"""
Emax stomatal conductance model from maespa.
Similar options for specialised stomatal conducance,
but using emax soil water methods.
"""
@Maespa struct EmaxModel{GS<:AbstractStomatalConductance,
                         SM<:EmaxSoilMethod{<:AbstractPotentialDependence}
                        } <: AbstractMaespaModel
    gsmodel::GS    | BallBerryStomatalConductance() | _ | _ | _ | _
    soilmethod::SM | EmaxSoilMethod()               | _ | _ | _ | _
end

@MaespaVars mutable struct EmaxVars{kPa,mMoM2S}
    minleafwp::kPa   | 0.1         | kPa                   | _
    emaxleaf::mMoM2S | 400.0       | mmol*m^-2*s^-1        | _
end

soilmoisture_conductance!(f::EmaxSoilMethod, v) = begin
    pd = potential_dependence(f.soildata, v.weightedswp)
    v.vcmax *= pd
    v.jmax *= pd
    v.fsoil = oneunit(v.fsoil)
end

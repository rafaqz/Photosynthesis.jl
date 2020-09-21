module Photosynthesis
@doc let
    path = joinpath(dirname(@__DIR__), "README.md")
    include_dependency(path)
    read(path, String)
end Photosynthesis

using DocStringExtensions,
      FieldDefaults,
      FieldDocTables,
      FieldMetadata,
      Mixers,
      SimpleRoots,
      Unitful

using Unitful: R, °C, K, Pa, kPa, MPa, J, W, kJ, kg, g, m, s, mol, mmol, μmol, σ

import FieldMetadata: @flattenable, @bounds, @default, @description, @units,
                      flattenable, bounds, default, description, units

@metadata units NoUnits


export enbal!,
       photosynthesis!,
       stomatal_conductance!,
       soil_water_conductance!,
       co2_compensation_point,
       rubisco_compensation_point,
       rubisco_regeneration,
       rubisco_limited_rate,
       transport_limited_rate,
       slope,
       max_electron_transport_rate,
       max_rubisco_activity,
       decoupling, 
       leaftemp,
       cmolar,
       shape_gs,
       respiration,
       factor_conductance,
       radiation_conductance,
       latent_heat_water_vapour,
       saturated_vapour_pressure,
       vapour_pressure_deficit,
       wang_radiation_conductance,
       boundary_conductance_free,
       boundary_conductance_forced,
       forced_boundary_conductance,
       arrhenius,
       grashof_number,
       gs_div_a,
       evapotranspiration,
       penman_monteith_evapotranspiration

export Compensation, BadgerCollatzCompensation, BernacchiCompensation

export AbstractJmax, Jmax

export AbstractVcmax, NoOptimumVcmax, OptimumVcmax

export AbstractFlux, Flux, DukeFlux, PotentialModifiedFlux

export AbstractRubiscoRegen, RubiscoRegen

export AbstractRespiration, Respiration, AcclimatizedRespiration

export StomatalConductanceShape, HardMinimum, HyperbolicMinimum

export AbstractRadiationConductance, WangRadiationConductance

export AbstractBoundaryConductance, BoundaryConductance

export AbstractDecoupling, McNaughtonJarvisDecoupling, NoDecoupling

export AbstractEvapotranspiration, PenmanMonteithEvapotranspiration

export AbstractSoilData, AbstractDeficitSoilData, AbstractContentSoilData,
       DeficitSoilData, ContentSoilData, SimulatedSoilData, PotentialSoilData, NoSoilData

export AbstractPotentialDependence, LinearPotentialDependence, 
       ZhouPotentialDependence, NoPotentialDependence

export AbstractSoilMethod, NoSoilMethod, VolumetricSoilMethod, ConstantSoilMethod,
       DeficitSoilMethod, PotentialSoilMethod, EmaxSoilMethod, TuzetSoilMethod


export AbstractJarvisCO2, JarvisNoCO2, JarvisLinearCO2, JarvisNonlinearCO2

export AbstractJarvisVPD, JarvisHyperbolicVPD, JarvisLohammerVPD,
       JarvisFractionDeficitVPD, JarvisLinearDeclineVPD

export AbstractJarvisLight, JarvisLight

export AbstractJarvisTemp, JarvisNoTemp, JarvisTemp1, JarvisTemp2


export AbstractStomatalConductanceSubModel, AbstractBallBerryStomatalConductanceSubModel, 
       BallBerryStomatalConductanceSubModel, LeuningStomatalConductanceSubModel, 
       MedlynStomatalConductanceSubModel, TuzetStomatalConductanceSubModel

export AbstractStomatalConductance, AbstractBallBerryStomatalConductance, 
       BallBerryStomatalConductance, TuzetStomatalConductance, 
       AbstractEmaxStomatalConductance, EmaxStomatalConductance, JarvisStomatalConductance

export AbstractPhotosynthesis, AbstractFvCBPhotosynthesis, FvCBPhotosynthesis

export AbstractEnergyBalance, AbstractMaespaEnergyBalance, MaespaEnergyBalance, 
       EmaxEnergyBalance, TuzetEnergyBalance

export BallBerryVars, EmaxVars, TuzetVars, JarvisVars


const FIELDDOCTABLE = FieldDocTable((:Default, :Units, :Description),
                                    (default, units, description);
                                    truncation=(40, 40, 100))

@chain columns @udefault_kw @units @bounds @description

include("constants.jl")
include("vars.jl")

include("biophysical/biophysical.jl")
include("biophysical/boundary_conductance.jl")
include("biophysical/decoupling.jl")
include("biophysical/evapotranspiration.jl")
include("biophysical/radiation_conductance.jl")
include("biophysical/shape.jl")

include("core/compensation.jl")
include("core/flux.jl")
include("core/non_stomatal_soilmoisture.jl")
include("core/stomatal_soilmoisture.jl")
include("core/respiration.jl")
include("core/rubisco_regen.jl")

include("interfaces/energy_balance.jl")
include("interfaces/photosynthesis.jl")
include("interfaces/stomatal_conductance.jl")

include("formulations/fvcb.jl")
include("formulations/maespa.jl")
include("formulations/jarvis.jl")
include("formulations/ballberry.jl")
include("formulations/medlyn.jl")
include("formulations/leuning.jl")
include("formulations/emax.jl")
include("formulations/tuzet.jl")

end # module

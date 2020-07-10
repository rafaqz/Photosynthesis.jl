module Photosynthesis

using Unitful,
      FieldDefaults,
      DocStringExtensions,
      Mixers,
      SimpleRoots,
      FieldMetadata

using Unitful: R, °C, K, Pa, kPa, MPa, J, W, kJ, kg, g, m, s, mol, mmol, μmol, σ

import FieldMetadata: @flattenable, @bounds, @default, @description, @units,
                      flattenable, bounds, default, description, units


export enbal!,
       run_enbal!,
       photosynthesis!,
       stomatal_conductance!,
       soil_water_conductance!,
       extremes!,
       co2_compensation_point,
       rubisco_compensation_point,
       rubisco_regeneration,
       rubisco_limited_rate,
       transport_limited_rate,
       slope,
       flux,
       max_electron_transport_rate,
       max_rubisco_activity,
       decoupling, leaftemp,
       gsdiva,
       cmolar,
       ball_berry!,
       jarvis!,
       shape_gs,
       psil_gs!,
       respiration,
       model_init!,
       model_update!,
       converge_tleaf!,
       vjmax_water,
       evapotranspiration,
       penman_monteith,
       factor_conductance,
       radiation_conductance,
       boundary_conductance_free,
       boundary_conductance_forced,
       grashof_number,
       latent_heat_water_vapour,
       saturated_vapour_pressure,
       vapour_pressure_deficit,
       yingping_radiation_conductance,
       forced_boundary_conductance,
       free_boundary_conductance,
       grashof_number,
       arrhenius,
       penman_monteith

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

export AbstractPotentialDependence, LinearPotentialDependence, ZhouPotentialDependence, NoPotentialDependence

export AbstractSoilMethod, NoSoilMethod, VolumetricSoilMethod, ConstantSoilMethod,
       DeficitSoilMethod, PotentialSoilMethod, EmaxSoilMethod, TuzetSoilMethod


export AbstractJarvisCO2, JarvisNoCO2, JarvisLinearCO2, JarvisNonlinearCO2

export AbstractJarvisVPD, JarvisHyperbolicVPD, JarvisLohammerVPD,
       JarvisFractionDeficitVPD, JarvisLinearDeclineVPD

export AbstractJarvisLight, JarvisLight

export AbstractJarvisTemp, JarvisNoTemp, JarvisTemp1, JarvisTemp2


export AbstractStomatalConductanceSubModel, BallBerryStomatalConductanceSubModel, LeuningStomatalConductanceSubModel,
       MedlynStomatalConductanceSubModel, ThreeParStomatalConductanceSubModel, TuzetStomatalConductanceSubModel

export AbstractStomatalConductance, AbstractBallBerryStomatalConductance, BallBerryStomatalConductance,
       TuzetStomatalConductance, AbstractEmaxStomatalConductance, EmaxStomatalConductance, JarvisStomatalConductance

export AbstractPhotosynthesis, AbstractFvCBPhotosynthesis, FvCBPhotosynthesis

export AbstractyFvCBEnergyBalance, FvCBEnergyBalance, EmaxEnergyBalance

export BallBerryVars, EmaxVars, TuzetVars, JarvisVars

@template TYPES =
    """
    $(TYPEDEF)
    $(FIELDS)
    """

@template (FUNCTIONS, METHODS, MACROS) =
    """
    $(SIGNATURES)
    """

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
include("formulations/jarvis/jarvis.jl")
include("formulations/ballberry/ballberry.jl")
include("formulations/ballberry/medlyn.jl")
include("formulations/ballberry/leuning.jl")
include("formulations/ballberry/maespa/maespa.jl")
include("formulations/ballberry/maespa/emax.jl")
include("formulations/ballberry/maespa/tuzet.jl")

end # module

# Photosynthesis

A Farquhar-von-Cammerer-Berry (FvCB) photosynthesis model.

This module at its core is a rewrite of the Maespa/Maestra models in Julia. It
has been written to have a modular structure that can be easily modified or added
to at any level. It is also designed to itself be a component of larger
model.

It aims to provide a comprehensive yet minimalist set of photosynthesis and leaf
energy balance models, eventually including as many as possible published
photosynthesis formulations: but nothing else. Growth models, 3D canopy models,
nutrient uptake etc. should be supplied from another package. The code for these
from Masepas/Maestra has *intentionally* not been included to improve modularity
and reduce complexity.

Only parameters that are actually used must be specified, and all
parameters have sensible defaults.

Photosynthesis.jl uses a hierarchical structure of parameters that are composed
of modular pieces of functionality. They contain the variables the formulation
requires and also are dispatched to the methods specified for the formulation.

All sub-model choices are specified by selecting sub-model types, not `if`
blocks as is common in older models. This means that the set of used parameters
and functions is always known to the compiler, which may allow more automation
of sensitivity analysis and parameter estimation.

# Main photosynthesis routines

These are the main functions that run overall photosynthesis processes.

```@docs
photosynthesis!
```

```@docs
phototranspiration!
```

```@docs
run_phototrans!
```

# Specific functions

Any of these functions can be overridden for a specific parameter types.

```@docs
stomatal_conductance!
```

```@docs
calc_km
```

```@docs
rubisco_compensation
```

```@docs
calc_slope
```

```@docs
calc_jmax
```

```@docs
calc_vcmax
```

```@docs
gsdiva
```

```@docs
cmolar
```

```@docs
respiration
```

```@docs
evapotranspiration
```

```@docs
penman_monteith
```

```@docs
soil_induced_conductance
```

```@docs
factor_conductance
```

```@docs
radiation_conductance
```

```@docs
boundary_conductance_free
```

```@docs
boundary_conductance_forced
```

```@docs
co2_compensation
```

```@docs
rubusico_compenstion
```

```@docs
saturated_vapour_pressure
```


# Type Hierarchy

## Main parameters

```@docs
PhotoParams
```
```@docs
PhotoVars
```

## Model selection

```@docs
AbstractPhotoModel
```
```@docs
AbstractBallBerryModel
```
```@docs
AbstractMaespaModel
```
```@docs
BallBerryModel
```
```@docs
TuzetModel
```
```@docs
EmaxModel
```
```@docs
JarvisModel
```

## CO2Compensation

```@docs
AbstractCO2Compensation
```
```@docs
BadgerCollatzCO2Compensation
```
```@docs
BernacchiCO2Compensation
```
```@docs
AbstractJarvisCO2
```
```@docs
JarvisNoCO2
```
```@docs
JarvisLinearCO2
```
```@docs
JarvisNonlinearCO2
```
```@docs
AbstractJarvisVPD
```
```@docs
JarvisHyperbolicVPD
```
```@docs
JarvisLohammerVPD
```
```@docs
JarvisFractionDeficitVPD
```
```@docs
JarvisLinearDeclineVPD
```

## Stomatal conductance

```@docs
AbstractStomatalConductance
```
```@docs
BallBerryStomatalConductance
```
```@docs
LeuningStomatalConductance
```
```@docs
MedlynStomatalConductance
```
```@docs
ThreeParStomatalConductance
```
```@docs
TuzetStomatalConductance
```

# Photosynthesis

```@docs
AbstractJmax
```
```@docs
AbstractVcmax
```
```@docs
AbstractVcJmax
```
```@docs
Jmax
```
```@docs
NoOptimumVcmax
```
```@docs
OptimumVcmax
```
```@docs
VcJmax
```
```@docs
DukeVcJmax
```
```@docs
AbstractRubiscoRegen
```
```@docs
RubiscoRegen
```
```@docs
AbstractRespiration
```
```@docs
Respiration
```
```@docs
AbstractRadiationConductance
```
```@docs
YingPingRadiationConductance
```
```@docs
AbstractBoundaryConductance
```
```@docs
BoundaryConductance
```

## Other

```@docs
AbstractDecoupling
```
```@docs
McNaughtonJarvisDecoupling
```
```@docs
NoDecoupling
```

## Soil water 

```@docs
AbstractSoilMethod
```
```@docs
VolumetricSoilMethod
```
```@docs
ConstantSoilMethod
```
```@docs
DeficitSoilMethod
```
```@docs
PotentialSoilMethod
```
```@docs
EmaxSoilMethod
```
```@docs
TuzetSoilMethod
```

### Data

Kinds of soil water data available.

```@docs
AbstractSoilData
```
```@docs
AbstractDeficitSoilData
```
```@docs
AbstractContentSoilData
```
```@docs
DeficitSoilData
```
```@docs
ContentSoilData
```
```@docs
SimulatedSoilData
```
```@docs
PotentialSoilData
```
```@docs
NoSoilData
```

### Dependence 

Dependence of vcmax and jmax on soil moisture. This is used in the Emax and
Tuzet models with the [`vjmax_water`](@ref) function.

```@docs
AbstractPotentialDependence
```
```@docs
LinearPotentialDependence
```
```@docs
ZhouPotentialDependence
```
```@docs
NoPotentialDependence
```

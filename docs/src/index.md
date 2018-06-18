
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

## Main photosynthesis routines

These are the main functions that run overall photosynthesis processes.

```@docs
photosynthesis!
phototranspiration!
run_phototrans!
```

## Specific functions

Any of these functions can be overridden for a specific parameter types.

```@docs
stomatal_conductance!
calc_km
rubisco_compensation
calc_slope
calc_jmax
calc_vcmax
gsdiva
cmolar
respiration
evapotranspiration
penman_monteith
soil_induced_conductance
factor_conductance
radiation_conductance
boundary_conductance_free
boundary_conductance_forced
co2_compensation
rubusico_compenstion
saturated_vapour_pressure
```


## Type Hierarchy

### Main parameters

```@docs
PhotoParams
PhotoVars
```

### Model selection

```@docs
FvCBPhoto
AbstractPhotoModel
AbstractBallBerryModel
AbstractMaespaModel
BallBerryModel
TuzetModel
EmaxModel
JarvisModel
```

### CO2 and Rubisco Compensation

```@docs
AbstractCompensation
BadgerCollatzCompensation
BernacchiCompensation
```

### Stomatal conductance

```@docs
AbstractStomatalConductance
BallBerryStomatalConductance
LeuningStomatalConductance
MedlynStomatalConductance
ThreeParStomatalConductance
TuzetStomatalConductance
```

### Photosynthesis

```@docs
AbstractJmax
AbstractVcmax
AbstractVcJmax
Jmax
NoOptimumVcmax
OptimumVcmax
VcJmax
DukeVcJmax
AbstractRubiscoRegen
RubiscoRegen
AbstractRespiration
Respiration
AbstractRadiationConductance
YingPingRadiationConductance
AbstractBoundaryConductance
BoundaryConductance
```

### Decoupling

```@docs
AbstractDecoupling
McNaughtonJarvisDecoupling
NoDecoupling
```

### Soil water 

```@docs
AbstractSoilMethod
VolumetricSoilMethod
ConstantSoilMethod
DeficitSoilMethod
PotentialSoilMethod
EmaxSoilMethod
TuzetSoilMethod
```

### Soil Data

Kinds of soil water data available.

```@docs
AbstractSoilData
AbstractDeficitSoilData
AbstractContentSoilData
DeficitSoilData
ContentSoilData
SimulatedSoilData
PotentialSoilData
NoSoilData
```

### Dependence 

Dependence of vcmax and jmax on soil moisture. This is used in the Emax and
Tuzet models with the [`vjmax_water`](@ref) function.

```@docs
AbstractPotentialDependence
LinearPotentialDependence
ZhouPotentialDependence
NoPotentialDependence
```

### Jarvis Model

```@docs
AbstractJarvisCO2
JarvisNoCO2
JarvisLinearCO2
JarvisNonlinearCO2
AbstractJarvisVPD
JarvisHyperbolicVPD
JarvisLohammerVPD
JarvisFractionDeficitVPD
JarvisLinearDeclineVPD
```


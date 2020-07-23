
```@docs
Photosynthesis
```


## Energy balance model

```@docs
Photosynthesis.AbstractEnergyBalance
Photosynthesis.AbstractLeuningEnergyBalance
LeuningEnergyBalance
enbal!
Photosynthesis.enbal_init!
Photosynthesis.enbal_update!
```

### Model variables

```@docs
BallBerryVars
```

## Biophysical components

```@docs
Photosynthesis.saturated_vapour_pressure
Photosynthesis.vapour_pressure_deficit
Photosynthesis.latent_heat_water_vapour
Photosynthesis.penman_monteith
Photosynthesis.arrhenius
Photosynthesis.leaftemp
```

### Boundary conductance

```@docs
AbstractBoundaryConductance
BoundaryConductance
Photosynthesis.boundary_conductance_free
Photosynthesis.boundary_conductance_forced
Photosynthesis.vapour_conductance! 
Photosynthesis.grashof_number
Photosynthesis.cmolar
Photosynthesis.free_boundary_conductance
Photosynthesis.forced_boundary_conductance
```

### Canopy-atmosphere decoupling

```@docs
AbstractDecoupling
McNaughtonJarvisDecoupling
NoDecoupling
Photosynthesis.decoupling
```

### Evapotranspiration

```@docs
AbstractEvapotranspiration 
PenmanMonteithEvapotranspiration
Photosynthesis.evapotranspiration
Photosynthesis.penman_monteith(ρa, Δ, lhv, Rn, Da, gh, gv)
Photosynthesis.slope
```

### Radiation conductance

```@docs
AbstractRadiationConductance
WangRadiationConductance
radiation_conductance
wang_radiation_conductance
```

## Photosynthesis model

```@docs
AbstractPhotosynthesis
AbstractFvCBPhotosynthesis
FvCBPhotosynthesis
photosynthesis!
```

### CO2 and Rubisco Compensation

```@docs
Compensation
BadgerCollatzCompensation
BernacchiCompensation
Photosynthesis.co2_compensation_point
Photosynthesis.rubisco_compensation_point
```

### Flux

```@docs
AbstractFlux
Flux
DukeFlux
PotentialModifiedFlux
Photosynthesis.flux
```

#### Electron transport rate

```@docs
AbstractJmax
Jmax
Photosynthesis.max_electron_transport_rate
```

### Rubisco activity 

```@docs
AbstractVcmax
OptimumVcmax
NoOptimumVcmax
Photosynthesis.max_rubisco_activity
```

### Rubisco regeneration 

```@docs
AbstractRubiscoRegen
RubiscoRegen
Photosynthesis.rubisco_regeneration
```

### Respiration

```@docs
AbstractRespiration
Respiration
AcclimatizedRespiration
Photosynthesis.respiration
```

### Non-stomatal soil water-potential dependence 

```@docs
AbstractPotentialDependence
LinearPotentialDependence
ZhouPotentialDependence
NoPotentialDependence
Photosynthesis.non_stomatal_potential_dependence
```

## Stomatal conductance models

```@docs
AbstractStomatalConductance
AbstractBallBerryStomatalConductance
BallBerryStomatalConductance
Photosynthesis.stomatal_conductance!
Photosynthesis.rubisco_limited_rate
Photosynthesis.transport_limited_rate
Photosynthesis.gs_init!
Photosynthesis.gs_update!
Photosynthesis.check_extremes!
Photosynthesis.update_extremes!
```

### Stomatal conductance sub-models

```@docs
AbstractStomatalConductanceSubModel
AbstractBallBerryStomatalConductanceSubModel
BallBerryStomatalConductanceSubModel
LeuningStomatalConductanceSubModel
MedlynStomatalConductanceSubModel
Photosynthesis.gs_div_a
```

### Stomatal conductance shape 

```@docs
StomatalConductanceShape 
HardMinimum 
HyperbolicMinimum
Photosynthesis.shape_gs
```


### Stomatal soil water dependence

```@docs
AbstractSoilMethod
VolumetricSoilMethod
ConstantSoilMethod
DeficitSoilMethod
PotentialSoilMethod
NoSoilMethod
Photosynthesis.soilmoisture_conductance
```

#### Soil Data

Kinds of soil water data used in [`DeficitSoilMethod`](@ref) and
[`VolumetricSoilMethod`](@ref).

```@docs
AbstractSoilData
AbstractDeficitSoilData
AbstractContentSoilData
DeficitSoilData
ContentSoilData
SimulatedSoilData
```

## EMAX Model

```@docs
EmaxEnergyBalance
EmaxStomatalConductance
EmaxSoilMethod
EmaxVars
```

## Tuzet Model

```@docs
TuzetEnergyBalance
TuzetStomatalConductance
TuzetStomatalConductanceSubModel
TuzetVars
Photosynthesis.leaf_water_potential_finder
Photosynthesis.fpsil
```

## Jarvis Model

```@docs
JarvisStomatalConductance
Photosynthesis.factor_conductance
AbstractJarvisLight
JarvisLight
Photosynthesis.light_factor

AbstractJarvisTemp
JarvisNoTemp
JarvisTemp1
JarvisTemp2
Photosynthesis.temp_factor

AbstractJarvisCO2
JarvisNoCO2
JarvisLinearCO2
JarvisNonlinearCO2
Photosynthesis.co2_factor

AbstractJarvisVPD
JarvisHyperbolicVPD
JarvisLohammerVPD
JarvisFractionDeficitVPD
JarvisLinearDeclineVPD
Photosynthesis.vpd_factor
```

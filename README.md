# Photosynthesis

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://rafaqz.github.io/Photosynthesis.jl/dev)
[![Build Status](https://travis-ci.org/rafaqz/Photosynthesis.jl.svg?branch=master)](https://travis-ci.org/rafaqz/Photosynthesis.jl)
[![codecov.io](http://codecov.io/github/rafaqz/Photosynthesis.jl/coverage.svg?branch=master)](http://codecov.io/github/rafaqz/Photosynthesis.jl?branch=master)

A Farquhar-von-Cammerer-Berry (FvCB) photosynthesis modelling framework.

This module at its core is a rewrite of the Maespa/Maestra photosynthesis models in Julia. 
It has been written to have a modular structure that can be easily modified or
added to at any level. It is also designed to itself be a component of larger model.

It aims to provide a comprehensive yet minimalist set of photosynthesis and leaf
energy balance models, eventually including as many as possible published
photosynthesis formulations. But nothing else. Growth models, 3D canopy models,
nutrient uptake etc. can be supplied from another package.

Included formulations of the basic FvCB model include Ball-Berry (and multiple
sub-variants), EMAX, Tuzet, and Jarvis formulations.

These are written out in interchangeable components: you can write you own
components in a script and swap them into the model - for _any_ part of this
package. You Just need to define the `struct` for parameters, and define the
corresponding method for the formulation you with to write.

This means the old problem of multiple forks with minor formulation
changes is essentially solved. Never edit the source code of this package,
just define new alternate components with your desired changes.


## Example

Here we will define a basic Ball-Berry model, and run it:

```julia
params = FvCBEnergyBalance(
    photosynthesis_model=FvCBPhotosynthesis(
        stomatal_conductance=BallBerryStomatalConductance(
            gs_submodel=BallBerryStomatalConductance(),
            soil_model=NoSoilMethod()
        ),
        flux_model=DukeFlux(),
        compensation_model=BadgerCollatzCompensation()
    )
)
vars = BallBerryVars()

enbal!(vars, params)
```

Notice we run `enbal!` not `photosynthesis!`, as `enbal!` runs
`photosynthesis!` in a loop to calculate both the temperature
and the assimilated C.


# Notes about strategies used in this package

This package was intentionally written as an exploration of how process model
components could generally be structured in better ways to traditional models -
to reduce the overheads of making small changes to any part of the model and
improve the potential for collaboration. This method also facilitates composing
formulations from multiple packages into new models without having to change the
codebases.

##  Units

This implementation takes advantage of
[Unitful.jl](https://github.com/PainterQubits/Unitful.jl) to provide unit
checking and conversion. This removes the possibility of a large class of
errors, especially when adding your own custom formulations, and has little or
no runtime cost in most cases.

## Nested parameters

Parameters are composed of nested `struct`s, so that components can be modified
and new components can be added without altering the original. This structure
should be liberating for those wanting to make modifications to existing
formulations, and to share them.

Formulations all follow the pattern of defining an abstract type and a
function at the top of the file. All formulation variants inherit from the
abstract type, and define a method of the function.

This allows easy addition of any model components without altering the package
code.

To convert model parameters to a vector for optimisation or interactive
interfaces, [Flatten.jl](https://github.com/rafaqz/Flatten.jl) can be used.

[FieldMetadata.jl](https://github.com/rafaqz/FieldMetadata.jl) allows adding
other metadata to parameters (all the `| XX |` in struct definitions), such as
descriptions, defaults, units and bounds. These can also be flattened with
Flatten.jl for used in an optimiser or interface without writing any custom
code.

Other tools like the (unpublished)
[Codify.jl](https://github.com/rafaqz/Codify.jl) turns a nested model struct
into the code for keyword argument constructors that would build it. 


# Issues

## Tests

Currently this package has _no real tests_, as there were none for the original
Fortran. Instead it tests against the original formulations. Many of these test
pass. But Tuzet, Jarvis, and EMAX models are as yet untested.

While very few photosynthesis models actually have tests, we should work towards
Photosynthesis.jl being properly tested so that we know that formulations
do what we intend them to.

## Non-functional variables struct 

Temp variables are written directly to a variables struct. This could be
replaced with a functional style where variables are rebuilt and explicitly
updated. This would mean that the package can be used on GPUs, and so the
updating of state flow is easier to follow.

# Credits

Thanks got to the creators of the Maestra Fortran libraries including Belinda
Medlyn and Remko Duursma, and all earlier contributors to Maestra and Maestro.

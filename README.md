# Photosynthesis

[![Build Status](https://travis-ci.org/rafaqz/Photosynthesis.jl.svg?branch=master)](https://travis-ci.org/rafaqz/Photosynthesis.jl)
[![codecov.io](http://codecov.io/github/rafaqz/Photosynthesis.jl/coverage.svg?branch=master)](http://codecov.io/github/rafaqz/Photosynthesis.jl?branch=master)

A Farquhar-von-Cammerer-Berry (FvCB) photosynthesis modelling framework.

This module at its core is a rewrite of the Maespa/Maestra/ models in Julia. It
has been written to have a modular structure that can be easily modified or
added to at any level. It is also designed to itself be a component of larger
model.

It aims to provide a comprehensive yet minimalist set of photosynthesis and leaf
energy balance models, eventually including as many as possible published
photosynthesis formulations: but nothing else. Growth models, 3D canopy models,
nutrient uptake etc. should be supplied from another package.

Included formulations of the basic FvCB model include Ball-Berry (and multiple
sub-variants), EMAX, Tuzet, and Jarvis.

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
codebases. Some of the structural ideas seem promising, others may need
rethinking.

##  Units

This implementation takes advantage of
[Unitful.jl](https://github.com/PainterQubits/Unitful.jl) to provide unit
checking and conversion. This removes the possibility of a large class of
errors, especially when adding your own custom formulations, and has little or
no runtime cost in most cases

## Nested parameters

Parameters are composed of nested `struct`s, so that components can be modified
and new components can be added without altering the original. This structure
should be liberating for those wanting to make modifications to existing
formulations.

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
into the code for kwarg constructors that would build it. 


# Issues

## Tests

Currently this package has _no real tests_, as there were none for the original
Fortran. Instead it tests against the original formulations. Many of these test
pass, but some still don't.

Obviously this situation needs to change. Unit tests and comprehensive
integration tests are required. But is a significant amount of work
and it would need a more specific application and larger user base to be
justifiable.

## Non-functional variable struct 

Temp variables are written directly to a variables struct. This should be
replaced with a functional style where variables are rebuilt and explicitly
updated. This would mean that the package can be used on GPUs, and so the
updating of state flow is easier to follow.

# Credits

Thanks got to the creators of the Maestra Fortran libraries including Belinda
Medlyn and Remko Duursma, and all earlier contributors to Maestra and Maestro.

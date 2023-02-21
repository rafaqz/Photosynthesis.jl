# Photosynthesis


## *This package is not registered and is not currently maintained*

I don't work on photosynthesis currently and have no personal need to use this. I also currently work on a lot of other packages and have very little time to spare to update this. 

There is a lot of useful code here, if you need it you will need to do some work for it a little. Likely a fraction of what I put into writing it. It's fast, well structured and mostly verified to be correct. But its not a polished experience - this was a proof of concept from a masters project. 

If you want examples of everything that currently works, *Read the tests in the test folder*. You can work through the functions in there and see how it all works.

Just ignore the code required to run the original fortran binaries I tested againsts.


[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://rafaqz.github.io/Photosynthesis.jl/dev)
[![Build Status](https://travis-ci.com/rafaqz/Photosynthesis.jl.svg?branch=master)](https://travis-ci.com/rafaqz/Photosynthesis.jl)
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
v = BallBerryVars()
p1 = MaespaEnergyBalance(
    photosynthesis=FvCBPhotosynthesis(
        stomatal_conductance=EmaxStomatalConductance(
            gs_submodel=BallBerryStomatalConductanceSubModel(),
            soilmethod=EmaxSoilMethod(),
        ),
        flux=DukeFlux(),
        compensation=BadgerCollatzCompensation(),
        respiration_model=Respiration(),
    ),
    atol=0.005K,
)
enbal!(v, p)
```


Formulations using the MAESPA energy balance and FvCB photosynthesis, 
all tested against MAESPA:

- BallBerry stomatal conducance types:
  - `BallBerry`
  - `Medlyn`
  - `Leuning`
- `Emax`


Other formulations, working but not tested
- `Tuzet` : PSILFIND is not yet tested against MAESPA. Use with caution.
- `Jarvis` : Not tested against Maestra, also use with caution.



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


# Tests

Currently this package has extensive tests against the original formaultions in the 
MAESPA package. The diffrences between 32bit and 64bit usually breaks parity in anwsers
when some value get to the extreme end of ranges, such as above 50 degrees C.

The Jarvis model is not tested, as it is not included in Maespa. The outer loops of the 
Tuzet model are also not yet tested. Please put in an issue if you would like to use
the Tuzet model.

There are, however _no real unit tests_, as there were none for the original Fortran. 
It would, of course be good to write them, but beyond the current scope of writing 
this package.


# Credits

Thanks got to the creators of the Maestra Fortran libraries including Belinda
Medlyn and Remko Duursma, and all earlier contributors to Maestra and Maestro.

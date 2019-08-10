# Photosynthesis

[![Build Status](https://travis-ci.org/rafaqz/Photosynthesis.jl.svg?branch=master)](https://travis-ci.org/rafaqz/Photosynthesis.jl)
[![Coverage Status](https://coveralls.io/repos/rafaqz/Photosynthesis.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/rafaqz/Photosynthesis.jl?branch=master)
[![codecov.io](http://codecov.io/github/rafaqz/Photosynthesis.jl/coverage.svg?branch=master)](http://codecov.io/github/rafaqz/Photosynthesis.jl?branch=master)

This package is a rewrite of core photosynthesis model from the Maespa and
Maestra Fortran libraries of Belinda Medlyn and Remko Duursma, as well as some
direction from Duursmas' R package YplantQMC. Many thanks go to them and early
authors of the original Maestro package.

The model formulation has not been changed. But it is just heavily refactored to
be as flexible and modular as possible. This version should give similar
performance to the original Fortran, but from significantly less lines of code.

It's possible to use custom methods for most aspects of photosynthesis without
changing the code this package, simply by writing your own version that dispatches 
on your custom parameter type.

There is also complete parameter modularity: you only have to set the parameters
that are actually used by the specific formulations you choose. On top of this,
all parameters come with default values that can be overridden by passing
keyword arguments to the constructor, shuch as:

```
params = PhotoParams(model=TuzetModel())
```

In future this package may be extended to include C4 and CAM photosynthesis and
alternate stomatal conductance models. But other kinds of functionality should be
added in external packages, to keep the focus specifically on the process of
photosynthesis. Efforts will be made to make this as widely interoperable as
possible.


# Notes about strategies used in this package

This package was intentionally written as an exploration of how process model
components could generally be structured differently to traditional models -
primarily to reduce the overheads of making small changes to any part of the
model. This method also facilitates composing formulations from multiple packages into new
models without having to change the codebases. Some of the structural
ideas seem promising, others may need rethinking.


##  Units
This implementation takes advantage of
[Unitful.jl](https://github.com/PainterQubits/Unitful.jl) to provide 
unit checking and conversion. This removes the possibility of a large
class of errors, especially when adding your own custom formulations, and 
has little or no runtime cost in most cases

## Nested parameters
Parameters are composed of nested structs so that components can be modified
and new components can be added using method dispatch without altering the 
original. This structure should be liberating for those wanting to make
modifications to existing formulations.

However, this means parameters aren't readily available as a vector for use in 
an optimiser. [Flatten.jl](https://github.com/rafaqz/Flatten.jl) mostly solves 
this problem, and flattens any nested struct to a tuple very efficiently. 
This also makes them available for building Interact interfaces.

[FieldMetadata.jl](https://github.com/rafaqz/FieldMetadata.jl) allows adding
other metadata to parameters (all the `| x |` in struct definitions), such as 
descriptions, defaults, units and bounds. These can also be flattened with 
Flatten.jl for used in an optimiser or interface without writing any custom code.

Other tools like the (unpublished) [Codify.jl](https://github.com/rafaqz/Codify.jl) 
turns a nested model struct into the kwarg code that would build it. Its a hack, 
but may be useful until a better solution exists.

## Interface / implementation

Separating out different implementations makes it much easier to see what a
specific implementation does/changes, but harder to follow the flow of the
models.

## Tests
This package has no real tests, as there were none for the original fortran. 
Instead it tests against the original formulations, many of which do actually work, 
some that still don't.

Obviously this situation needs to change, but is a significant amount of work
that needs a more specific application to be justifiable.


# Structural problems

This package was started after six months writing Julia, and while it has seen 
a few iterations, on reflection there are things that could be better

## Direct field access
The main non-julian style and limit to modularity used here is a lot of direct
field access. These should really be restricted to leaf nodes in nested 
parameter struct, for one-off parameter access in methods that dispatch
on the containing struct type.

Otherwise all field access should be replaced with getter methods.

## Non-functional variable struct 
Temp variables are written directly to a variables struct.
This should be replaced with a functional style where variables are rebuilt and
explicitly updated: so that the package can be used on GPUs, and so the updating
of state flow is easier to follow.

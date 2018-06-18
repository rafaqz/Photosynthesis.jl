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

This implementation also takes advantage of Unitful.jl to provide comprehensive
unit checking and conversion *with basically no runtime cost*. All inputs and
outputs must have the correct units. This removes the possibility of a large
class of errors, especially when adding your own custom formulations.

It's possible to use custom methods for most aspects of photosynthesis without
changing the code this package, simply by importing the function you need to
change with `import` and writing your own version that dispatches on your custom
parameter type. Pull requests for of any additional models would be greatly
appreciated.

There is also complete parameter modularity: you only have to set the parameters
that are actually used by the specific formulations you choose. On top of this,
all parameters come with default values that can be overridden by passing
keyword arguments to the constructor, shuch as:

```
params = PhotoParams(model=TuzetModel())
```

In future this package will be extended to include C4 and CAM photosynthesis and
alternate stomatal conductance models. But other kinds of functionality will be
added in external packages, to keep the focus specifically on the process of
photosynthesis. Efforts will be made to make this as widely interoperable as
possible.


# COSMOAccelerators.jl

## Description
COSMOAccelerators defines an abstract `AbstractAccelerator` type that can be used as an interface to write accelerator methods for algorithms based on non-expansive operators, e.g. fixed-point iterations. 
The acceleration methods in this package were originally developed for the convex conic solver [COSMO.jl](https://github.com/oxfordcontrol/COSMO.jl), but can also be used with other fixed-point methods, e.g. algorithms based on the [ProximalOperators.jl](https://github.com/kul-forbes/ProximalOperators.jl) package.

## Installation

`COSMOAccelerators` can be added via the Julia package manager (type `]`): `pkg> add COSMOAccelerators`

## Usage
The `AbstractAccelerator` interface requires implementation of the following methods:





## Related package
- [COSMO.jl](https://github.com/oxfordcontrol/COSMO.jl) - an ADMM based solver for convex conic optimisation problems.
- [ProximalOperators.jl](https://github.com/kul-forbes/ProximalOperators.jl) - a set of operator building blocks for proximal algorithms

## Credits
COSMOAccelerators.jl is developed by [Michael Garstka](www.migarstka.com) at University of Oxford.


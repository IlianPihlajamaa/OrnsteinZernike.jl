# Systems

A `System` object specifies what type of system needs to be solved. This includes the dimensionality, density and temperature, as well as the interaction potential as defined according to the previous page. Right now, this package exports only one type of `System`: `SimpleLiquid`. 

## `SimpleLiquid`

The `SimpleLiquid` system assumes that the potentials involved are rotationally symmetric, and therefore depends only on center-to-center distance. Additionally, it assumes that there is no present external field. It can be used both in the single component case, as well as for multi-component systems (mixtures).

For single component systems, the density `ρ` is assumed to be a scalar, while for multicomponent systems, it should be a vector. Internally, the package converts it into a diagonal matrix.

Example 1: a 2-dimensional Lennard-Jones system

```@example 1
using OrnsteinZernike 
potential = SingleComponentLennardJones(1.0, 1.0)
kBT = 1.0
ρ = 0.5
dims = 2
system = SimpleLiquid(dims, ρ, kBT, potential)
```

Example 1: a 3-dimensional 10-component hard-sphere system

```@example 1
using OrnsteinZernike # hide 
D = 0.1:0.1:1.0
potential = HardSpheres(D)
kBT = 1.0
ρ = ones(10)/10
dims = 3
system = SimpleLiquid(dims, ρ, kBT, potential)
```

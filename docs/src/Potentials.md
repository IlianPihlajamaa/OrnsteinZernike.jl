# Interaction Potentials

This package defines several potentials that can be used out of the box. It is straightforward to implement your own potential, see [Defining your own potentials](@ref). To evaluate the potential, call `OrnsteinZernike.evaluate_potential(potential, r)`. For example:

```@example
using OrnsteinZernike, Plots
r = 0.9:0.01:4.0
potential = SingleComponentLennardJones(1.0, 1.0)
u = OrnsteinZernike.evaluate_potential(potential, r)
plot(r, u, xlabel="r", ylabel="u(r)", ylims=(-1,1), label=nothing)
```

## Implemented interaction potentials
Below is a list of implemented closures. We use the notation shown in the [Theory](@ref) section.

```@autodocs
Modules = [OrnsteinZernike]
Filter = t -> (typeof(t) === DataType || typeof(t) === UnionAll) && t <: OrnsteinZernike.Potential
Private = false
```
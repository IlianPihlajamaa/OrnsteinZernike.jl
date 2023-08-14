# Closures

This package defines several closure relations that can be used out of the box. It is straightforward to implement your own closure, see [Defining your own closure](@ref). 

```@example
using OrnsteinZernike, Plots
closure = PercusYevick()
```

## Implemented Closures
Below is a list of implemented closures. We use the notation shown in the [Theory](@ref) section.

```@autodocs
Modules = [OrnsteinZernike]
Filter = t -> ((typeof(t) === DataType || typeof(t) === UnionAll) && t <: OrnsteinZernike.Closure)
Private = false
```
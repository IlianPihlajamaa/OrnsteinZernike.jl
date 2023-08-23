# Closures

This package defines several closure relations that can be used out of the box. It is straightforward to implement your own closure, see [Defining your own closure](@ref). 

## Implemented Closures
Below is an alphabetical list of implemented closures. We use the notation shown in the [Theory](@ref) section.

```@index
Modules = [OrnsteinZernike]
Pages = ["Closures.md"]
Order   = [:type]
Private = false
```

```@autodocs
Modules = [OrnsteinZernike]
Filter = t -> ((typeof(t) === DataType || typeof(t) === UnionAll) && t <: OrnsteinZernike.Closure)
Private = false
```
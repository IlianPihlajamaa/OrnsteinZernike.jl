# Closures

This package defines several closure relations that can be used out of the box. It is straightforward to implement your own closure, see [Defining your own closure](@ref). 

Some closures, for example the Rogers-Young closure, include free parameters that may be fixed by the requirement of thermodynamic consistency. See the [Thermodynamic Consistency](@ref) page for an example. 

Some closures use a renormalized indirect correlation function $\gamma^*(r) = \gamma(r) - u_{LR}(r)$ instead of the standard one. Here, $u_{LR}(r)$ is the long range tail of the interaction potential. There are several ways in which the interaction potential can be split into a short-range and long range part, the most common one is the Weeks-Chandler-Andersen construction. In order to use these ... 

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
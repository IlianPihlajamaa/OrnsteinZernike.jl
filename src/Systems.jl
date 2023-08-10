
"""
    System

Abstract type for holding information about the system that needs to be solved
"""
abstract type System end


"""
    SimpleLiquid{dims, ...} <: System

Holds information about a homogeneous, isotropic system with radially symmetric pair interactions. `dims` is the dimensionality.

Construct using

`SimpleLiquid(dims, ρ, kBT, potential)`

Fields:

- ρ: number density, must be either a `Number` in case of a single component system, or a `Vector` in case of a mixture. In the latter case, each element contains the number density of the respective component.
- kBT: thermal energy
- potential::Potential: the interaction potential.  

Examples:

```julia
ρ = 0.5; kBT = 1.1; dims = 3
pot = SingleComponentHardSpheres()
system = SimpleLiquid(dims, ρ, kBT, pot)
```

```julia
ρ = [0.5, 0.1]; kBT = 5.2; dims = 3
pot = MultiComponentHardSpheres([1.0, 0.8])
system = SimpleLiquid(dims, ρ, kBT, pot)
```
"""
mutable struct SimpleLiquid{dims, species, T1, T2, P} <: System
    ρ::T1
    kBT::T2
    potential::P
end

function SimpleLiquid(dims, ρ::Number, kBT, potential)
    Tρ = typeof(ρ)
    TkT = typeof(kBT)
    Tpot = typeof(potential)
    return SimpleLiquid{dims, 1, Tρ, TkT, Tpot}(ρ, kBT, potential)
end

function SimpleLiquid(dims, ρ::AbstractVector, kBT, potential)
    TkT = typeof(kBT)
    Tpot = typeof(potential)
    Ns = length(ρ)
    ρ = Diagonal(SVector{Ns}(ρ))
    Tρ = typeof(ρ)
    return SimpleLiquid{dims, Ns, Tρ, TkT, Tpot}(ρ, kBT, potential)
end

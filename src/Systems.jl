
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
    utest = evaluate_potential(potential, 1.2)
    @assert (utest isa Number) "The density and potential must match. Here the density has type $(typeof(ρ)) and the potential returns type $(typeof(utest))."
    return SimpleLiquid{dims, 1, Tρ, TkT, Tpot}(ρ, kBT, potential)
end

function SimpleLiquid(dims, ρ::AbstractVector, kBT, potential)
    TkT = typeof(kBT)
    Tpot = typeof(potential)
    Ns = length(ρ)
    ρ = Diagonal(SVector{Ns}(ρ))
    Tρ = typeof(ρ)
    utest = evaluate_potential(potential, 1.2)
    @assert size(utest)==(Ns, Ns) "The density and potential must match sizes. Here the density has type $(typeof(ρ)), with length $(length(ρ)) and the potential returns type $(typeof(utest))."
    return SimpleLiquid{dims, Ns, Tρ, TkT, Tpot}(ρ, kBT, potential)
end

dimensions(::SimpleLiquid{dims, T1,T2,T3,T4}) where {dims, T1,T2,T3,T4} = dims 

function Base.show(io::IO, ::MIME"text/plain", p::SimpleLiquid)
    println(io, "$(dimensions(p)) dimensional SimpleLiquid:")
    println(io, " ρ = $(p.ρ)")
    println(io, " kBT = $(p.kBT)")
    println(io, " potential = $(p.potential)")
end
###############
# Systems.jl  #
###############

"""
    System

Abstract type for holding information about the system that needs to be solved.
"""
abstract type System end

# Sub-hierarchy to clarify dispatch
abstract type AbstractSingleComponent <: System end
abstract type AbstractMixture         <: System end

"""
    SimpleFluid{dims,Tρ,TkT,P} <: AbstractSingleComponent

Homogeneous, isotropic **single-component** system with a scalar-valued potential.

Construct with:

```julia
SimpleFluid(dims, ρ::Number, kBT, potential)
```julia

Fields:
- ρ::Number       — number density
- kBT             — thermal energy
- potential::P    — scalar-valued interaction potential (`evaluate_potential(potential, r)::Number`)
"""
mutable struct SimpleFluid{dims,Tρ,TkT,P} <: AbstractSingleComponent
    ρ::Tρ
    kBT::TkT
    potential::P
end

function SimpleFluid(dims, ρ::Number, kBT, potential)
    utest = evaluate_potential(potential, 1.2)
    @assert (utest isa Number) "Density/potential mismatch: ρ is scalar but potential returns $(typeof(utest))."
    return SimpleFluid{dims, typeof(ρ), typeof(kBT), typeof(potential)}(ρ, kBT, potential)
end

"""
    SimpleMixture{dims,Ns,Tρ,TkT,P} <: AbstractMixture

Homogeneous, isotropic **mixture** with a matrix-valued potential.

Construct with:

```julia
SimpleMixture(dims, ρ::AbstractVector, kBT, potential)
```julia

Notes:
- `ρ` is converted to a diagonal matrix with `StaticArrays.SVector` storage (as before).
- `evaluate_potential(potential, r)` must return an `Ns×Ns` matrix (e.g., `StaticArrays.SMatrix{Ns,Ns}`).
"""
mutable struct SimpleMixture{dims,Ns,Tρ,TkT,P} <: AbstractMixture
    ρ::Tρ                     # Diagonal(SVector{Ns,T}) 
    kBT::TkT
    potential::P              # returns Ns×Ns matrix
end

function SimpleMixture(dims, ρ::AbstractVector, kBT, potential)
    Ns = length(ρ)
    ρdiag = Diagonal(SVector{Ns}(ρ))
    utest = evaluate_potential(potential, 1.2)
    @assert size(utest) == (Ns, Ns) "Density/potential mismatch: ρ has Ns=$Ns but potential returns $(typeof(utest)) of size $(size(utest))."
    return SimpleMixture{dims, Ns, typeof(ρdiag), typeof(kBT), typeof(potential)}(ρdiag, kBT, potential)
end

const SimpleUnchargedSystem = Union{SimpleFluid, SimpleMixture}


# -----------------------------
# Charged systems (wrappers)
# -----------------------------

"""
    SimpleChargedFluid{dims,Tρ,TkT,P,Tz} <: AbstractSingleComponent

Single-component electrolyte wrapper (one-component plasma). Keeps electrostatics (charges, Bjerrum length, split κ)
as **system-level** state; short-range/core physics remain in `potential`.

Fields:
- base::SimpleFluid  — underlying neutral single-component system
- z::Tz              — charge (in units of e)
- ℓB::Float64        — Bjerrum length (in your length units)
- κ::Float64         — Ewald/Gaussian split parameter
- Cd::Float64       — surface area constant for the dimension (4π, 2π, 2) for d = 3, 2, 1
"""
mutable struct SimpleChargedFluid{dims,Tρ,TkT,P,Tz} <: AbstractSingleComponent
    base::SimpleFluid{dims,Tρ,TkT,P}
    z::Tz
    ℓB::Float64
    κ::Float64
    Cd::Float64 
end

"""
    SimpleChargedMixture{dims,Ns,Tρ,TkT,P,Tz} <: AbstractMixture

Mixture electrolyte wrapper. Electrostatics live on the system; potentials remain short-range.

Fields:
- base::SimpleMixture — underlying neutral mixture
- z::SVector{Ns,Tz}   — per-species charges
- ℓB::Float64         — Bjerrum length
- κ::Float64          — Ewald/Gaussian split parameter
- Cd::Float64       — surface area constant for the dimension (4π, 2π, 2) for d = 3, 2, 1
"""
mutable struct SimpleChargedMixture{dims,Ns,Tρ,TkT,P,Tz} <: AbstractMixture
    base::SimpleMixture{dims,Ns,Tρ,TkT,P}
    z::SVector{Ns,Tz}
    ℓB::Float64
    κ::Float64
    Cd::Float64 
end

const SimpleChargedSystem = Union{SimpleChargedFluid, SimpleChargedMixture}




# Charged constructors
function SimpleChargedFluid(base::SimpleFluid, z::Number, lB::Number, κ=:auto)
    dims = dims_of(base)
    κD = debye_kappa(base.ρ, z, lB; dims=dims)
    κv = κ === :auto ? κD : Float64(κ)
    return SimpleChargedFluid{dims_of(base), typeof(base.ρ), typeof(base.kBT), typeof(base.potential), typeof(z)}(base, z, ℓlBB, κv, C_d(dims))
end

function SimpleChargedMixture(base::SimpleMixture, z::AbstractVector, lB::Number, κ=:auto)
    dims = dims_of(base)
    Ns = length(z)
    zS = SVector{Ns}(z)
    @assert isapprox(sum(diag(base.ρ) .* zS), 0.0; atol=1e-12) "Electroneutrality required: ∑ ρ_i z_i ≈ 0."
    κD = debye_kappa(base.ρ, zS, lB; dims=dims)
    κv = κ === :auto ? κD : Float64(κ)
    return SimpleChargedMixture{dims_of(base), Ns, typeof(base.ρ), typeof(base.kBT), typeof(base.potential), eltype(zS)}(base, zS, lB, κv, C_d(dims))
end


# -----------------------------
# Common utilities / UI
# -----------------------------

dims_of(::SimpleFluid{dims,Tρ,TkT,P})          where {dims,Tρ,TkT,P} = dims
dims_of(::SimpleMixture{dims,Ns,Tρ,TkT,P})        where {dims,Ns,Tρ,TkT,P} = dims
dims_of(::SimpleChargedFluid{dims,Tρ,TkT,P,Tz})   where {dims,Tρ,TkT,P,Tz} = dims
dims_of(::SimpleChargedMixture{dims,Ns,Tρ,TkT,P,Tz}) where {dims,Ns,Tρ,TkT,P,Tz} = dims

number_of_species(::SimpleFluid) = 1
number_of_species(::SimpleMixture{dims,Ns,Tρ,TkT,P})        where {dims,Ns,Tρ,TkT,P} = Ns
number_of_species(::SimpleChargedFluid) = 1
number_of_species(::SimpleChargedMixture{dims,Ns,Tρ,TkT,P,Tz}) where {dims,Ns,Tρ,TkT,P,Tz} = Ns


# base carrier (potential, kBT, ρ)
base_of(sys::SimpleFluid)          = sys
base_of(sys::SimpleMixture)        = sys
base_of(sys::SimpleChargedFluid)   = sys.base
base_of(sys::SimpleChargedMixture) = sys.base

ρ_of(sys::SimpleFluid)           = sys.ρ                      # Number
ρ_of(sys::SimpleMixture)         = sys.ρ                      # Diagonal
ρ_of(sys::SimpleChargedFluid)    = sys.base.ρ
ρ_of(sys::SimpleChargedMixture)  = sys.base.ρ

# charge data: default 'uncharged'
zvec_of(::SimpleFluid)             = SVector{1, Float64}(0.0)
zvec_of(::SimpleMixture{dims,Ns,Tρ,TkT,P}) where {dims,Ns,Tρ,TkT,P} = zero(SVector{dims, TkT}) 
zvec_of(sys::SimpleChargedFluid)   = SVector{1, Float64}(float(sys.z))
zvec_of(sys::SimpleChargedMixture) = sys.z

ℓB_of(::SimpleFluid)              = 0.0
ℓB_of(::SimpleMixture)            = 0.0
ℓB_of(sys::SimpleChargedFluid)    = sys.ℓB
ℓB_of(sys::SimpleChargedMixture)  = sys.ℓB

κsplit_of(::SimpleFluid)              = 1.0
κsplit_of(::SimpleMixture)            = 1.0
κsplit_of(sys::SimpleChargedFluid)    = sys.κ
κsplit_of(sys::SimpleChargedMixture)  = sys.κ

has_coulomb(sys::System) = any(!iszero, zvec_of(sys)) && ℓB_of(sys) > 0

Base.show(io::IO, ::MIME"text/plain", s::SimpleFluid) = begin
    println(io, "$(dims_of(s))D SimpleFluid (single component):")
    println(io, "  ρ = $(s.ρ)")
    println(io, "  kBT = $(s.kBT)")
    println(io, "  potential = $(s.potential)")
end

Base.show(io::IO, ::MIME"text/plain", s::SimpleMixture) = begin
    println(io, "$(dims_of(s))D SimpleMixture (Ns=$(number_of_species(s))):")
    println(io, "  ρ = $(s.ρ)")
    println(io, "  kBT = $(s.kBT)")
    println(io, "  potential = $(s.potential)")
end

Base.show(io::IO, ::MIME"text/plain", s::SimpleChargedFluid) = begin
    println(io, "$(dims_of(s))D SimpleChargedFluid (single component):")
    println(io, "  z = $(s.z), ℓB = $(s.ℓB), κ = $(s.κ)")
    show(io, MIME"text/plain"(), s.base)
end

Base.show(io::IO, ::MIME"text/plain", s::SimpleChargedMixture) = begin
    println(io, "$(dims_of(s))D SimpleChargedMixture (Ns=$(number_of_species(s))):")
    println(io, "  z = $(s.z), ℓB = $(s.ℓB), κ = $(s.κ)")
    show(io, MIME"text/plain"(), s.base)
end


"""
    SimpleLiquid

DEPRECATED constructor alias kept for backward compatibility.

    SimpleLiquid(dims, ρ::Number,  kBT, potential) → SimpleFluid
    SimpleLiquid(dims, ρ::Vector,  kBT, potential) → SimpleMixture

Use `SimpleFluid` or `SimpleMixture` instead.
"""
function SimpleLiquid(dims, ρ::Number, kBT, potential)
    Base.depwarn("`SimpleLiquid(dims, ρ::Number, kBT, potential)` is deprecated; use `SimpleFluid(...)`.", :SimpleLiquid)
    return SimpleFluid(dims, ρ, kBT, potential)
end

function SimpleLiquid(dims, ρ::AbstractVector, kBT, potential)
    Base.depwarn("`SimpleLiquid(dims, ρ::Vector, kBT, potential)` is deprecated; use `SimpleMixture(...)`.", :SimpleLiquid)
    return SimpleMixture(dims, ρ, kBT, potential)
end
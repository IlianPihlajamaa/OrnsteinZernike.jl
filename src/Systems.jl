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

Single-component electrolyte wrapper. Keeps electrostatics (charges, Bjerrum length, split κ)
as **system-level** state; short-range/core physics remain in `potential`.

Fields:
- base::SimpleFluid  — underlying neutral single-component system
- z::Tz              — charge (in units of e)
- ℓB::Float64        — Bjerrum length (in your length units)
- κ::Float64         — Ewald/Gaussian split parameter
"""
mutable struct SimpleChargedFluid{dims,Tρ,TkT,P,Tz} <: AbstractSingleComponent
    base::SimpleFluid{dims,Tρ,TkT,P}
    z::Tz
    ℓB::Float64
    κ::Float64
end

"""
    SimpleChargedMixture{dims,Ns,Tρ,TkT,P,Tz} <: AbstractMixture

Mixture electrolyte wrapper. Electrostatics live on the system; potentials remain short-range.

Fields:
- base::SimpleMixture — underlying neutral mixture
- z::SVector{Ns,Tz}   — per-species charges
- ℓB::Float64         — Bjerrum length
- κ::Float64          — Ewald/Gaussian split parameter
"""
mutable struct SimpleChargedMixture{dims,Ns,Tρ,TkT,P,Tz} <: AbstractMixture
    base::SimpleMixture{dims,Ns,Tρ,TkT,P}
    z::SVector{Ns,Tz}
    ℓB::Float64
    κ::Float64
end

# ---- Dimension constant for the Coulomb kernel's FT:  v̂(k) = (C_d / k^2) / (ε_r β) ----
# (Matches the surface area of the unit sphere: C_3=4π, C_2=2π, C_1=2.)
C_d(d::Integer) = d == 3 ? 4π : d == 2 ? 2π : d == 1 ? 2 : error("Unsupported dimension d=$d")

# ---- Bjerrum length: ℓB = 1 / (C_d * εr * kBT)  (assumes charges measured in |e| units) ----
bjerrum_length(kBT, εr; dims::Integer=3) = 1.0 / (C_d(dims) * εr * kBT)

# ---- Debye screening parameter κ_D  (scalar, single-species) ----
# κ_D^2 = C_d * ℓB * ρ * z^2  ⇒  κ_D = sqrt(C_d * ℓB * ρ * z^2)
debye_kappa(ρ::Number, z::Number, ℓB; dims::Integer=3) =
    sqrt(C_d(dims) * ℓB * ρ * z^2)

# ---- Debye κ_D for mixtures (vector densities) ----
# Accept either a Diagonal of densities or a plain vector; z is the vector of valences.
debye_kappa(ρdiag::Diagonal, z::AbstractVector, ℓB; dims::Integer=3) =
    sqrt(C_d(dims) * ℓB * sum(diag(ρdiag) .* (z .^ 2)))

debye_kappa(ρ::AbstractVector, z::AbstractVector, ℓB; dims::Integer=3) =
    sqrt(C_d(dims) * ℓB * sum(ρ .* (z .^ 2)))


# Charged constructors
function SimpleChargedFluid(base::SimpleFluid; z::Number, εr=78.4, κ=:auto)
    dims = dimensions(base)
    ℓB = bjerrum_length(base.kBT, εr; dims=dims)
    κD = debye_kappa(base.ρ, z, ℓB; dims=dims)
    κv = κ === :auto ? κD : Float64(κ)
    return SimpleChargedFluid{dimensions(base), typeof(base.ρ), typeof(base.kBT), typeof(base.potential), typeof(z)}(base, z, ℓB, κv)
end

function SimpleChargedMixture(base::SimpleMixture; z::AbstractVector, εr=78.4, κ=:auto)
    dims = dimensions(base)
    Ns = length(z)
    zS = SVector{Ns}(z)
    @assert isapprox(sum(diag(base.ρ) .* zS), 0.0; atol=1e-12) "Electroneutrality required: ∑ ρ_i z_i ≈ 0."
    ℓB = bjerrum_length(base.kBT, εr; dims=dims)
    κD = debye_kappa(base.ρ, zS, ℓB; dims=dims)
    κv = κ === :auto ? κD : Float64(κ)
    return SimpleChargedMixture{dimensions(base), Ns, typeof(base.ρ), typeof(base.kBT), typeof(base.potential), eltype(zS)}(base, zS, ℓB, κv)
end


# -----------------------------
# Common utilities / UI
# -----------------------------

dimensions(::SimpleFluid{dims,Tρ,TkT,P})          where {dims,Tρ,TkT,P} = dims
dimensions(::SimpleMixture{dims,Ns,Tρ,TkT,P})        where {dims,Ns,Tρ,TkT,P} = dims
dimensions(::SimpleChargedFluid{dims,Tρ,TkT,P,Tz})   where {dims,Tρ,TkT,P,Tz} = dims
dimensions(::SimpleChargedMixture{dims,Ns,Tρ,TkT,P,Tz}) where {dims,Ns,Tρ,TkT,P,Tz} = dims

number_of_species(::SimpleFluid) = 1
number_of_species(::SimpleMixture{dims,Ns,Tρ,TkT,P})        where {dims,Ns,Tρ,TkT,P} = Ns
number_of_species(::SimpleChargedFluid) = 1
number_of_species(::SimpleChargedMixture{dims,Ns,Tρ,TkT,P,Tz}) where {dims,Ns,Tρ,TkT,P,Tz} = Ns



Base.show(io::IO, ::MIME"text/plain", s::SimpleFluid) = begin
    println(io, "$(dimensions(s))D SimpleFluid (single component):")
    println(io, "  ρ = $(s.ρ)")
    println(io, "  kBT = $(s.kBT)")
    println(io, "  potential = $(s.potential)")
end

Base.show(io::IO, ::MIME"text/plain", s::SimpleMixture) = begin
    println(io, "$(dimensions(s))D SimpleMixture (Ns=$(number_of_species(s))):")
    println(io, "  ρ = $(s.ρ)")
    println(io, "  kBT = $(s.kBT)")
    println(io, "  potential = $(s.potential)")
end

Base.show(io::IO, ::MIME"text/plain", s::SimpleChargedFluid) = begin
    println(io, "$(dimensions(s))D SimpleChargedFluid (single component):")
    println(io, "  z = $(s.z), ℓB = $(s.ℓB), κ = $(s.κ)")
    show(io, MIME"text/plain"(), s.base)
end

Base.show(io::IO, ::MIME"text/plain", s::SimpleChargedMixture) = begin
    println(io, "$(dimensions(s))D SimpleChargedMixture (Ns=$(number_of_species(s))):")
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
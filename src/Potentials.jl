
"""
    Potential

Abstract potential type
"""
abstract type Potential end


"""
    HardSpheres


Implements the hard-sphere pair interaction for single component systems \$ u(r) = \\infty\$ for \$r < 1\$ and \$u(r) = 0\$ for \$r > 1\$,
or \$u_{ij}(r) = \\infty\$ for \$r < D_{ij}\$ and \$u_{ij}(r) = 0\$ for \$r > D_{ij}\$ for mixtures.

For mixtures expects either a vector \$D_i\$ of diameters for each of the species in which case an additive mixing rule is used \$\\left(D_{ij} = (D_{i}+D_{j})/2\\right)\$ 
or a matrix \$D_ij\$ of pair diameters.

Example:
```example 1
potential = HardSpheres(1.0)
```

Example:
```example 2
potential = HardSpheres([0.8, 0.9, 1.0])
```

```example 3
Dij = rand(3,3)
potential = HardSpheres(Dij)
```

"""
struct HardSpheres{T} <: Potential 
    D::T

    HardSpheres(D::Number) = new{typeof(D)}(D)

    function HardSpheres(D::AbstractVector{T}) where T 
        Ds = SVector{length(D),T}(D)
        Dij = (Ds .+ Ds')/2
        T2 = typeof(Dij)
        new{T2}(Dij)
    end
    
    function HardSpheres(D::AbstractMatrix{T}) where T 
        @assert size(D, 1) == size(D, 2) "matrix of pair diameters must be square"
        Ns = size(D, 1)
        Ds = SMatrix{Ns, Ns}(D)
        T2 = typeof(Ds)
        new{T2}(Ds)
    end
end

function evaluate_potential(potential::HardSpheres{T}, r::Number) where T
    Dij = potential.D
    pot = @. ifelse(r < Dij, Inf64, 0.0)
    return pot
end



"""
    LennardJones

Implements the Lennard-Jones pair interaction \$u(r) = 4\\epsilon [(\\sigma/r)^{12} - (\\sigma/r)^6]\$.

Expects values `ϵ` and `σ`, which respecively are the strength of the potential and particle size. 

Example:
```julia
potential = LennardJones(1.0, 2.0)
```
"""
struct LennardJones{T1, T2} <: Potential 
    ϵ::T1
    σ::T2
end


function evaluate_potential(potential::LennardJones, r::Number)
    ϵ = potential.ϵ
    σ = potential.σ
    return @. 4ϵ * ((σ/r)^12 - (σ/r)^6)
end


"""
    PowerLaw

Implements the power law pair interaction \$u(r) = \\epsilon (\\sigma/r)^{n}\$.

Expects values `ϵ`, `σ`, and `n`, which respecively are the strength of the potential and particle size. 

Example:
```julia
potential = PowerLaw(1.0, 2.0, 8)
```
"""
struct PowerLaw{T1, T2, T3} <: Potential 
    ϵ::T1
    σ::T2
    n::T3
end

function evaluate_potential(potential::PowerLaw, r::Number)
    ϵ = potential.ϵ
    σ = potential.σ
    n = potential.n
    return @. ϵ * (σ/r)^n 
end

"""
    WCA

Implements the Weeks-Chandler-Andersen pair interaction \$u(r) = 4\\epsilon [(\\sigma/r)^{12} - (\\sigma/r)^6 - 1]\$ for \$r>2^{1/2}\\sigma\$ and \$0\$ otherwise.

Expects values `ϵ`, `σ`, and `n`, which respecively are the strength of the potential and particle size. 

Example:
```julia
potential = WCA(1.0, 2.0)
```
"""
struct WCA{T1, T2} <: Potential 
    ϵ::T1
    σ::T2
end

function evaluate_potential(potential::WCA, r::Number)
    ϵ = potential.ϵ
    σ = potential.σ
    return @. ifelse(r > σ * 2^(1/6), zero(ϵ), 4ϵ * ((σ/r)^12 - (σ/r)^6 - one(ϵ)) )
end

"""
    CustomPotential

Implements a potential that evaluates a user defined function.

Expects values `f`, and `p`, which respecively are a callable and a list of parameters.
The function should be called `f(r::Number, p)` and it should produce either a `Number`,
in the case of a single-component system, or an `SMatrix`, in the case of a multicomponent system. 

Example:
```julia
f = (r, p) -> 4*p[1]*((p[2]/r)^12 -  (p[2]/r)^6)
potential = CustomPotential(f, (1.0, 1.0))
```
"""
struct CustomPotential{T1, T2} <: Potential 
    f::T1
    p::T2
end

function evaluate_potential(potential::CustomPotential, r::Number)
    return potential.f(r, potential.p)
end


"""
exp(- beta * u) - 1.
"""
function find_mayer_f_function(system::SimpleLiquid{dims, species, T1, T2, P}, r::Number) where {dims, species, T1, T2, P}
    U = evaluate_potential(system.potential, r)
    β = 1/system.kBT
    f = @. exp(-β*U) - 1.0
    return f
end

function find_mayer_f_function(system::SimpleLiquid{dims, species, T1, T2, P}, r::AbstractArray) where {dims, species, T1, T2, P}
    return find_mayer_f_function.((system, ), r)
end

function evaluate_potential(potential::Potential, r::AbstractArray)
    return evaluate_potential.((potential, ), r)
end

evaluate_potential_derivative(potential::HardSpheres, ::Number) = zero(typeof(potential.D))

function evaluate_potential_derivative(potential::Potential, r::AbstractVector)
    return evaluate_potential_derivative.((potential, ), r)
end

function evaluate_potential_derivative(potential::Potential, r::Number)
    #check for discontinuities
    ϵ = sqrt(eps(r))


    discs = discontinuities(potential)
    for discontinuity in discs
        if any(abs.(discontinuity .- r) .< ϵ)
            error("Trying to evaluate the derivative of the potential at the discontinuity. To fix, define a specialized method for `evaluate_potential_derivative(potential::MyPotential, r)`")
        end
    end 

    u2 = evaluate_potential(potential, r+ϵ)
    u1 = evaluate_potential(potential, r-ϵ)
    if isinf(u2) && isinf(u1)
        return zero(u2)
    end
    return (u2-u1)/(2ϵ)
end

function discontinuities(::Potential)
    return Float64[]
end

function discontinuities(p::HardSpheres{T}) where T<:AbstractArray
    return p.D[:]
end
function discontinuities(p::HardSpheres{T}) where T<:Number
    return [p.D]
end

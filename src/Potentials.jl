
"""
Potential

Abstract potential type
"""
abstract type Potential end


"""
    SingleComponentHardSpheres

Implements the hard-sphere pair interaction \$ u(r) = \\infty\$ for \$r < 1\$ and \$u(r) = 0\$ for \$r > 1\$.

Example:
```julia
closure = SingleComponentHardSpheres()
```
"""
struct SingleComponentHardSpheres <: Potential end

function evaluate_potential(::SingleComponentHardSpheres, r::Number)
    return ifelse(r < oneunit(r), Inf64, 0.0)*oneunit(r)
end


"""
    MultiComponentHardSpheres

Implements the hard-sphere pair interaction \$u_{ij}(r) = \\infty\$ for \$r < D_{ij}\$ and \$u_{ij}(r) = 0\$ for \$r > D_{ij}\$.

Expects a vector \$D_i\$ of diameters for each of the species. An additive mixing rule is used \$\\left(D_{ij} = (D_{i}+D_{j})/2\\right)\$.

Example:
```julia
closure = MultiComponentHardSpheres([0.8, 0.9, 1.0])
```
"""
struct MultiComponentHardSpheres{N_components, T} <: Potential 
    D::SVector{N_components, T}
end

function MultiComponentHardSpheres(D::AbstractVector{T}) where T<:Number
    N = length(D)
    return MultiComponentHardSpheres{N,T}(SVector{N,T}(D))
end

function evaluate_potential(potential::MultiComponentHardSpheres{N,T}, r::Number) where {N,T}
    pot = MMatrix{N, N, T, N*N}(undef)
    D = potential.D
    for species1 = 1:N
        for species2 = 1:N
            pot[species2, species1] = ifelse(r < (D[species1]+D[species2])/2, Inf64, 0.0)*oneunit(r)
        end
    end
    return SMatrix(pot)
end


"""
    SingleComponentLennardJones

Implements the Lennard-Jones pair interaction \$u(r) = 4\\epsilon [(\\sigma/r)^{12} - (\\sigma/r)^6]\$.

Expects values `ϵ` and `σ`, which respecively are the strength of the potential and particle size. 

Example:
```julia
closure = SingleComponentLennardJones(1.0, 2.0)
```
"""
struct SingleComponentLennardJones{T1, T2} <: Potential 
    ϵ::T1
    σ::T2
end


function evaluate_potential(potential::SingleComponentLennardJones, r::Number)
    ϵ = potential.ϵ
    σ = potential.σ
    return 4ϵ * ((σ/r)^12 - (σ/r)^6)
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

evaluate_potential_derivative(potential::SingleComponentHardSpheres, r) = 0.0

function evaluate_potential_derivative(potential::Potential, r)
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

function discontinuities(::SingleComponentHardSpheres)
    return [1.0]
end

function discontinuities(potential::MultiComponentHardSpheres)
    D = potential.D
    T = eltype(D)
    disc = MMatrix{N, N, T, N*N}(undef)
    for species1 = 1:N
        for species2 = 1:N
            disc[species2, species1] =  (D[species1]+D[species2])/2
        end
    end
    return [SMatrix(disc)]
end

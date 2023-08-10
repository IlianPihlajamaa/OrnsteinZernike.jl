
"""
Potential

Abstract potential type
"""
abstract type Potential end


"""
    SingleComponentHardSpheres

Implements the hard-sphere pair interaction u(r) = inf r<1 and u(r) = 0 r>1.

Example:
`closure = SingleComponentHardSpheres()`
"""
struct SingleComponentHardSpheres <: Potential end

function evaluate_potential(::SingleComponentHardSpheres, r::Number)
    return ifelse(r < oneunit(r), Inf64, 0.0)*oneunit(r)
end


"""
    MultiComponentHardSpheres

Implements the hard-sphere pair interaction uᵢⱼ(r) = inf r<Dᵢⱼ and uᵢⱼ(r) = 0 r>Dᵢⱼ for a multicomponent system.

Expects a vector Dᵢ of diameters for each of the species. An additive mixing rule is used (Dᵢⱼ = (Dᵢ+Dⱼ)/2).

Example:
`closure = MultiComponentHardSpheres([0.8, 0.9, 1.0])`
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
exp(- beta * u) - 1.
"""
function find_mayer_f_function(system::SimpleLiquid{dims, species, T1, T2, P}, r::Number, β::Number) where {dims, species, T1, T2, P}
    U = evaluate_potential(system.potential, r)
    return @. exp(-β*U) - 1.0
end

function find_mayer_f_function(system::SimpleLiquid{dims, species, T1, T2, P}, r::AbstractArray, β::Number) where {dims, species, T1, T2, P}
    return find_mayer_f_function.((system, ), r, β)
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
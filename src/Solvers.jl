abstract type Method end

struct Exact <: Method 
    M::Int
    dr::Float64
end

function Exact(;  M=2^10, dr = sqrt(π/(M+1))/(2π))
    return Exact(M, dr)
end


struct FourierIteration <: Method 
    M::Int
    dr::Float64
    mixing_parameter::Float64
    max_iterations::Int64
    tolerance::Float64
    verbose::Bool
end

function FourierIteration(; mixing_parameter=0.5, max_iterations=10^5, tolerance=10^-6, verbose=true, M=2^10, dr=sqrt(π/(M+1))/(2π))
    @assert max_iterations > 0 
    @assert tolerance > 0 
    @assert 0 <= mixing_parameter <= 1
    @assert M > 0
    FourierIteration(M, dr, mixing_parameter, max_iterations, tolerance, verbose)
end

struct NgIteration <: Method 
    M::Int
    dr::Float64
    N_stages::Int64
    max_iterations::Int64
    tolerance::Float64
    verbose::Bool
end

function NgIteration(; N_stages=3, max_iterations=10^3, tolerance=10^-6, verbose=true, M=2^10, dr=sqrt(π/(M+1))/(2π))
    @assert max_iterations > 0 
    @assert tolerance > 0 
    @assert N_stages > 0
    @assert M > 0
    NgIteration(M, dr, N_stages, max_iterations, tolerance, verbose)
end

struct DensityRamp{T<:Method, T2<:AbstractVector} <: Method 
    method::T
    densities::T2
    verbose::Bool
end

function DensityRamp(method, densities; verbose=true)
    DensityRamp(method, densities, verbose)
end

defaultsolver() = NgIteration()

function solve(system::SimpleLiquid, closure::Closure)
    solve(system, closure, defaultsolver())
end

ndims(::SimpleLiquid{Ndims, species, T1, T2, P}) where {Ndims, species, T1, T2, P} = Ndims

function solve(liquid::SimpleLiquid, closure::Closure, ::Exact)
    P = typeof(liquid.potential)
    error("The potential $(P) with closure $(typeof(closure)) has no implemented exact solution in $(ndims(liquid)) dimensions.")
end

function compute_error(y1::Vector{T}, y2::Vector{T}) where T
    error = zero(T)
    @assert length(y1) == length(y2)
    for i in eachindex(y1)
        error += (y1[i] - y2[i]) .^ 2
    end
    return sqrt(sum(error))
end

function find_g_from_c_and_Γ(c::Vector{T}, Γ::Vector{T}, r) where T<:Number
    return  @. c + Γ / r + one(T)
end

function find_g_from_c_and_Γ(c::Vector{T}, Γ::Vector{T}, r) where T<:AbstractMatrix
    g = copy(c)
    for i in eachindex(c, Γ)
        g[i] = c[i] + Γ[i] / r[i] + ones(T)
    end 
    return g
end

function find_S_from_ĉ_and_ρ(ĉ::Vector{T}, ρ) where T<:Number
    return  @. one(T) + ρ * ĉ / (one(T) - ρ * ĉ)
end

function find_S_from_ĉ_and_ρ(ĉ::Vector{T}, ρ) where T<:AbstractMatrix
    S = copy(ĉ)
    ρtot = sum(ρ.diag)
    x = ρ.diag/ρtot
    for i in eachindex(ĉ)
        Si = inv(one(T) ./ x .- sum(ρ)*ĉ[i])
        S[i] = Si
    end 
    return S
end


function construct_r_and_k_grid(::SimpleLiquid{3, species, T1, T2, P}, method::Union{Exact, NgIteration, FourierIteration}) where {species, T1, T2, P}
    M = method.M
    dr = method.dr
    dk = π/(M+1)/dr
    r = collect((1:M)*dr)
    k = collect((1:M)*dk)
    return r, k
end


include("Solvers/Exact.jl")
include("Solvers/FourierIteration.jl")
include("Solvers/NgIteration.jl")
include("Solvers/DensityRamp.jl")




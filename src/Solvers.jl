"""
Method

Abstract method type
"""
abstract type Method end


"""
    Exact <: Method

Solves the system exactly. This is only implemented for specific systems.

Construct using `Exact(;  M=2^10, dr = sqrt(π/(M+1))/(2π))`
Here, `M` is the number of points that the exact solution is evaluated on, and `dr` is the grid spacing.
These are used to perform fourier transformations.

Examples

`method = Exact()`
`method = Exact(M=1000)`
`method = Exact(M=1000, dr=0.01)`
"""
struct Exact <: Method 
    M::Int
    dr::Float64
end

function Exact(;  M=2^10, dr = sqrt(π/(M+1))/(2π))
    return Exact(M, dr)
end

"""
FourierIteration <: Method

Solves the system by recursive iteration in Fourier Space. Essentially, the algorithm is:

1. guess an initial γ(r)
2. find c(r) using the closure relation
3. fourier transform to get ĉ(k)
4. find γ(k) using the OZ-eq in k-space
5. compute γ(r) with a inverse fourier transform
6. compare with previous value, if not converged go to 2.

Optionally, a mixing rule is used to mix the new and previous iteration of c(r) in step 2. 

Arguments:
- M::Int: number of points discretize the solution on 
- dr::Float64: grid spacing in real space
- mixing_parameter::Float64: mixing parameter for iteration mixing. A value of 1 is no mixing. Must be between 0 and 1. 
- max_iterations::Int64: maximal number of iterations 
- tolerance::Float64: tolerance to be reached
- verbose::Bool: whether or not to print convergence information

default: FourierIteration(; mixing_parameter=0.5, max_iterations=10^5, tolerance=10^-6, verbose=true, M=2^10, dr=sqrt(π/(M+1))/(2π))
"""
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

"""
NgIteration <: Method

Solves the system by recursive iteration in Fourier Space, and uses the Ng acceleration method. Essentially, the algorithm is:

1. guess an initial γ(r)
2. find c(r) using the closure relation
3. fourier transform to get ĉ(k)
4. find γ(k) using the OZ-eq in k-space
5. compute γ(r) with a inverse fourier transform
6. use Ng's method to provide a next guess for γ
7. compare with previous value, if not converged go to 2.

Arguments:
- M::Int: number of points discretize the solution on 
- dr::Float64: grid spacing in real space
- N_stages::Int: Number of previous values to take into account for step 6. A higher number should lead to faster convergence, yet more computation time per iteration.
- max_iterations::Int64: maximal number of iterations 
- tolerance::Float64: tolerance to be reached
- verbose::Bool: whether or not to print convergence information

default: NgIteration(; N_stages=3, max_iterations=10^3, tolerance=10^-6, verbose=true, M=2^10, dr=sqrt(π/(M+1))/(2π))
"""
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

"""
DensityRamp <: Method

Solves the system by iteratively solving systems of increasing density, using the previous solution as initial guess at a higher density. This may help deal with convergence issues

Arguments
- method: method by which to solve the system for individual densities.
- densities: densities to consider. Must be a vector of increasing values.
- verbose: whether to print information.

Example:
`DensityRamp(NgIteration(), [0.1, 0.3, 0.4]; verbose=false)`

"""
struct DensityRamp{T<:Method, T2<:AbstractVector} <: Method 
    method::T
    densities::T2
    verbose::Bool
end

function DensityRamp(method, densities; verbose=true)
    @assert issorted(sum.(densities))
    DensityRamp(method, densities, verbose)
end

defaultsolver() = NgIteration()


"""
    solve(system::SimpleLiquid, closure::Closure, method::Method)

Solves the system `system` using the closure `closure` with method `method`.

    solve(system::SimpleLiquid, closure::Closure)

Solves the system `system` using the closure `closure` with the default method `NgIteration()`.
"""
function solve end

function solve(system::SimpleLiquid, closure::Closure)
    solve(system, closure, defaultsolver())
end

ndims(::SimpleLiquid{Ndims, species, T1, T2, P}) where {Ndims, species, T1, T2, P} = Ndims

function solve(system::SimpleLiquid, closure::Closure, ::Exact)
    P = typeof(system.potential)
    error("The potential $(P) with closure $(typeof(closure)) has no implemented exact solution in $(ndims(system)) dimensions.")
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




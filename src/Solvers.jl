
"""
    Method

Abstract method type
"""
abstract type Method end


"""
    Exact <: Method

Solves the system exactly. This is only implemented for specific systems.

Right now, the implemented exact methods are
- three-dimensional single-component system of hard spheres with the Percus Yevick closure [1]
- three-dimensional multi-component system of additive hard spheres with the Percus Yevick closure [2]
- one-dimensional single-component system of hard spheres with the Percus Yevick closure [3]
- five-dimensional single-component system of hard spheres with the Percus Yevick closure [3]

Construct using 

```Exact(;  M=1000, dr = 10.0/M)```

Here, `M` is the number of points that the exact solution is evaluated on, and `dr` is the grid spacing.
These are used to perform Fourier transformations.

Examples

`method = Exact()`

`method = Exact(M=1000)`

`method = Exact(M=1000, dr=0.01)`

References:
1. Wertheim, M. S. "Exact solution of the Percus-Yevick integral equation for hard spheres." Physical Review Letters 10.8 (1963): 321.
2. Baxter, R. J. "Ornstein–Zernike relation and Percus–Yevick approximation for fluid mixtures." The Journal of Chemical Physics 52.9 (1970): 4559-4562.
3. Leutheusser, E. "Exact solution of the Percus-Yevick equation for a hard-core fluid in odd dimensions." Physica A: Statistical Mechanics and its Applications 127.3 (1984): 667-676.
"""
struct Exact <: Method
    M::Int
    dr::Float64
end

function Exact(; M=1000, dr=10.0 / M)
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
- `M::Int`: number of points discretize the solution on 
- `dr::Float64`: grid spacing in real space
- `mixing_parameter::Float64`: mixing parameter for iteration mixing. A value of 1 is no mixing. Must be between 0 and 1. 
- `max_iterations::Int64`: maximal number of iterations 
- `tolerance::Float64`: tolerance to be reached
- `verbose::Bool`: whether or not to print convergence information

Default: `FourierIteration(; mixing_parameter=0.5, max_iterations=10^5, tolerance=10^-6, verbose=true, M=1000, dr=10.0/M)`
"""
struct FourierIteration <: Method
    M::Int
    dr::Float64
    mixing_parameter::Float64
    max_iterations::Int64
    tolerance::Float64
    verbose::Bool
end

function FourierIteration(; mixing_parameter=0.5, max_iterations=10^5, tolerance=10^-10, verbose=true, M=1000, dr=10.0 / M)
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
- `M::Int`: number of points discretize the solution on 
- `dr::Float64`: grid spacing in real space
- `N_stages::Int`: Number of previous values to take into account for step 6. A higher number should lead to faster convergence, yet more computation time per iteration.
- `max_iterations::Int64`: maximal number of iterations 
- `tolerance::Float64`: tolerance to be reached
- `verbose::Bool`: whether or not to print convergence information

Default: `NgIteration(; N_stages=3, max_iterations=10^3, tolerance=10^-6, verbose=true, M=1000, dr=10.0/M)`

References:
Ng, K. C. (1974). Hypernetted chain solutions for the classical one‐component plasma up to Γ= 7000. The Journal of Chemical Physics, 61(7), 2680-2689.
"""
struct NgIteration <: Method
    M::Int
    dr::Float64
    N_stages::Int64
    max_iterations::Int64
    tolerance::Float64
    verbose::Bool
end

function NgIteration(; N_stages=3, max_iterations=10^3, tolerance=10^-10, verbose=true, M=1000, dr=10.0 / M)
    @assert max_iterations > 0
    @assert tolerance > 0
    @assert N_stages > 0
    @assert M > 0
    NgIteration(M, dr, N_stages, max_iterations, tolerance, verbose)
end

"""
    DensityRamp <: Method

Solves the system by iteratively solving systems of increasing density, using the previous solution as initial guess at a higher density. This may help deal with convergence issues.

Arguments
- `method`: method by which to solve the system for individual densities.
- `densities`: densities to consider. Must be a vector of increasing values.
- `verbose`: whether to print information.

Example:
`DensityRamp(NgIteration(), [0.1, 0.3, 0.4]; verbose=false)`

"""
struct DensityRamp{T<:Method,T2<:AbstractVector} <: Method
    method::T
    densities::T2
    verbose::Bool
end

function DensityRamp(method::Method, densities::AbstractVector{T}; verbose=true) where {T}
    if method isa DensityRamp
        error("Nesting DensityRamp methods is not supported.")
    end
    if T <: Number
        @assert issorted(densities)
        return DensityRamp(method, densities, verbose)
    elseif T <: AbstractVector
        @assert issorted(sum.(densities))
        Ns = length(densities[1])
        @assert allequal(length.(densities))
        densities = Diagonal.(SVector{Ns}.(densities))
        return DensityRamp(method, densities, verbose)
    elseif T <: Diagonal
        @assert issorted(sum.(densities))
        return DensityRamp(method, densities, verbose)
    end
    error("Invalid type for densities.")
end

"""
    TemperatureRamp <: Method

Solves the system by iteratively solving systems of decreasing Temperature, using the previous solution as initial guess at a lower temperature. This may help deal with convergence issues.

Arguments
- `method`: method by which to solve the system for individual temperatures.
- `temperatures`: temperatures to consider. Must be a vector of increasing values.
- `verbose`: whether to print information.

Example:
`TemperatureRamp(NgIteration(), [0.1, 0.3, 0.4]; verbose=false)`

"""
struct TemperatureRamp{T<:Method,T2<:AbstractVector} <: Method
    method::T
    temperatures::T2
    verbose::Bool
end

function TemperatureRamp(method, temperatures::AbstractVector{T}; verbose=true) where {T}
    @assert issorted(temperatures, rev=true)
    return TemperatureRamp(method, temperatures, verbose)
end


defaultsolver() = NgIteration()


"""
    solve(system::System, closure::Closure, method::Method)

Solves the system `system` using the closure `closure` with method `method`.

    solve(system::System, closure::Closure)

Solves the system `system` using the closure `closure` with the default method `NgIteration()`.
"""
function solve end

function solve(system::System, closure::Closure; gamma_0=nothing)
    solve(system, closure, defaultsolver(); gamma_0=gamma_0)
end

function solve(system::System, closure::Closure, ::Exact)
    P = typeof(system.potential)
    error("The potential $(P) with closure $(typeof(closure)) has no implemented exact solution in $(dimensions(system)) dimensions.")
end

function compute_error(y1::AbstractArray, y2::AbstractArray) 
    T = eltype(y1)
    @assert T == eltype(y2)
    error = zero(T)
    @assert length(y1) == length(y2)
    for i in eachindex(y1)
        error += (y1[i] - y2[i]) .^ 2
    end
    return sqrt(sum(error))
end

struct OZSolverCache{T, S, F}
    mayer_f::Vector{T}
    fourierplan::F
    r::Vector{S}
    k::Vector{S}
    βu_dispersion_tail::Vector{T}
    βu::Vector{T}
    Γhat::Vector{T}
    C::Vector{T}
    Ĉ::Vector{T}
    Γ_new::Vector{T}
end

function OZSolverCache(system, method)
    r = method.dr * (1:method.M) |> collect
    βu1, _ = evaluate_long_range_potential(system.potential, system.kBT, r[1])
    elementtype = typeof(r[1] .* system.kBT .* system.ρ .* βu1)
    mayer_f = zeros(elementtype, length(r))
    fourierplan = get_fourier_plan(system, method, mayer_f)
    r .= fourierplan.r # in the case that dims != 3, we need to use the right grid
    k = fourierplan.k
    βu, βu_dispersion_tail = evaluate_long_range_potential(system.potential, system.kBT, r)
    mayer_f .= find_mayer_f_function.((system,), βu)

    T = eltype(mayer_f);  TT = eltype(T)
    Γhat = copy(mayer_f); Γ_new = copy(mayer_f)
    C = copy(mayer_f); Ĉ = copy(mayer_f)
    return OZSolverCache(mayer_f, fourierplan, r, k, βu_dispersion_tail, βu, Γhat, C, Ĉ, Γ_new)
end

# gets one iteration step of the OZ equation in k-space
# computes Γ_new from C
# overwrites Ĉ, Γhat, Γ_new
function oz_iteration_step!(C, Ĉ, Γhat, Γ_new, ρ, plan)
    k = plan.k
    fourier!(Ĉ, C, plan)
    @inbounds for ik in eachindex(Γhat, Ĉ)
        Γhat[ik] = (k[ik] * I - Ĉ[ik] * ρ) \ (Ĉ[ik] * ρ * Ĉ[ik])
    end
    inverse_fourier!(Γ_new, Γhat, plan)
    return 
end




include("Solvers/Exact.jl")
include("Solvers/FourierIteration.jl")
include("Solvers/NgIteration.jl")
include("Solvers/DensityRamp.jl")
include("Solvers/TemperatureRamp.jl")

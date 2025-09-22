
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

function Base.show(io::IO, ::MIME"text/plain", p::HardSpheres)
    println(io, "HardSpheres($(p.D))")
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

function find_mayer_f_function(βu::Number)
    return exp(-βu) - one(βu)
end

function find_mayer_f_function(βu::AbstractArray)
    return @. exp(-βu) - 1.0
end
"""
exp(- beta * u) - 1.
"""
function find_mayer_f_function(::SimpleUnchargedSystem, βU)
    f = @. exp(-βU) - 1.0
    return f
end

function find_mayer_f_function(system::SimpleUnchargedSystem, r::AbstractArray, βU::AbstractArray)
    return find_mayer_f_function.((system, ), r, βU)
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

########################################
# Yukawa / Screened-Coulomb (mixtures) #
########################################

"""
    Yukawa

Screened-Coulomb / Yukawa pair interaction

    u(r) = A * exp(-κ r) / r

Constructors:
- `Yukawa(A::Number, κ::Number)` — single-component
- `Yukawa(q::AbstractVector, κ::Number)` — mixture from "charges": A_ij = q_i q_j
- `Yukawa(Aij::AbstractMatrix, κ::Number)` — mixture with explicit pair amplitudes

Example (single component):
```julia
potential = Yukawa(1.0, 2.0)  # A=1, κ=2
```

Example (mixture from charges):
```julia
q  = [1.0, -1.0, 2.0]
pot = Yukawa(q, 1.5)
```

Example (mixture with explicit A_ij):
```julia
Aij = [1.0 0.2; 0.2 0.5]
pot = Yukawa(Aij, 1.0)
```
"""
struct Yukawa{TA,TK} <: Potential
    A::TA
    κ::TK
    Yukawa(A::Number, κ::Number) = new{typeof(A),typeof(κ)}(A, κ)

    function Yukawa(q::AbstractVector{T}, κ::Number) where {T<:Number}
        qs = SVector{length(q),T}(q)
        Aij = qs * qs'
        return new{typeof(Aij),typeof(κ)}(Aij, κ)
    end

    function Yukawa(Aij::AbstractMatrix{T}, κ::Number) where {T<:Number}
        Ns = size(Aij, 1); @assert size(Aij,1) == size(Aij,2)
        A = SMatrix{Ns,Ns,T}(Aij)
        return new{typeof(A),typeof(κ)}(A, κ)
    end
end

function evaluate_potential(p::Yukawa, r::Number)
    A, κ = p.A, p.κ
    return A .* exp.(-κ .* r) ./ r
end

# du/dr = A e^{-κ r} * (-(κ r + 1)) / r^2   (elementwise w.r.t. pair table A)
function evaluate_potential_derivative(p::Yukawa, r::Number)
    A, κ = p.A, p.κ
    return A .* exp.(-κ .* r) .* (-(κ .* r .+ 1.0)) ./ (r.^2)
end

discontinuities(::Yukawa) = Float64[]


#############################################
# Gaussian Core Model (GCM) — with mixtures #
#############################################

"""
    GaussianCore

Gaussian core (ultrasoft) pair interaction

    u(r) = ϵ * exp(-(r/σ)^2)

Constructors:
- `GaussianCore(ϵ::Number, σ::Number)` — single-component
- `GaussianCore(ϵ::AbstractVector, σ::AbstractVector)` — mixture with
   σ_ij = (σ_i + σ_j)/2 (additive), ϵ_ij = √(ϵ_i ϵ_j) (geometric mean)
- `GaussianCore(ϵij::AbstractMatrix, σij::AbstractMatrix)` — explicit pair tables

Example:
```julia
potential = GaussianCore(1.0, 1.5)
```

Example (mixture with mixing rules):
```julia
eps = [1.0, 0.8]
sig = [1.0, 1.2]
potential = GaussianCore(eps, sig)
```
"""
struct GaussianCore{TE,TS} <: Potential
    ϵ::TE
    σ::TS

    GaussianCore(ϵ::Number, σ::Number) = new{typeof(ϵ),typeof(σ)}(ϵ, σ)

    function GaussianCore(ϵ::AbstractVector{T1}, σ::AbstractVector{T2}) where {T1<:Number,T2<:Number}
        @assert length(ϵ) == length(σ)
        N  = length(ϵ)
        ϵs = SVector{N,T1}(ϵ)
        σs = SVector{N,T2}(σ)
        ϵij = sqrt.(ϵs * ϵs')
        σij = (σs .+ σs') / 2
        return new{typeof(ϵij), typeof(σij)}(ϵij, σij)
    end

    function GaussianCore(ϵij::AbstractMatrix{T1}, σij::AbstractMatrix{T2}) where {T1<:Number,T2<:Number}
        @assert size(ϵij,1)==size(ϵij,2)==size(σij,1)==size(σij,2)
        N  = size(ϵij,1)
        ϵM = SMatrix{N,N,T1}(ϵij)
        σM = SMatrix{N,N,T2}(σij)
        return new{typeof(ϵM), typeof(σM)}(ϵM, σM)
    end
end

function evaluate_potential(p::GaussianCore, r::Number)
    ϵ, σ = p.ϵ, p.σ
    return ϵ .* exp.(- (r ./ σ).^2)
end

# du/dr = u(r) * (-2 r / σ^2)
function evaluate_potential_derivative(p::GaussianCore, r::Number)
    ϵ, σ = p.ϵ, p.σ
    u = ϵ .* exp.(- (r ./ σ).^2)
    return u .* (-2 .* r ./ (σ.^2))
end

discontinuities(::GaussianCore) = Float64[]


##############################
# Square-Well (with mixtures)#
##############################

"""
    SquareWell

Square-well pair interaction:

    u(r) = ∞                     if r < σ
         = -ϵ                    if σ ≤ r ≤ λσ
         = 0                     if r > λσ

Constructors:
- `SquareWell(σ::Number, ϵ::Number, λ::Number)` — single-component
- `SquareWell(σ::AbstractVector, ϵ::AbstractVector, λ::Number)` — mixture with
    σ_ij = (σ_i + σ_j)/2, ϵ_ij = √(ϵ_i ϵ_j)
- `SquareWell(σij::AbstractMatrix, ϵij::AbstractMatrix, λ::Number)` — explicit pair tables

Example:
```julia
potential = SquareWell(1.0, 1.0, 1.5)
```

Example (mixture with mixing rules):
```julia
sig = [0.9, 1.1, 1.0]
eps = [1.0, 0.8, 1.2]
potential = SquareWell(sig, eps, 1.5)
```
"""
struct SquareWell{Tσ,Tϵ,Tλ} <: Potential
    σ::Tσ
    ϵ::Tϵ
    λ::Tλ
end

# Single-component
SquareWell(σ::Number, ϵ::Number, λ::Number) = SquareWell{typeof(σ),typeof(ϵ),typeof(λ)}(σ, ϵ, λ)

# Mixture from vectors
function SquareWell(σ::AbstractVector{Tσ}, ϵ::AbstractVector{Tϵ}, λ::Number) where {Tσ<:Number,Tϵ<:Number}
    @assert length(σ) == length(ϵ)
    N  = length(σ)
    σs = SVector{N,Tσ}(σ)
    ϵs = SVector{N,Tϵ}(ϵ)
    σij = (σs .+ σs') / 2
    ϵij = sqrt.(ϵs * ϵs')
    return SquareWell{typeof(σij),typeof(ϵij),typeof(λ)}(σij, ϵij, λ)
end

# Mixture from explicit pair tables
function SquareWell(σij::AbstractMatrix{Tσ}, ϵij::AbstractMatrix{Tϵ}, λ::Number) where {Tσ<:Number,Tϵ<:Number}
    @assert size(σij,1) == size(σij,2) == size(ϵij,1) == size(ϵij,2)
    N  = size(σij,1)
    σM = SMatrix{N,N,Tσ}(σij)
    ϵM = SMatrix{N,N,Tϵ}(ϵij)
    return SquareWell{typeof(σM),typeof(ϵM),typeof(λ)}(σM, ϵM, λ)
end

function evaluate_potential(p::SquareWell, r::Number)
    σ, ϵ, λ = p.σ, p.ϵ, p.λ
    # broadcasted piecewise: works for scalar or pair tables
    return ifelse.(r .< σ, Inf, ifelse.(r .<= λ .* σ, -ϵ, 0.0))
end

# Derivative is zero away from discontinuities
evaluate_potential_derivative(::SquareWell, ::Number) = 0.0

# Discontinuities at σ and λσ (return all pairwise values if matrix)
function discontinuities(p::SquareWell)
    σ, λ = p.σ, p.λ
    if σ isa AbstractArray
        return [vec(float.(σ)); vec(float.(λ .* σ))]  # concat
    else
        return [float(σ), float(λ*σ)]
    end
end


##############################
# Morse (with mixtures)      #
##############################

"""
    Morse

Morse potential:

    u(r) = ϵ * (exp(-2α(r-σ)) - 2 exp(-α(r-σ)))

Constructors:
- `Morse(ϵ::Number, σ::Number, α::Number)` — single-component
- `Morse(ϵ::AbstractVector, σ::AbstractVector, α::AbstractVector)` — mixture with
    σ_ij = (σ_i + σ_j)/2, ϵ_ij = √(ϵ_i ϵ_j), α_ij = (α_i + α_j)/2
- `Morse(ϵij::AbstractMatrix, σij::AbstractMatrix, αij::AbstractMatrix)` — explicit pair tables

Example:
```julia
potential = Morse(1.0, 1.0, 2.0)
```

Example (mixture with mixing rules):
```julia
eps = [1.0, 0.8]
sig = [1.0, 1.2]
alp = [2.0, 1.5]
potential = Morse(eps, sig, alp)
```
"""
struct Morse{Tϵ,Tσ,Tα} <: Potential
    ϵ::Tϵ
    σ::Tσ
    α::Tα
end

# Single-component
Morse(ϵ::Number, σ::Number, α::Number) = Morse{typeof(ϵ),typeof(σ),typeof(α)}(ϵ, σ, α)

# Mixture from vectors (common symmetric rules)
function Morse(ϵ::AbstractVector{Te}, σ::AbstractVector{Ts}, α::AbstractVector{Ta}) where {Te<:Number,Ts<:Number,Ta<:Number}
    @assert length(ϵ)==length(σ)==length(α)
    N  = length(ϵ)
    ϵs = SVector{N,Te}(ϵ)
    σs = SVector{N,Ts}(σ)
    αs = SVector{N,Ta}(α)
    ϵij = sqrt.(ϵs * ϵs')
    σij = (σs .+ σs') / 2
    αij = (αs .+ αs') / 2
    return Morse{typeof(ϵij),typeof(σij),typeof(αij)}(ϵij, σij, αij)
end

# Mixture from explicit pair tables
function Morse(ϵij::AbstractMatrix{Te}, σij::AbstractMatrix{Ts}, αij::AbstractMatrix{Ta}) where {Te<:Number,Ts<:Number,Ta<:Number}
    @assert size(ϵij,1)==size(ϵij,2)==size(σij,1)==size(σij,2)==size(αij,1)==size(αij,2)
    N  = size(ϵij,1)
    ϵM = SMatrix{N,N,Te}(ϵij)
    σM = SMatrix{N,N,Ts}(σij)
    αM = SMatrix{N,N,Ta}(αij)
    return Morse{typeof(ϵM),typeof(σM),typeof(αM)}(ϵM, σM, αM)
end

function evaluate_potential(p::Morse, r::Number)
    ϵ, σ, α = p.ϵ, p.σ, p.α
    x  = r .- σ
    e1 = exp.(-α .* x)
    return ϵ .* (e1.^2 .- 2 .* e1)
end

# du/dr = ϵ * (-2α e^{-2αx} + 2α e^{-αx})
function evaluate_potential_derivative(p::Morse, r::Number)
    ϵ, σ, α = p.ϵ, p.σ, p.α
    x  = r .- σ
    e1 = exp.(-α .* x)
    return ϵ .* (-2 .* α .* (e1.^2) .+ 2 .* α .* e1)
end

discontinuities(::Morse) = Float64[]


#######################
# Tabulated Potential #
#######################

"""
    TabulatedPotential

Piecewise-linear potential defined on a sorted grid `r_grid` with values `u_grid`.

- Linear interpolation within tabulated range.
- Extrapolation behavior controlled by `extrapolation`:
   - `:error` (default) — throw if r is outside [r_min, r_max]
   - `:flat` — clamp to end values
   - `:linear` — extend linearly using end slope

Example:
```julia
r = range(0.8, 6.0; length=500) |> collect
u = @. 4.0*((1.0/r)^12 - (1.0/r)^6)  # LJ shape as a table
pot = TabulatedPotential(r, u, :flat)
```
"""
struct TabulatedPotential{TR,TU} <: Potential
    r::TR
    u::TU
    extrapolation::Symbol
    function TabulatedPotential(r::AbstractVector{Tr}, u::AbstractVector{Tu}, extrapolation::Symbol=:error) where {Tr<:Real,Tu<:Real}
        @assert length(r) == length(u) "r and u must be the same length"
        @assert issorted(r) "r grid must be sorted ascending"
        return new{typeof(r), typeof(u)}(r, u, extrapolation)
    end
end

@inline function _bracket(r::AbstractVector{<:Real}, x::Real)
    return searchsortedlast(r, x)
end

@inline function _interp_linear(x, x1, y1, x2, y2)
    t = (x - x1) / (x2 - x1)
    return (1 - t) * y1 + t * y2
end

@inline function _slope(x1, y1, x2, y2)
    return (y2 - y1) / (x2 - x1)
end

function evaluate_potential(p::TabulatedPotential, r::Number)
    rg, ug, mode = p.r, p.u, p.extrapolation
    rmin, rmax = rg[1], rg[end]

    if r < rmin || r > rmax
        if mode === :error
            error("TabulatedPotential: r=$(r) outside [$(rmin), $(rmax)]")
        elseif mode === :flat
            return r < rmin ? ug[1] : ug[end]
        elseif mode === :linear
            if r < rmin
                m = _slope(rg[1], ug[1], rg[2], ug[2])
                return ug[1] + m*(r - rg[1])
            else
                m = _slope(rg[end-1], ug[end-1], rg[end], ug[end])
                return ug[end] + m*(r - rg[end])
            end
        else
            error("Unknown extrapolation mode $(mode)")
        end
    end

    i = _bracket(rg, r)
    if i == length(rg)
        return ug[end]
    elseif rg[i] == r
        return ug[i]
    else
        return _interp_linear(r, rg[i], ug[i], rg[i+1], ug[i+1])
    end
end

function evaluate_potential_derivative(p::TabulatedPotential, r::Number)
    rg, ug, mode = p.r, p.u, p.extrapolation
    rmin, rmax = rg[1], rg[end]

    # Out-of-range behavior
    if r < rmin || r > rmax
        if mode === :error
            error("TabulatedPotential derivative: r=$(r) outside [$(rmin), $(rmax)]")
        elseif mode === :flat
            return zero(promote_type(eltype(rg), eltype(ug)))
        elseif mode === :linear
            if r < rmin
                return _slope(rg[1], ug[1], rg[2], ug[2])
            else
                return _slope(rg[end-1], ug[end-1], rg[end], ug[end])
            end
        else
            error("Unknown extrapolation mode $(mode)")
        end
    end

    # In-range: piecewise-constant slope; at nodes use right slope (left at last node)
    i = _bracket(rg, r)  # largest index with rg[i] ≤ r
    if i == length(rg)             # at the last grid point
        return _slope(rg[end-1], ug[end-1], rg[end], ug[end])
    elseif rg[i] == r              # exactly on a node (not the last): use right slope
        return _slope(rg[i], ug[i], rg[i+1], ug[i+1])
    else                           # inside cell i → i+1
        return _slope(rg[i], ug[i], rg[i+1], ug[i+1])   # <-- fixed: ug[i] on the left
    end
end


discontinuities(::TabulatedPotential) = Float64[]

"""
    Closure

Abstract closure type
"""
abstract type Closure end

"""
    PercusYevick

Implements the Percus-Yevick closure c(r) = f(r)*(1+γ(r)), or equivalently b(r) = ln(1 + γ(r)) - γ(r).

Example:
`closure = PercusYevick()`
"""
struct PercusYevick <: Closure end

"""
    HypernettedChain

Implements the Hypernetted Chain closure c(r) = (f(r)+1)*exp(γ(r)) - γ(r) - 1, or equivalently b(r) = 0.

Example:
`closure = HypernettedChain()`
"""
struct HypernettedChain <: Closure end

"""
    MeanSphericalApproximation

Implements the MSA closure c(r) = -βu(r), or equivalently b(r) = ln(γ(r) - βu(r) + 1) - γ(r) +  βu(r).

Example:
`closure = MeanSphericalApproximation()`
"""
struct MeanSphericalApproximation <: Closure end


function bridge_function(::HypernettedChain, γ, _, _)
    zerounit = zero.(γ)
    B = zerounit
    return B
end

function bridge_function(::MeanSphericalApproximation, γ, u_long_range, _) 
    oneunit = one.(γ)
    s = @. γ - u_long_range # temp s = γ*
    B = @. log(oneunit + s) - s
    return B
end

function c_closure_from_γ(closure, r, mayer_f, γ, u_long_range)
    B = bridge_function(closure, γ, u_long_range, r)
    myone = one.(B)
    c = @. -myone - γ + (mayer_f + myone)*exp(γ)*real(exp(B))
    return c
end

function cmulr_closure_from_Γmulr(closure::Closure, r, mayer_f, Γmulr, u_long_range)
    γ = Γmulr/r
    return r*c_closure_from_γ(closure, r, mayer_f, γ, u_long_range) 
end


function cmulr_closure_from_Γmulr(::PercusYevick, r, mayer_f::T, Γmulr::T, u_long_range::T) where T
    return  @. mayer_f*(r + Γmulr)
end



"""
    Closure

Abstract closure type
"""
abstract type Closure end

"""
    PercusYevick

Implements the Percus-Yevick closure \$c(r) = f(r)(1+\\gamma(r))\$, or equivalently \$b(r) = \\ln(1 + \\gamma(r)) - γ(r)\$.

Example:
```julia
closure = PercusYevick()
```
"""
struct PercusYevick <: Closure end

"""
    HypernettedChain

Implements the Hypernetted Chain closure \$c(r) = (f(r)+1)\\exp(\\gamma(r)) - \\gamma(r) - 1\$, or equivalently \$b(r) = 0\$.

Example:
```julia
closure = HypernettedChain()
```
"""
struct HypernettedChain <: Closure end

"""
    MeanSphericalApproximation

Implements the MSA closure \$c(r) = -\\beta u(r)\$, or equivalently \$b(r) = \\ln(\\gamma(r) - \\beta u(r) + 1) - γ(r) +  \\beta u(r)\$.

Example:
```julia
closure = MeanSphericalApproximation()
```
"""
struct MeanSphericalApproximation <: Closure end


function bridge_function(::HypernettedChain, _, _, γ, _)
    zerounit = zero.(γ)
    B = zerounit
    return B
end

function bridge_function(::MeanSphericalApproximation, _, _, γ, u_long_range) 
    oneunit = one.(γ)
    s = @. γ - u_long_range # temp s = γ*
    B = @. log(oneunit + s) - s
    return B
end

function closure_c_from_gamma(closure, r, mayer_f, γ, u_long_range)
    B = bridge_function(closure, r, γ, u_long_range, r)
    myone = one.(B)
    c = @. -myone - γ + (mayer_f + myone)*exp(γ)*real(exp(B))
    return c
end

function closure_cmulr_from_gammamulr(closure::Closure, r, mayer_f, Γmulr, u_long_range)
    γ = Γmulr/r
    return r*closure_c_from_gamma(closure, r, mayer_f, γ, u_long_range) 
end


function closure_cmulr_from_gammamulr(::PercusYevick, r, mayer_f::T, Γmulr::T, _::T) where T
    return  @. mayer_f*(r + Γmulr)
end



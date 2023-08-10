abstract type Closure end

struct PercusYevick <: Closure end
struct HypernettedChain <: Closure end
struct MeanSpherical <: Closure end



# function bridge_function(::PercusYevick, γ, _, _)
#     oneunit = one.(γ)
#     B = @. log(Complex(oneunit + γ)) - γ
#     return B
# end

function bridge_function(::HypernettedChain, γ, _, _)
    zerounit = zero.(γ)
    B = zerounit
    return B
end

function bridge_function(::MeanSpherical, γ, u_long_range, _) 
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



"""
    Closure

Abstract closure type
"""
abstract type Closure end

function closure_c_from_gamma(closure::Closure, r, mayer_f, γ, βuLR)
    B = bridge_function(closure, r, mayer_f, γ .- βuLR)
    myone = one.(B)
    c = @. -myone - γ + (mayer_f + myone)*exp(γ)*real(exp(B))
    return c
end

function closure_cmulr_from_gammamulr(closure::Closure, r, mayer_f, Γmulr, βuLR)
    γ = Γmulr/r
    return r*closure_c_from_gamma(closure, r, mayer_f, γ, βuLR) 
end

"""
    PercusYevick

Implements the Percus-Yevick closure \$c(r) = f(r)(1+\\gamma(r))\$, or equivalently \$b(r) = \\ln(1 + \\gamma(r)) - γ(r)\$.

Example:
```julia
closure = PercusYevick()
```
"""
struct PercusYevick <: Closure end

function closure_cmulr_from_gammamulr(::PercusYevick, r::Number, mayer_f::T, Γmulr::T) where T
    return  @. mayer_f*(r + Γmulr)
end

function bridge_function(::PercusYevick, _, _, γ)
    B  = @. log1p(γ) - γ
    return B
end

"""
    HypernettedChain

Implements the Hypernetted Chain closure \$c(r) = (f(r)+1)\\exp(\\gamma(r)) - \\gamma(r) - 1\$, or equivalently \$b(r) = 0\$.

Example:
```julia
closure = HypernettedChain()
```
"""
struct HypernettedChain <: Closure end

function bridge_function(::HypernettedChain, _, _, γ)
    zerounit = zero.(γ)
    B = zerounit
    return B
end

"""
    MeanSpherical

Implements the MSA closure \$c(r) = -\\beta u(r)\$, or equivalently \$b(r) = \\ln(\\gamma(r) - \\beta u(r) + 1) - γ(r) +  \\beta u(r)\$.

Example:
```julia
closure = MeanSpherical()
```
"""
struct MeanSpherical <: Closure end

function bridge_function(::MeanSpherical, _, mayer_f, γ) 
    oneunit = one.(γ)
    βu = @. log(mayer_f+oneunit)
    s = @. γ - βu 
    B = @. log1p(s) - s
    return B
end

"""
    Verlet <: Closure

Implements the Verlet Closure \$b(r) = -\\frac{A\\gamma^2(r)/2}{1+B \\gamma(r)/2}\$, where by default \$A=1\$ and \$B=8/5\$. These values are tuned by the virial coefficients of the 3d hard sphere liquid.

Example:
```julia
closure = Verlet()
closure = Verlet(A=3.0, B=4.0)
```
References:

Verlet, Loup. "Integral equations for classical fluids: I. The hard sphere case." Molecular Physics 41.1 (1980): 183-190.

"""
struct Verlet{T} <: Closure 
    A::T
    B::T
end

function Verlet(; A=1.0, B=8.0/5.0)
    return Verlet(A, B)
end

function bridge_function(closure::Verlet, _, _, γ)
    A = closure.A; B = closure.B
    oneunit = one.(γ)
    return @. -(A*γ^2/2)/(oneunit + B*γ/2)
end


"""
    MartynovSarkisov <: Closure

Implements the Martynov-Sarkisov Closure \$b(r) = -\\sqrt{1+2\\gamma}-1-\\gamma, which is constructed for the 3d hard-sphere liquid.

Example:
```julia
closure = MartynovSarkisov()
```
References:

Martynov, G. A., and G. N. Sarkisov. "Exact equations and the theory of liquids. V." Molecular Physics 49.6 (1983): 1495-1504.

"""
struct MartynovSarkisov <: Closure end

function bridge_function(::MartynovSarkisov, _, _, γ)
    oneunit = one.(γ)
    return @. sqrt(oneunit+2γ) - oneunit - γ
end



"""
    SoftCoreMeanSpherical <: Closure

Implements the soft core mean spherical closure \$b(r) = \\ln(\\gamma^*(r) + 1) - \\gamma^*(r)\$. Here \$\\gamma^* = \\gamma - u_{LR}\$ , in which \$u_{LR}\$ is the long range tail of the potential.

Example:
```julia
closure = SoftCoreMeanSpherical()
```

References:


"""
struct SoftCoreMeanSpherical <: Closure end

function bridge_function(::SoftCoreMeanSpherical, _, _, γ)
    γstar = γ
    return @. -γstar + log1p(γstar)
end

"""
    RogersYoung <: Closure

Implements the Rogers-Young closure \$b(r) = \\ln(\\frac{\\exp(f(r)\\gamma(r)) - 1}{f(r)} + 1) - γ(r) \$. 
Here \$f(r)=1-\\exp(-\\alpha r)\$, in which \$\\alpha\$ is a free parameter, 
that may be chosen such that thermodynamic consistency is achieved.
Example:
```julia
closure = RogersYoung(0.5)
```

References:


"""
struct RogersYoung{T} <: Closure 
    α::T
end

function bridge_function(closure::RogersYoung, r, _, γ)
    oneunit = one.(γ)
    α = closure.α
    @assert α > 0 
    f = @. 1.0 - exp(-α*r)
    b = @. -γ + log1p((exp(f*γ)-oneunit)/f)
    return b
end

"""
ZerahHansen <: Closure

Implements the Zerah-Hansen (HMSA) (HNC-SMSA) closure \$b(r) = \\ln(\\frac{\\exp(f(r)\\gamma^*(r)) - 1}{f(r)} + 1) - γ^*(r) \$.
Here  \$\\gamma^* = \\gamma - u_{LR}\$ , in which \$ u_{LR}\$ is the long range tail of the potential, 
and \$f(r)=1-\\exp(-\\alpha r)\$, in which \$\\alpha\$ is a free parameter, 
that may be chosen such that thermodynamic consistency is achieved.

Example:
```julia
closure = ZerahHansen(0.5)
```

References:

Zerah, Gilles, and Jean‐Pierre Hansen. "Self‐consistent integral equations for fluid pair distribution functions: Another attempt." The Journal of chemical physics 84.4 (1986): 2336-2343.

"""
struct ZerahHansen{T<:Number} <: Closure 
    α::T
end

function bridge_function(closure::ZerahHansen, r, _, γstar)
    α = closure.α
    f = 1.0 - exp(-α*r)
    return @. -γstar + log1p((exp(f*γstar)-1)/f)
end

"""
    DuhHaymet <: Closure

Implements the Duh-Haymet closure \$b(r) = \\frac{-\\gamma^*(r)^2}{2\\left[1+\\left(\\frac{5\\gamma^*(r)+11}{7\\gamma^*(r)+9}\\right)\\gamma^*(r)\\right]} \$ for \$\\gamma^*(r) >0\$, and \$b(r)=-\\gamma^*(r)^2/2\$ otherwise.
Here  \$\\gamma^* = \\gamma - u_{LR}\$ , in which \$ u_{LR}\$ is the long range tail of the potential, 

Example:
```julia
closure = DuhHaymet()
```

References:


"""
struct DuhHaymet <: Closure end

function bridge_function(::DuhHaymet, _, _, γstar)
    oneunit = one.(γ)
    b = @. ifelse(γstar>0,  (-γstar^2)  / (2 * (oneunit + (5γstar+11oneunit)/(7γstar+9oneunit) * γstar)), -γstar^2/2)
    return b
end

"""
    Lee <: Closure

Implements the Lee closure \$b(r) = -\\frac{\\zeta\\gamma^*(r)^2}{2} \\left( 1- \\frac{\\phi \\alpha \\gamma^*(r)}{1 + \\alpha\\gamma^*(r)} \\right) \$.
Here  \$\\gamma^* = \\gamma + ρ f(r)/2\$ , in which \$ f(r)\$ is the Mayer-f function and \$\\rho\$ the density. Additionally, \$\\zeta\$,\$\\phi\$, and, \$\\alpha\$ are free parameters that 
can be determined with thermodynamic consistency or zero-separation theorems.

Example:
```julia
closure = Lee(1.073, 1.816, 1.0, 0.4) # ζ, ϕ, α, ρ
```

References:

Lee, Lloyd L. "An accurate integral equation theory for hard spheres: Role of the zero‐separation theorems in the closure relation." The Journal of chemical physics 103.21 (1995): 9388-9396.

"""
struct Lee{T} <: Closure 
    ζ::T
    ϕ::T
    α::T
end

function bridge_function(closure::Lee, _, mayerf, γ)
    oneunit = one.(γ)
    ρ = closure.ρ
    γstar = @. γ - ρ*mayerf/2
    ζ = closure.ζ
    α = closure.α
    ϕ = closure.ϕ
    return @. -(ζ*γstar^2/2)*(oneunit - (ϕ*α*γstar)/(oneunit + α*γstar))
end

"""
    ChoudhuryGhosh <: Closure

Implements the Choudhury-Ghosh closure \$b(r) = -\\frac{-\\gamma^*(r)^2}{2(1+\\alpha \\gamma^*(r))} \$ for \$\\gamma^*(r) >0\$, and \$b(r)=-\\gamma^*(r)^2/2\$ otherwise.
Here  \$\\gamma^* = \\gamma - u_{LR}\$ , in which \$ u_{LR}\$ is the long range tail of the potential, and \$\\alpha\$ is a free parameter that 
is determined by an empirical relation.

Example:
```julia
α(ρ) = 1.01752 - 0.275ρ # see the reference for this empirical relation
closure = ChoudhuryGhosh(α(0.4))
```

References:

Choudhury, Niharendu, and Swapan K. Ghosh. "Integral equation theory of Lennard-Jones fluids: A modified Verlet bridge function approach." The Journal of chemical physics 116.19 (2002): 8517-8522.
"""
struct ChoudhuryGhosh{T} <: Closure 
    α::T
end

function bridge_function(closure::ChoudhuryGhosh, _, _, γstar)
    oneunit = one.(γ)
    α = closure.α
    return @. ifelse(γstar>0,  (-γstar^2)  / (2 * (oneunit + α * γstar)), -γstar^2/2)
end

"""
    BallonePastoreGalliGazzillo <: Closure

Implements the Ballone-Pastore-Galli-Gazzillo closure \$b(r) = (1 + s \\gamma(r))^{1/s} - \\gamma(r) - 1 \$.
Here, \$s\$ is a free parameter that can be determined with thermodynamic consistency.


Example:
```julia
closure = BallonePastoreGalliGazzillo(1.5)
```

References:

Ballone, P., et al. "Additive and non-additive hard sphere mixtures: Monte Carlo simulation and integral equation results." Molecular Physics 59.2 (1986): 275-290.
"""
struct BallonePastoreGalliGazzillo{T} <: Closure 
    s::T
end

function bridge_function(closure::BallonePastoreGalliGazzillo, _, _, γ)
    oneunit = one.(γ)
    s = closure.s
    return @. (oneunit + s*γ)^(1/s) - γ - oneunit
end

"""
    VompeMartynov <: Closure

Implements the Vompe-Martynov closure \$b(r) = \\sqrt{1+2\\gamma^*(r) } - \\gamma^*(r) - 1 \$.
Here  \$\\gamma^* = \\gamma - u_{LR}\$ , in which \$ u_{LR}\$ is the long range tail of the potential, and \$\\alpha\$ is a free parameter that can be determined with thermodynamic consistency.

Example:
```julia
closure = VompeMartynov()
```

References:

Vompe, A. G., and G. A. Martynov. "The bridge function expansion and the self‐consistency problem of the Ornstein–Zernike equation solution." The Journal of chemical physics 100.7 (1994): 5249-5258.
"""
struct VompeMartynov <: Closure end

function bridge_function(::VompeMartynov, _, _, γstar)
    oneunit = one.(γstar)
    return @. sqrt(oneunit+2γstar) - oneunit - γstar
end

"""
    CharpentierJackse <: Closure

Implements the Charpentier-Jackse closure \$b(r) = \\frac{1}{2\\alpha}\\left(\\sqrt{1+4\\alpha\\gamma^*(r) } - 2\\alpha\\gamma^*(r) - 1\\right) \$.
Here  \$\\gamma^* = \\gamma - u_{LR}\$ , in which \$ u_{LR}\$ is the long range tail of the potential, and \$\\alpha\$ is a free parameter that  can be determined with thermodynamic consistency.

Example:
```julia
closure = CharpentierJackse(0.5)
```

References:

Charpentier, I., and N. Jakse. "Exact numerical derivatives of the pair-correlation function of simple liquids using the tangent linear method." The Journal of Chemical Physics 114.5 (2001): 2284-2292.
"""
struct CharpentierJackse{T} <: Closure 
    α::T
end

function bridge_function(closure::CharpentierJackse, _, _, γstar)
    oneunit = one.(γ)
    α = closure.α
    return @. (sqrt(oneunit+4α*γstar) - oneunit - 2α*γstar)/(2α)
end

"""
    BomontBretonnet <: Closure

Implements the Bomont-Bretonnet closure \$b(r) = \\sqrt{1+2\\gamma^*(r) + f \\gamma^*(r)^2} - \\gamma^*(r) - 1  \$.
Here  \$\\gamma^* = \\gamma - u_{LR}\$ , in which \$ u_{LR}\$ is the long range tail of the potential, and \$f\$ is a free parameter that can be determined with thermodynamic consistency.

Example:
```julia
closure = BomontBretonnet(0.5)
```

References:

Bomont, J. M., and J. L. Bretonnet. "A self-consistent integral equation: Bridge function and thermodynamic properties for the Lennard-Jones fluid." The Journal of chemical physics 119.4 (2003): 2188-2191.
"""
struct BomontBretonnet{T} <: Closure 
    f::T
end

function bridge_function(closure::BomontBretonnet, _, _, γstar)
    oneunit = one.(γ)
    f = closure.f
    return @. sqrt(oneunit + 2γstar + f*γstar^2) - oneunit - γstar
end


"""
    Khanpour <: Closure

Implements the Khanpour closure \$b(r) = \\frac{1}{\\alpha}\\ln(1+\\alpha\\gamma(r)) - \\gamma  \$.
Here \$\\alpha\$ is a free parameter that can be determined with thermodynamic consistency.

Example:
```julia
closure = Khanpour(0.5)
```

References:

Khanpour, Mehrdad. "A unified derivation of Percus–Yevick and hyper-netted chain integral equations in liquid state theory." Molecular Physics 120.5 (2022): e2001065.
"""
struct Khanpour{T} <: Closure 
    α::T
end

function bridge_function(closure::Khanpour, _, _, γ)
    α = closure.α
    return @. log1p(α*γ)/α - γ
end

"""
    ReferenceHypernettedChain <: Closure

Implements the Reference Hypernetted Chain closure \$b(r) = b_{HS}(r) \$. Here \$b_{HS}(r)=\\left((a_1+a_2x)(x-a_3)(x-a_4)/(a_3 a_4)\\right)^2\$ for \$x<a_4\$ and \$b_{HS}(r)=\\left(A_1 \\exp(-a_5(x-a_4))\\sin(A_2(x-a_4))/r\\right)^2\$ is the hard sphere bridge function found in Malijevský & Labík.
The parameters are defined as

\$x = r-1\$

\$A_1 = (a_1+a_2 a_4)(a_4-a_3)(a_4+1)/(A_2 a_3 a_4)\$

\$A_2 =  \\pi / (a_6 - a_4 - 1)\$

\$a_1 = \\eta (1.55707 - 1.85633\\eta) / (1-\\eta)^2\$

\$a_2 = \\eta (1.28127 - 1.82134\\eta) / (1-\\eta)\$

\$a_3 =  (0.74480 - 0.93453\\eta)\$

\$a_4 = (1.17102 - 0.68230\\eta)\$

\$a_5 = 0.15975/\\eta^2\$

\$a_6 = (2.69757 - 0.86987\\eta)\$

and \$\\eta\$ is the volume fraction of the hard sphere reference system. This closure only works for single component systems in three dimensions.

Example:
```julia
closure = ReferenceHypernettedChain(0.4)
```

References:

Lado, F. "Perturbation correction for the free energy and structure of simple fluids." Physical Review A 8.5 (1973): 2548.

Malijevský, Anatol, and Stanislav Labík. "The bridge function for hard spheres." Molecular Physics 60.3 (1987): 663-669.
"""
struct ReferenceHypernettedChain{T} <: Closure 
    η::T
end

function bridge_function(closure::ReferenceHypernettedChain, r, _, _)
    oneunit = one(r)
    η = closure.η
    x = r-oneunit
    a_1 = η * (1.55707 - 1.85633η) / (1-η)^2
    a_2 = η * (1.28127 - 1.82134η) / (1-η)
    a_3 = (0.74480 - 0.93453η) 
    a_4 = (1.17102 - 0.68230η) 
    a_5 = 0.15975/η^3
    a_6 = (2.69757 - 0.86987η)
    A_2 =  π / (a_6 - a_4 - 1)
    A_1 = (a_1+a_2*a_4)*(a_4-a_3)*(a_4+1)/(A_2*a_3*a_4)
    b = -ifelse(x<a_4, 
        ((a_1+a_2*x)*(x-a_3)*(x-a_4)/(a_3* a_4))^2,
    	(A_1*exp(-a_5*(x-a_4))*sin(A_2*(x-a_4))/r)^2
    )
    return b
end


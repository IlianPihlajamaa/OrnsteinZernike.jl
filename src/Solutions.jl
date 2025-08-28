
"""
    OZSolution

Holds the solution of an Ornstein Zernike problem. 

Fields:
- r: vector of distances
- k: vector of wave numbers
- gr: radial distribution function    
- Sk: static structure factor
- ck: direct correlation function in k space
- cr: direct correlation function in real space
- gamma_r: indirect correlation function in real space
- gamma_k: indirect correlation function in k space

if the system was a single-component system, `gr`, `Sk`, etc, are vectors. 
If instead the system was a multicomponent one, they are three dimensional vectors, 
where the first dimension contains the values along r, and the second and third dimension
contain the data for the species.
"""
struct OZSolution{T1,T2}
    r::T1
    k::T1
    gr::T2
    Sk::T2
    cr::T2
    ck::T2
    gamma_r::T2
    gamma_k::T2
end

function convert_vecofmat_to_3darr(a::Vector{T}) where {T<:AbstractMatrix}
    elT = eltype(eltype(a))
    Ns1, Ns2 = size(a[1])
    Nr = length(a)
    a = reinterpret(reshape, elT, a)
    a = reshape(a, Ns1, Ns2, Nr)
    a = permutedims(a, (3, 1, 2))
    return a
end

convert_vecofmat_to_3darr(a::Array{T,3}) where {T} = a


function find_g_from_c_and_γ(c::Vector{T}, γ::Vector{T}) where {T}
    g = copy(c)
    O = T <: SMatrix ? ones(T) : one(T)
    for i in eachindex(c, γ)
        g[i] = c[i] + γ[i] + O
    end
    return g
end

function find_S_from_ĉ_and_ρ(ĉ::Vector{T}, ρ) where {T<:Number}
    return @. one(T) + ρ * ĉ / (one(T) - ρ * ĉ)
end

function find_S_from_ĉ_and_ρ(ĉ::Vector{T}, ρ) where {T<:AbstractMatrix}
    S = copy(ĉ)
    ρtot = sum(ρ.diag)
    x = ρ.diag / ρtot
    for i in eachindex(ĉ)
        Si = inv(one(T) ./ x .- sum(ρ) * ĉ[i])
        S[i] = Si
    end
    return S
end


function construct_solution(r::T1, k::T1, cr::T, ck::T, gamma_r::T, gamma_k::T, ρ) where {T1,T<:Vector{<:AbstractMatrix}}
    gr = find_g_from_c_and_γ(cr, gamma_r)
    Sk = find_S_from_ĉ_and_ρ(ck, ρ)
    gr = convert_vecofmat_to_3darr(gr)
    Sk = convert_vecofmat_to_3darr(Sk)
    ck = convert_vecofmat_to_3darr(ck)
    cr = convert_vecofmat_to_3darr(cr)
    gamma_r = convert_vecofmat_to_3darr(gamma_r)
    gamma_k = convert_vecofmat_to_3darr(gamma_k)
    OZSolution(r, k, gr, Sk, cr, ck, gamma_r, gamma_k)
end

function construct_solution(r::T1, k::T1, cr::T, ck::T, gamma_r::T, gamma_k::T, ρ) where {T1,T<:Vector{<:Number}}
    gr = find_g_from_c_and_γ(cr, gamma_r)
    Sk = find_S_from_ĉ_and_ρ(ck, ρ)
    OZSolution(r, k, gr, Sk, cr, ck, gamma_r, gamma_k)
end


function Base.show(io::IO, ::MIME"text/plain", p::OZSolution{T1,T2}) where {T1,T2<:AbstractVector}
    println(io, "$(typeof(p)):")
    print(io, " r = ")
    show(io, p.r')
    print(io, "'\n k = ")
    show(io, p.k')
    print(io, "'\n cr = ")
    show(io, p.cr')
    print(io, "'\n gr = ")
    show(io, p.gr')
    print(io, "'\n ck = ")
    show(io, p.ck')
    print(io, "'\n Sk = ")
    show(io, p.Sk')
    print("' \n")
end

function Base.show(io::IO, ::MIME"text/plain", p::OZSolution{T1,T2}) where {T1,T2<:AbstractArray}
    println(io, "$(typeof(p)):")
    print(io, " r = ")
    show(io, p.r')
    print(io, "'\n k = ")
    show(io, p.k')
    print(io, "'\n cr = ")
    show(io, p.cr)
    print(io, "\n gr = ")
    show(io, p.gr)
    print(io, "\n ck = ")
    show(io, p.ck)
    print(io, "\n Sk = ")
    show(io, p.Sk)
    print("\n")
end
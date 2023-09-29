
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

if the system was a single-component system, `gr`, `Sk`, `ck` and `cr` are vectors. 
If instead the system was a multicomponent one, they are three dimensional vectors, 
where the first dimension contains the values along r, and the second and third dimension
contain the data for the species.
"""
struct OZSolution{T1, T2}
    r::T1
    k::T1
    gr::T2
    Sk::T2
    ck::T2
    cr::T2
end

function convert_vecofmat_to_3darr(a)
    elT = eltype(eltype(a))
    Ns1, Ns2 = size(a[1])
    Nr = length(a)
    a = reinterpret(reshape, elT, a)
    a = reshape(a, Ns1, Ns2, Nr)
    a = permutedims(a, (3,1,2))
    return a
end


function OZSolution(r::T1, k::T1, gr::T, Sk::T, ck::T, cr::T) where {T1,T<:Vector{<:AbstractMatrix}}
    gr = convert_vecofmat_to_3darr(gr)
    Sk = convert_vecofmat_to_3darr(Sk)
    ck = convert_vecofmat_to_3darr(ck)
    cr = convert_vecofmat_to_3darr(cr)
    OZSolution(r, k, gr, Sk, ck, cr)
end

function Base.show(io::IO, ::MIME"text/plain", p::OZSolution{T1,T2}) where {T1, T2<:AbstractVector}
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
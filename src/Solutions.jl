struct OZSolution{T1, T2}
    r::T1
    k::T1
    gr::T2
    Sk::T2
    Ck::T2
    Cr::T2
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


function OZSolution(r::T1, k::T1, gr::T, Sk::T, Ck::T, Cr::T) where {T1,T<:Vector{<:AbstractMatrix}}
    gr = convert_vecofmat_to_3darr(gr)
    Sk = convert_vecofmat_to_3darr(Sk)
    Ck = convert_vecofmat_to_3darr(Ck)
    Cr = convert_vecofmat_to_3darr(Cr)
    OZSolution(r, k, gr, Sk, Ck, Cr)
end

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

function solve(system::SimpleUnchargedSystem, closure::Closure, method::NgIteration; gamma_0=nothing, init=nothing)
    if !isnothing(init)
        @warn "The `init` keyword argument is deprecated. Please use `gamma_0` instead."
        isnothing(gamma_0) && gamma_0 = init
    end
    N_stages = method.N_stages
    ρ = system.ρ

    cache = OZSolverCache(system, method)
    mayer_f, fourierplan, r, k, βu_long_range, βu, Γhat, C, Ĉ, Γ_new = 
        cache.mayer_f, cache.fourierplan, cache.r, cache.k, cache.βu_long_range, cache.βu, cache.Γhat, cache.C, cache.Ĉ, cache.Γ_new

    T = eltype(mayer_f)
    TT = eltype(T)
    A = zeros(TT, N_stages, N_stages)
    b = zeros(TT, N_stages)

    # if T <: AbstractMatrix, we need to store gi as a vector of matrices
    # if T <: Number, we need to store gi as a vectors of numbers 
    gn_full = [copy(mayer_f) for _ in 1:N_stages+1] # first element is g_n, second is g_{n-1} etc
    fn_full = [copy(mayer_f) for _ in 1:N_stages+1]
    dn_full = [copy(mayer_f) for _ in 1:N_stages+1]
    d0n_full = [copy(mayer_f) for _ in 1:N_stages] # first element is d01 second is d02 etc

    # flattened arrays!!
    # if T <: Number, flatten_array is identity
    gn_flat = [flatten_array(gn_full[i]) for i in eachindex(gn_full)]
    fn_flat = [flatten_array(fn_full[i]) for i in eachindex(fn_full)]
    dn_flat = [flatten_array(dn_full[i]) for i in eachindex(dn_full)]
    d0n_flat = [flatten_array(d0n_full[i]) for i in eachindex(d0n_full)]
    Γ_new_flat = flatten_array(Γ_new)

    if !(isnothing(gamma_0))
        for i = eachindex(fn_full[end])
            fn_full[end][i] = gamma_0[i] * r[i]
        end
    else
        fill!(fn_full[end], zero(T))
    end

    max_iterations = method.max_iterations
    tolerance = method.tolerance

    err = tolerance * 2
    iteration = 0
    # first bootstrapping steps 
    for stage = reverse(1:N_stages)
        C .= closure_cmulr_from_gammamulr.((closure,), r, mayer_f, fn_full[stage+1], βu_long_range)
        oz_iteration_step!(C, Ĉ, Γhat, gn_full[stage+1], ρ, fourierplan)
        fn_flat[stage] .= gn_flat[stage+1]
    end

    while err > tolerance
        if iteration > max_iterations
            error("Recursive iteration did not converge within $iteration steps. Current error = $err.")
        end
        C .= closure_cmulr_from_gammamulr.((closure,), r, mayer_f, fn_full[1], βu_long_range)
        oz_iteration_step!(C, Ĉ, Γhat, Γ_new, ρ, fourierplan)
        gn_flat[1] .= Γ_new_flat
        err = compute_error(gn_flat[1], fn_flat[1])
        if method.verbose && iteration % 10 == 0
            println("After iteration $iteration, the error is $(round(err, digits=ceil(Int, 1-log10(tolerance)))).")
        end
        for stage = 1:N_stages+1
            @. dn_flat[stage] = gn_flat[stage] - fn_flat[stage]
        end
        for stage = 1:N_stages
            @. d0n_flat[stage] = dn_flat[1] - dn_flat[stage+1]
        end

        coefficients = find_Ng_coefficients(A, b, N_stages, d0n_flat, dn_flat, r)
        update_Γ_new_Ng!(Γ_new_flat, coefficients, gn_flat)

        for stage = reverse(1:N_stages)
            gn_flat[stage+1] .= gn_flat[stage]
            fn_flat[stage+1] .= fn_flat[stage]
        end
        fn_flat[1] .= Γ_new_flat
        iteration += 1
    end

    if method.verbose
        print("Converged after $iteration iterations, ")
        println("the error is $(round(err, digits=ceil(Int, 1-log10(err)))).")
    end
    c = C ./ r
    ĉ = Ĉ ./ k
    γ = Γ_new ./ r
    γ̂ = Γhat ./ k
    return construct_solution(r, k, c, ĉ, γ, γ̂, ρ)
end


flatten_array(a::Vector{T}) where T<:Number = a
function flatten_array(a::Vector{T}) where T<:SMatrix
    L1, L2 = size(a[1])
    @assert L1 == L2
    N = length(a)
    return reshape(reinterpret(reshape, eltype(T), a), L1 * L2 * N)
end

unflatten_array(a::Vector{S}, args...) where {S<:Number} = a
function unflatten_array(a::Base.ReshapedArray{T, 1, R}, args...) where {T<:Number, R<:Base.ReinterpretArray}
    return a.parent.parent
end

function update_Γ_new_Ng!(Γ_new::AbstractVector{T}, coeffs, gn::AbstractVector{<:AbstractVector{T}}) where {T<:Number}
    N_stages = length(gn) - 1
    Γ_new .= (one(T) - sum(coeffs)) .* gn[1]
    for stage = 1:N_stages
        Γ_new .+= coeffs[stage] .* gn[stage+1]
    end
end

function update_Γ_new_Ng!(Γ_new::Vector{T}, coeffs, gn::Vector{Vector{T}}) where {T<:AbstractMatrix}
    #coeffs is vector of matrix
    N_stages = length(gn) - 1
    for i in eachindex(Γ_new)
        Γ_new[i] = (ones(T) .- sum(coeffs)) .* gn[1][i]
        for stage = 1:N_stages
            Γ_new[i] += coeffs[stage] .* gn[stage+1][i]
        end
    end
end


function find_Ng_coefficients(A::Matrix{T}, b::Vector{T}, N_stages, d0n::Vector{S}, dn::Vector{S}, r) where {T<:Number, S<:AbstractArray}
    @assert eltype(S) == T
    for stage1 = 1:N_stages
        for stage2 = stage1:N_stages
            Aij = inner(d0n[stage1], d0n[stage2], r)
            A[stage1, stage2] = Aij
            A[stage2, stage1] = Aij
        end
        b[stage1] = inner(dn[1], d0n[stage1], r)
    end
    coeffs = A \ b
    return coeffs
end

function inner(u::AbstractArray, v::AbstractArray, r) 
    @assert eltype(u) <: Number
    T = eltype(u)
    @assert length(u) == length(v)
    if length(r) == length(u) #singlecomponent
        S = zero(T)
        for i in firstindex(u):(lastindex(u)-1)
            S += u[i] * v[i] * (r[i+1] - r[i])
        end
        return S
    else # multicomponent
        S = zero(T)
        Ns2 = length(u) ÷ length(r)
        for iss = 1:Ns2
            for i in firstindex(r):(lastindex(r)-1)
                idx = (i - 1) * Ns2 + iss
                S += u[idx] * v[idx] * (r[i+1] - r[i])
            end
        end
        return S
    end
end

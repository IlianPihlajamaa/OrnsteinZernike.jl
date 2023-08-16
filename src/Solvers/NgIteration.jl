function solve(system::SimpleLiquid{dims, 1, T1, T2, P}, closure::Closure, method::NgIteration; init=nothing) where {dims, T1, T2, P}
    N_stages = method.N_stages
    ρ = system.ρ


    r = method.dr * (1:method.M) |> collect
    r = method.dr * (1:method.M) |> collect
    mayer_f = find_mayer_f_function(system, r)
    elementtype = promote_type(eltype(r), typeof(system.kBT), typeof(system.ρ), eltype(mayer_f))

    mayer_f = elementtype.(mayer_f)

    fourierplan = get_fourier_plan(system, method, mayer_f)
    r .= fourierplan.r # in the case that dims != 3, we need to use the right grid
    k = fourierplan.k
    mayer_f .= find_mayer_f_function(system, r)
    u_long_range = copy(mayer_f)*0.0

    Γhat = copy(mayer_f)
    C = copy(mayer_f)
    Ĉ = copy(mayer_f)
    Γ_new = copy(mayer_f)
    gn = [copy(mayer_f) for _ in 1:N_stages+1] # first element is g_n, second is g_{n-1} etc
    fn = [copy(mayer_f) for _ in 1:N_stages+1]
    dn = [copy(mayer_f) for _ in 1:N_stages+1]
    d0n = [copy(mayer_f) for _ in 1:N_stages] # first element is d01 second is d02 etc
    if !(isnothing(init))
        fn[end] .= init.*r
    end
    T = eltype(mayer_f)
    A = zeros(T, N_stages, N_stages)
    b = zeros(T, N_stages)

    max_iterations = method.max_iterations
    tolerance = method.tolerance
    
    err = tolerance*2
    iteration = 0

    # first bootstrapping steps 
    for stage = reverse(1:N_stages)
        C .= closure_cmulr_from_gammamulr.((closure, ), r, mayer_f, fn[stage+1], u_long_range)
        fourier!(Ĉ, C, fourierplan)
        @. Γhat = (k - Ĉ*ρ) \ (Ĉ * ρ * Ĉ)
        inverse_fourier!(gn[stage+1], Γhat, fourierplan)
        fn[stage] .= gn[stage+1]
        # @show error = compute_error(gn[stage+1], fn[stage+1])
    end

    while err > tolerance
        if iteration > max_iterations
            error("Recursive iteration did not converge within $iteration steps. Current error = $err.")
        end
        Γold = fn[1]  
        C .= closure_cmulr_from_gammamulr.((closure, ), r, mayer_f, Γold, u_long_range)
        fourier!(Ĉ, C, fourierplan)
        @. Γhat = (k - Ĉ*ρ) \ (Ĉ * ρ * Ĉ)

        inverse_fourier!(Γ_new, Γhat, fourierplan)
        gn[1] .= Γ_new
        err = compute_error(gn[1], fn[1])
        if method.verbose && iteration % 10 == 0
            println("After iteration $iteration, the error is $(round(err, digits=ceil(Int, 1-log10(tolerance)))).")
        end
        for stage = 1:N_stages+1
            @. dn[stage] = gn[stage] - fn[stage]
        end
        for stage = 1:N_stages
            @. d0n[stage] = dn[1] - dn[stage+1]
        end
        
        coefficients = find_Ng_coefficients(A, b, N_stages, d0n, dn, r)
        update_Γ_new_Ng!(Γ_new, coefficients, gn)

        for stage = reverse(1:N_stages)
            gn[stage+1] .= gn[stage]
            fn[stage+1] .= fn[stage]
        end
        fn[1] .= Γ_new
        iteration += 1

    end

    if method.verbose 
        print("Converged after $iteration iterations, ")
        println("the error is $(round(err, digits=ceil(Int, 1-log10(tolerance)))).")
    end
    c = C ./ r
    g = find_g_from_c_and_Γ(c, Γ_new, r)
    ĉ = Ĉ ./ k
    Sk =  find_S_from_ĉ_and_ρ(ĉ, ρ)
    return OZSolution(r, k, g, Sk, ĉ, c)
end 

function initialize_vector_of_vectors(::Type{T}, N, L) where T
    return [zeros(T, L) for _ = 1:N]
end

function reinterpret_vector_of_vectors(vecofvec, ::Type{T}, L1, L2) where T
    return [reinterpret(reshape, T, reshape(vecofvec[i], (L1, L2))) for i in eachindex(vecofvec)]
end

function solve(system::SimpleLiquid{dims, species, T1, T2, P}, closure::Closure, method::NgIteration; init=nothing) where {dims, species, T1, T2, P}
    N_stages = method.N_stages
    ρ = system.ρ

    r = method.dr * (1:method.M) |> collect
    mayer_f = find_mayer_f_function(system, r)
    elementtype = promote_type(typeof(eltype(r)), typeof(system.kBT), typeof(system.ρ), eltype(mayer_f))
    mayer_f = elementtype.(mayer_f)
    fourierplan = get_fourier_plan(system, method, mayer_f .* (ρ * system.kBT))
    r .= fourierplan.r # in the case that dims != 3, we need to use the right grid
    k = fourierplan.k
    mayer_f .= find_mayer_f_function(system, r) 

    Ns = length(ρ.diag)
    Nr = length(r)

    u_long_range = copy(mayer_f)*0.0
    T = eltype(mayer_f)
    TT = eltype(T)
    Γhat = copy(mayer_f)
    C = copy(mayer_f)
    Ĉ = copy(mayer_f)
    Γ_new = copy(mayer_f)

    gn = initialize_vector_of_vectors(TT, N_stages+1, Ns*Ns*Nr) # first element is g_n, second is g_{n-1} etc
    fn = initialize_vector_of_vectors(TT, N_stages+1, Ns*Ns*Nr)
    dn = initialize_vector_of_vectors(TT, N_stages+1, Ns*Ns*Nr)
    d0n = initialize_vector_of_vectors(TT, N_stages, Ns*Ns*Nr) # first element is d01 second is d02 etc
    if !(isnothing(init))
        fn[end] .= init.*r
    end
    A = zeros(TT, N_stages, N_stages)
    b = zeros(TT, N_stages)
    Γ_new_full = reshape(reinterpret(reshape, TT, Γ_new), Ns*Ns*Nr) # for going back and forth between vec{float} and vec{Smat}
    gn_red = reinterpret_vector_of_vectors(gn, T, Ns*Ns, Nr)#[reinterpret(reshape, T, reshape(gn[i], (Ns*Ns, Nr))) for i in eachindex(gn)]
    fn_red = reinterpret_vector_of_vectors(fn, T, Ns*Ns, Nr)#[reinterpret(reshape, T, reshape(fn[i], (Ns*Ns, Nr))) for i in eachindex(fn)]

    max_iterations = method.max_iterations
    tolerance = method.tolerance
    
    err = tolerance*2
    iteration = 0
    # first bootstrapping steps 
    for stage = reverse(1:N_stages)
        C .= closure_cmulr_from_gammamulr.((closure, ), r, mayer_f, fn_red[stage+1], u_long_range)
        fourier!(Ĉ, C, fourierplan)
        for ik in eachindex(Γhat, Ĉ)
            Γhat[ik] = (k[ik] * I - Ĉ[ik]*ρ) \ (Ĉ[ik] * ρ * Ĉ[ik])
        end
        inverse_fourier!(gn_red[stage+1], Γhat, fourierplan)
        fn[stage] .= gn[stage+1]
    end

    while err > tolerance
        if iteration > max_iterations
            error("Recursive iteration did not converge within $iteration steps. Current error = $err.")
        end
        C .= closure_cmulr_from_gammamulr.((closure, ), r, mayer_f, fn_red[1], u_long_range)
        fourier!(Ĉ, C, fourierplan)
        for ik in eachindex(Γhat, Ĉ)
            Γhat[ik] = (k[ik] * I - Ĉ[ik]*ρ) \ (Ĉ[ik] * ρ * Ĉ[ik])
        end
        inverse_fourier!(Γ_new, Γhat, fourierplan)
        gn[1] .= Γ_new_full
        err = compute_error(gn[1], fn[1])
        if method.verbose && iteration % 10 == 0
            println("After iteration $iteration, the error is $(round(err, digits=ceil(Int, 1-log10(tolerance)))).")
        end
        for stage = 1:N_stages+1
            @. dn[stage] = gn[stage] - fn[stage]
        end
        for stage = 1:N_stages
            @. d0n[stage] = dn[1] - dn[stage+1]
        end
        
        coefficients = find_Ng_coefficients(A, b, N_stages, d0n, dn, r)
        update_Γ_new_Ng!(Γ_new_full, coefficients, gn)

        for stage = reverse(1:N_stages)
            gn[stage+1] .= gn[stage]
            fn[stage+1] .= fn[stage]
        end
        fn[1] .= Γ_new_full
        iteration += 1

    end

    if method.verbose 
        print("Converged after $iteration iterations, ")
        println("the error is $(round(err, digits=ceil(Int, 1-log10(tolerance)))).")
    end
    c = C ./ r
    g = find_g_from_c_and_Γ(c, Γ_new, r)
    ĉ = Ĉ ./ k
    Sk =  find_S_from_ĉ_and_ρ(ĉ, ρ)
    return OZSolution(r, k, g, Sk, ĉ, c)
end 

function update_Γ_new_Ng!(Γ_new::AbstractVector{T}, coeffs, gn::AbstractVector{<:AbstractVector{T}}) where T<:Number
    N_stages = length(gn)-1
    Γ_new .= (one(T) - sum(coeffs)).*gn[1] 
    for stage = 1:N_stages
        Γ_new .+= coeffs[stage].*gn[stage+1]
    end
end

function update_Γ_new_Ng!(Γ_new::Vector{T}, coeffs, gn::Vector{Vector{T}}) where T<:AbstractMatrix
    #coeffs is vector of matrix
    N_stages = length(gn)-1
    for i in eachindex(Γ_new)
        Γ_new[i] = (ones(T) .- sum(coeffs)) .* gn[1][i] 
        for stage = 1:N_stages
            Γ_new[i] += coeffs[stage] .* gn[stage+1][i]
        end
    end
end
# function update_Γ_new_Ng!(Γ_new::Vector{T}, coeffs, gn::Vector{Vector{T}}) where T<:AbstractMatrix
#     #coeffs is vector of matrix
#     N_stages = length(gn)-1
#     for i in eachindex(Γ_new)
#         Γ_new[i] = (1.0 - sum(coeffs)) * gn[1][i] 
#         for stage = 1:N_stages
#             Γ_new[i] += coeffs[stage] * gn[stage+1][i]
#         end
#     end
# end

function find_Ng_coefficients(A::Matrix{T}, b::Vector{T}, N_stages, d0n::Vector{Vector{T}}, dn::Vector{Vector{T}}, r) where T<:Number
    for stage1= 1:N_stages
        for stage2 = stage1:N_stages
            Aij = inner(d0n[stage1],d0n[stage2],r)
            A[stage1, stage2] = Aij
            A[stage2, stage1] = Aij
        end
        b[stage1] = inner(dn[1], d0n[stage1], r)
    end
    coeffs = A\b
    return coeffs
end


# function find_Ng_coefficients(A::Matrix{T}, b::Vector{T}, N_stages, d0n::Vector{Vector{T}}, dn::Vector{Vector{T}}, r) where T<:AbstractMatrix
#     Ns = size(d0n[1][1], 1)
#     for stage1= 1:N_stages
#         for stage2 = stage1:N_stages
#             Aij = inner(d0n[stage1], d0n[stage2],r)
#             println("3")

#             A[stage1, stage2] = Aij
#             A[stage2, stage1] = Aij
#         end
#         println("4")

#         b[stage1] = inner(dn[1], d0n[stage1], r)
#     end
#     #coeffs vector of matrix, each element is a matrix of the per species coeffs of that stage
#     coeffs = zeros(Ns, Ns, N_stages)
#     for species2 = 1:Ns
#         for species1 = 1:Ns
#             coeffs12 = getindex.(A, species1, species2)\getindex.(b, species1, species2)
#             coeffs[species1, species2, :] .= coeffs12
#         end
#     end
#     coeffs = reshape(coeffs, (Ns*Ns, N_stages))
#     coeffs = reinterpret(reshape, SMatrix{Ns, Ns, eltype(T), Ns*Ns}, coeffs)
#     coeffs = Vector(coeffs)
#     return coeffs
# end


function inner(u::Vector{T},v::Vector{T}, r) where T<:Number
    @assert length(u) == length(v)
    if length(r) == length(u) #singlecomponent
        S = zero(T)
        for i in firstindex(u):(lastindex(u)-1)
            S += u[i]*v[i]*(r[i+1]-r[i])
        end
        return S
    else # multicomponent
        S = zero(T)
        Ns2 =  length(u)÷length(r)
        for iss = 1:Ns2
            for i in firstindex(r):(lastindex(r)-1)
                idx = (i-1)*Ns2 + iss
                S += u[idx]*v[idx]*(r[i+1]-r[i])
            end
        end
        return S
    end
end

# function inner(u::Vector{T},v::Vector{T}, r) where T<:AbstractMatrix
#     S = zero(T)
#     for i in firstindex(u):(lastindex(u)-1)
#         S = S + u[i] .* v[i] .* (r[i+1]-r[i])
#     end
#     return S
# end
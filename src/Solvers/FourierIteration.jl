function solve(system::SimpleUnchargedSystem, closure::Closure, method::FourierIteration; gamma_0=nothing, init=nothing)
    if !isnothing(init)
        @warn "The `init` keyword argument is deprecated. Please use `gamma_0` instead."
        if isnothing(gamma_0); gamma_0 = init; end
    end
    ρ = system.ρ

    cache = OZSolverCache(system, method)
    mayer_f, fourierplan, r, k, βu_long_range, βu, Γhat, C, Ĉ, Γ_new = 
        cache.mayer_f, cache.fourierplan, cache.r, cache.k, cache.βu_long_range, cache.βu, cache.Γhat, cache.C, cache.Ĉ, cache.Γ_new
    
    C_old = copy(mayer_f); Γ_old = copy(mayer_f)

    if isnothing(gamma_0)
        fill!(Γ_old, zero(eltype(Γ_old)))
    else
        for i = eachindex(Γ_old)
            Γ_old[i] = gamma_0[i] * r[i]
        end
    end

    max_iterations = method.max_iterations
    tolerance = method.tolerance
    mixing_parameter = method.mixing_parameter

    err = tolerance * 2
    iteration = 0

    while err > tolerance
        if iteration > max_iterations
            error("Recursive iteration did not converge within $iteration steps. Current error = $err.")
        end
        C .= closure_cmulr_from_gammamulr.((closure, ), r, mayer_f, Γ_old, βu_long_range) 
        if iteration != 0 
            @. C = mixing_parameter * C + (1.0 - mixing_parameter) * C_old 
        end 
        oz_iteration_step!(C, Ĉ, Γhat, Γ_new, ρ, fourierplan)
        err = compute_error(Γ_new, Γ_old)

        if method.verbose && iteration % 100 == 0
            println("After iteration $iteration, the error is $(round(err, digits=ceil(Int, 5-log10(tolerance)))).")
        end
        Γ_old .= Γ_new
        C_old .= C
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

function solve(system::SimpleChargedSystem, closure::Closure, method::FourierIteration; gamma_0=nothing, CoulombSplitting=NoCoulombSplitting())
    ρ = system.ρ

    # construct solver cache from base system. Does not include Coulomb parts
    cache = OZSolverCache(base_of(system), method)
    mayer_f, fourierplan, r, k, βu_LR_disp, βu, γ_SR_hat, c_SR, c_SR_hat, γ_SR_new = 
        cache.mayer_f, cache.fourierplan, cache.r, cache.k, cache.βu_long_range, cache.βu, cache.Γhat, cache.C, cache.Ĉ, cache.Γ_new
    
    # definitions:
    # βu = βu_short_range + βu_LR_disp + βu_SR_coul + βu_LR_coul (first two are WCA splitting latter are coul splitting)
    # c = c_short_range + Φ (long range part is just Φ=-βu_LR_coul)
    # γ = γ_short_range + γ_long_range 
    # h = h_short_range + q (defined through long range oz relation: q = Φ + Φ ρ q (fourier))

    
    βu_SR_coul, βu_LR_coul = coulomb_splitting(system, r, CoulombSplitting)
    βu = βu .+ βu_SR_coul .+ βu_LR_coul

    Φ = -βu_LR_coul # long range part of c 
    

    Φ_hat = fourier(Φ .* r, fourierplan) ./ k # note that the FT work on the "mulr" functions
    q_hat = similar(Φ_hat)

    for i = eachindex(q_hat)
        if k[i] == 0.0
            Φ_hat[i] = zero(eltype(Φ_hat))
            q_hat[i] = zero(eltype(q_hat))
        else
            q_hat[i] =  (I - Φ_hat[i] * ρ) \ Φ_hat[i]
        end
    end
    q = inverse_fourier(q_hat .* k, fourierplan) ./ r

    c_SR_old = copy(mayer_f); γ_SR_old = copy(mayer_f)
    c = copy(mayer_f); c_hat = copy(mayer_f)
    γ = copy(mayer_f); γ_hat = copy(mayer_f)


    # initial guess for short ranged gamma    
    if isnothing(gamma_0)
        fill!(γ_SR_old, zero(eltype(γ_SR_old)))
    else
        for i = eachindex(γ_SR_old)
            γ_SR_old[i] = gamma_0[i]
        end
    end

    max_iterations = method.max_iterations
    tolerance = method.tolerance
    mixing_parameter = method.mixing_parameter
    err = tolerance * 2
    iteration = 0

    while err > tolerance
        c_SR .= closure_c_from_gamma_short_ranged.((closure, ), r, βu, γ_SR_old, q, βu_LR_disp, βu_LR_coul) 
        if iteration != 0 
            @. c_SR = mixing_parameter * c_SR + (1.0 - mixing_parameter) * c_SR_old 
        end    
        c_SR .*= r
        fourier!(c_SR_hat, c_SR, fourierplan)
        c_SR_hat ./= k
        c_hat .= c_SR_hat .+ Φ_hat
        @. γ_hat = (I - c_hat * ρ) \ (c_hat * ρ * c_hat) 
        γ_SR_hat .= γ_hat .- (q_hat - Φ_hat) 
        γ_SR_hat .*= k
        inverse_fourier!(γ_SR_new, γ_SR_hat, fourierplan)
        γ_SR_old .*= r
        err = compute_error(γ_SR_new, γ_SR_old) # compute error on mul_r functions, consistent with noncharged case
        γ_SR_new ./= r

        if method.verbose && iteration % 100 == 0
            println("After iteration $iteration, the error is $(round(err, digits=ceil(Int, 5-log10(tolerance)))).")
        end

        γ_SR_old .= γ_SR_new
        c_SR_old .= c_SR

        iteration += 1
    end
    if method.verbose
        print("Converged after $iteration iterations, ")
        println("the error is $(round(err, digits=ceil(Int, 1-log10(err)))).")
    end
    c .= c_SR .+ Φ
    c_hat .= c_SR_hat .+ Φ_hat
    γ .= γ_SR_new .+ q .- Φ
    γ_hat .= γ_SR_hat .+ q_hat .- Φ_hat
    return construct_solution(r, k, c, c_hat, γ, γ_hat, ρ)
end



"""
βu_LR_disp is the long ranged part of the base potential only, e.g. r^-6 for LJ, not including Coulomb)
γ_SR is the short ranged part of γ, i.e. with the long coulomb tail subtracted
q is the long ranged part of h, i.e.  the long coulomb tail
the mayer function includes the full potential
"""
function closure_c_from_gamma_short_ranged(closure::Closure, r::Number, βu, γ_SR, q, βu_LR_disp, βu_long_range_coul)
    mayer_f = find_mayer_f_function.(βu)
    B = bridge_function(closure, r, mayer_f, γ_SR .- βu_LR_disp) # same as in uncharged case
    myone = one.(B)
    c = @. -myone - γ_SR - q + exp(-(βu- βu_long_range_coul) + γ_SR + q + real.(B))
    return c
end
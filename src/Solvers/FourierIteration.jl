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


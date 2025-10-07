function solve(system::SimpleUnchargedSystem, closure::Closure, method::FourierIteration; gamma_0=nothing, init=nothing)
    if !isnothing(init)
        @warn "The `init` keyword argument is deprecated. Please use `gamma_0` instead."
        if isnothing(gamma_0); gamma_0 = init; end
    end
    ρ = system.ρ

    ensure_dispersion_support(closure, system.potential)

    cache = OZSolverCache(system, method)
    mayer_f, fourierplan, r, k, βu_disp_tail, βu, Γhat, C, Ĉ, Γ_new = 
        cache.mayer_f, cache.fourierplan, cache.r, cache.k, cache.βu_dispersion_tail, cache.βu, cache.Γhat, cache.C, cache.Ĉ, cache.Γ_new
    
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
    ctx = UnchargedClosureEvalContext(r, mayer_f, Γ_old, βu, βu_disp_tail)

    while err > tolerance 
        if iteration > max_iterations
            error("Recursive iteration did not converge within $iteration steps. Current error = $err.")
        end
        closure_apply!(C, closure, ctx)
        oz_iteration_step!(C, Ĉ, Γhat, Γ_new, ρ, fourierplan)
        err = compute_error(Γ_new, Γ_old)
        if isnan(err)
            error("Encountered NaN in the solution after $iteration iterations.")
        end
        if method.verbose && iteration % 100 == 0
            println("After iteration $iteration, the error is $(round(err, digits=ceil(Int, 5-log10(tolerance)))).")
        end
        @. Γ_new = mixing_parameter * Γ_new + (1.0 - mixing_parameter) * Γ_old
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

# definitions:
# βu = βu_short_range + βu_LR_disp + βu_SR_coul + βu_LR_coul (first two are WCA splitting latter are coul splitting)
# c = c_short_range + Φ (long range part is just Φ=-βu_LR_coul)
# γ = γ_short_range + γ_long_range 
# h = h_short_range + q (defined through long range oz relation: q = Φ + Φ ρ q (fourier))
function solve(system::SimpleChargedSystem, closure::Closure, method::FourierIteration; gamma_0=nothing, coulombsplitting=NoCoulombSplitting())
    ρ = ρ_of(system)

    ensure_dispersion_support(closure, base_of(system).potential)

    # construct solver cache from base system. Does not include Coulomb parts
    cache = OZSolverCache(base_of(system), method)
    mayer_f, fourierplan, r, k, βu_LR_disp, βu, Γ_SR_hat, C_SR, C_SR_hat, Γ_SR_new = 
        cache.mayer_f, cache.fourierplan, cache.r, cache.k, cache.βu_dispersion_tail, cache.βu, cache.Γhat, cache.C, cache.Ĉ, cache.Γ_new
    
    βu_SR_coul, βu_LR_coul = split_coulomb_potential(r, system, coulombsplitting)
    βu = βu .+ βu_SR_coul .+ βu_LR_coul
    mayer_f .= find_mayer_f_function.(βu)

    φ = -βu_LR_coul # long range part of c 
    Φ_hat = fourier(φ .* r, fourierplan)  # note that the FT work on the "mulr" functions
    φ_hat = Φ_hat ./ k

    q_hat = [(I - φ_hat[i] * ρ) \ φ_hat[i] for i in eachindex(φ_hat)] # q_hat = (I - φ_hat * ρ)^(-1) * φ_hat
    Q_hat = q_hat .* k # mul_r function
    Q = inverse_fourier(Q_hat, fourierplan)
    q = Q ./ r

    # capital letters are multiplied by r or k
    Γ_SR_old = copy(mayer_f); C_hat = 0.0copy(mayer_f); Γ_hat = 0.0copy(mayer_f) 

    # initial guess for short ranged gamma    
    if isnothing(gamma_0)
        fill!(Γ_SR_old, zero(eltype(Γ_SR_old)))
    else
        Γ_SR_old .= gamma_0 .* r .- Q .+ Φ # gamma_0 is full gamma, we need short ranged part
    end

    max_iterations = method.max_iterations
    tolerance = method.tolerance
    mixing_parameter = method.mixing_parameter
    err = tolerance * 2;  iteration = 0
    ctx = ChargedClosureEvalContext(r, mayer_f, Γ_SR_old, βu, βu_LR_disp, βu_LR_coul, q)

    while err > tolerance 
        if iteration > max_iterations || isnan(err)
            error("Recursive iteration did not converge within $iteration steps. Current error = $err.")
        end
        closure_apply!(C_SR, closure, ctx)
        fourier!(C_SR_hat, C_SR, fourierplan)
        C_hat .= C_SR_hat .+ Φ_hat
        for ik in eachindex(Γ_hat)
            Γ_hat[ik] = (I*k[ik] - C_hat[ik] * ρ) \ (C_hat[ik] * ρ * C_hat[ik]) 
        end
        Γ_SR_hat .= Γ_hat .- (Q_hat .- Φ_hat) 
        inverse_fourier!(Γ_SR_new, Γ_SR_hat, fourierplan)
        err = compute_error(Γ_SR_new, Γ_SR_old) # compute error on mul_r functions, consistent with noncharged case
        if isnan(err)
            error("Encountered NaN in the solution after $iteration iterations.")
        end
        # mixing 
        @. Γ_SR_new = mixing_parameter * Γ_SR_new + (1.0 - mixing_parameter) * Γ_SR_old 

        if method.verbose && iteration % 100 == 0
            println("After iteration $iteration, the error is $(round(err, digits=ceil(Int, 5-log10(tolerance)))).")
        end
        Γ_SR_old .= Γ_SR_new
        iteration += 1
    end

    if method.verbose
        print("Converged after $iteration iterations, ")
        println("the error is $(round(err, digits=ceil(Int, 1-log10(err)))).")
    end

    C = @. C_SR/r + φ
    @. C_hat = C_SR_hat/r + φ_hat
    Γ = @. Γ_SR_new/r + q - φ
    @. Γ_hat = Γ_SR_hat/k + q_hat - φ_hat

    return construct_solution(r, k, C, C_hat, Γ, Γ_hat, ρ)
end

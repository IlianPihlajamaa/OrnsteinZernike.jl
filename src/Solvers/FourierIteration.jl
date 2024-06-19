
function solve(system::SimpleLiquid{dims, species, T1, T2, P}, closure::Closure, method::FourierIteration; init=nothing) where {dims, species, T1, T2, P}
    ρ = system.ρ

    r = method.dr * (1:method.M) |> collect

    βu1, _ = evaluate_long_range_potential(system.potential, system.kBT, r[1])
    elementtype = typeof(r[1] .* system.kBT .* system.ρ .* βu1)
    mayer_f = zeros(elementtype, length(r))
    fourierplan = get_fourier_plan(system, method, mayer_f)
    r .= fourierplan.r # in the case that dims != 3, we need to use the right grid
    k = fourierplan.k
    βu, βu_long_range = evaluate_long_range_potential(system.potential, system.kBT, r)
    mayer_f .= find_mayer_f_function.((system,), βu)

    Ĉ = copy(mayer_f) #ĉ*k
    Γ_new = copy(mayer_f) #γ̂ *k
    Γ_old = copy(mayer_f)*0.0 #γ*r
    if !(isnothing(init))
        Γ_old .= init.*r
    else
        Γ_old .= (zero(eltype(Γ_old)), )
    end
    Γhat = copy(mayer_f)
    C = copy(mayer_f) #c*r
    C_old = copy(mayer_f)
        
    max_iterations = method.max_iterations
    tolerance = method.tolerance
    mixing_parameter = method.mixing_parameter

    err = tolerance*2
    iteration = 0

    while err > tolerance
        if iteration > max_iterations
            error("Recursive iteration did not converge within $iteration steps. Current error = $err.")
        end
        C .= closure_cmulr_from_gammamulr.((closure, ), r, mayer_f, Γ_old, βu_long_range)
        if iteration != 0
            @. C = mixing_parameter * C + (1.0 - mixing_parameter) * C_old
        end

        fourier!(Ĉ, C, fourierplan)
        for ik in eachindex(Γhat, Ĉ)
            Γhat[ik] = (k[ik] * I - Ĉ[ik]*ρ) \ (Ĉ[ik] * ρ * Ĉ[ik])
        end
        inverse_fourier!(Γ_new, Γhat, fourierplan)
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
    g = find_g_from_c_and_Γ(c, Γ_new, r)
    ĉ = Ĉ ./ k
    Sk =  find_S_from_ĉ_and_ρ(ĉ, ρ)
    return OZSolution(r, k, g, Sk, ĉ, c)
end 


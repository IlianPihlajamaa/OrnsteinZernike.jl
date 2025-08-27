
function solve(system::SimpleUnchargedSystem, closure::Closure, method::FourierIteration; init=nothing) 
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



# === Generic helpers ===========================================================

"""
    add_numeric_LR_hat!(Ĉ_LR, βu_LR, r, plan)

Accumulate the k-space long-range piece using your FFT convention:
- Input: βu_LR(r) = long-range βu in real space.
- We set c_LR(r) = -βu_LR(r), then transform (c_LR*r) → (ĉ_LR*k) with `fourier!`.
- Output is accumulated into Ĉ_LR (which stores ĉ*k).
"""
function add_numeric_LR_hat!(Ĉ_LR, βu_LR::AbstractVector, r::AbstractVector, plan)
    tmp = similar(Ĉ_LR)
    @. tmp = (-βu_LR) * r          # C_LR = c_LR(r)*r
    fourier!(tmp, tmp, plan)       # tmp := (ĉ_LR*k)
    @. Ĉ_LR += tmp
    return Ĉ_LR
end

"""
    coulomb_split_hat!(Ĉ_LR, z, ℓB, κ, k)

Accumulate the analytic Coulomb long-range piece (3D Gaussian/Ewald split):
ĉ_LR(k) = -4π z^2 ℓB * exp(-k^2/(4κ^2)) / k^2  ⇒ we store (ĉ_LR * k).
The k=0 bin is set to 0 (neutralizing background removes the zero mode).
"""
function coulomb_split_hat!(Ĉ_LR, z, ℓB, κ, k::AbstractVector)
    @inbounds for i in eachindex(k)
        ki = k[i]
        Ĉ_LR[i] += iszero(ki) ? zero(eltype(Ĉ_LR)) :
                    -(4π * z^2 * ℓB) * exp(-(ki^2) / (4κ^2)) / ki
    end
    return Ĉ_LR
end

"""
    coulomb_real_parts(z, ℓB, κ, r)

Return (βu_Coul_sr, βu_Coul_lr_real) for the Gaussian/Ewald split in real space.
"""
function coulomb_real_parts(z, ℓB, κ, r::AbstractVector)
    β = one(eltype(r)) # β is already included in ℓB; keep types consistent
    βu_Coul_sr      = (z^2 * ℓB) .* erfc.(κ .* r) ./ r
    βu_Coul_lr_real = (z^2 * ℓB) .* erf.(κ .* r)  ./ r
    return βu_Coul_sr, βu_Coul_lr_real
end

# === Charged solver for arbitrary base potentials ==============================

"""
    solve(system::SimpleChargedFluid, closure::Closure, method::FourierIteration;
          init=nothing, enforce_SL=true, use_numeric_coulomb=false)

General OZ–closure iteration with *arbitrary* base potential + Coulomb background.

- The closure only sees the **short-range** potential (base_sr + Coul_sr).
- In k-space we add LR pieces:  Ĉ_total = Ĉ_sr + Ĉ_LR_base + Ĉ_LR_coul.
- For the base potential, the LR k-piece is obtained **numerically** from βu_base_LR (safe & generic).
- For Coulomb we use the **analytic** expression (set `use_numeric_coulomb=true` to force numeric FT).
"""
function solve(system::SimpleChargedFluid, closure::Closure, method::FourierIteration;
               init=nothing, enforce_SL::Bool=true, use_numeric_coulomb::Bool=false)

    ρ = system.base.ρ

    # Grids & Fourier plan (inherits dims grid from plan)
    r = method.dr * (1:method.M) |> collect
    βu_probe, _ = evaluate_long_range_potential(system.base.potential, system.base.kBT, r[1])
    T = typeof(r[1] * system.base.kBT * ρ * βu_probe)
    work = zeros(T, length(r))
    plan = get_fourier_plan(system.base, method, work)
    r  .= plan.r
    k   =  plan.k

    # Base potential LR/SR via existing hook
    βu_base_total, βu_base_LR = evaluate_long_range_potential(system.base.potential, system.base.kBT, r)
    βu_base_SR = βu_base_total .- βu_base_LR

    # Coulomb split (system-level)
    βu_Coul_SR, βu_Coul_LR = coulomb_real_parts(system.z, system.ℓB, system.κ, r)

    # Mayer f and the long-range subtraction passed to the closure
    βu_total_for_mayer = βu_base_total .+ βu_Coul_SR .+ βu_Coul_LR
    βu_LR_for_closure  = βu_base_LR     .+               βu_Coul_LR
    mayer_f = similar(work);  mayer_f .= find_mayer_f_function.((system.base,), βu_total_for_mayer)

    # Work arrays (convention: C=c*r, Ĉ=ĉ*k, Γ=γ*r, Γ̂=γ̂*k)
    Ĉ    = copy(mayer_f)
    Γ_new = copy(mayer_f)
    Γ_old = zero.(mayer_f)
    if init !== nothing
        @. Γ_old = init * r
    end
    Γhat  = copy(mayer_f)
    C     = copy(mayer_f)
    Cold  = copy(mayer_f)

    # Prebuild total LR in k-space: base (numeric) + Coulomb (analytic or numeric)
    Ĉ_LR = zero.(mayer_f)                       # accumulates ĉ_LR*k
    add_numeric_LR_hat!(Ĉ_LR, βu_base_LR, r, plan)
    if use_numeric_coulomb
        add_numeric_LR_hat!(Ĉ_LR, βu_Coul_LR, r, plan)
    else
        coulomb_split_hat!(Ĉ_LR, system.z, system.ℓB, system.κ, k)
    end

    maxit = method.max_iterations
    tol   = method.tolerance
    mix   = method.mixing_parameter

    err, it = tol*2, 0
    while err > tol
        it > maxit && error("Charged OZ iteration did not converge in $it steps (err=$err).")

        # Closure on SR potential (base_SR + Coul_SR); LR provided to subtract inside closure
        C .= closure_cmulr_from_gammamulr.((closure,), r, mayer_f, Γ_old, βu_LR_for_closure)
        if it != 0
            @. C = mix*C + (1 - mix)*Cold
        end

        # Fourier & add LR
        fourier!(Ĉ, C, plan)             # ⇒ Ĉ_sr = (ĉ_sr*k)
        @. Ĉ = Ĉ + Ĉ_LR                 # ⇒ Ĉ_total

        # OZ in k-space (γ̂*k)
        @inbounds for i in eachindex(Γhat, Ĉ)
            Γhat[i] = (k[i] * I - ρ*Ĉ[i]) \ (ρ * Ĉ[i] * Ĉ[i])
        end

        # Optional Stillinger–Lovett (only meaningful if Coulomb present)
        if enforce_SL && (abs(system.z) > 0)
            κD2 = 4π * system.ℓB * (system.z^2) * ρ
            i1  = first(eachindex(k))
            k1  = k[i1]
            if k1 > 0
                S_target      = (k1^2) / κD2
                Hhat_target_k = (S_target - 1) * k1 / ρ
                Γhat[i1] = Hhat_target_k - Ĉ[i1]
            end
        end

        inverse_fourier!(Γ_new, Γhat, plan)  # back to real space

        err = compute_error(Γ_new, Γ_old)
        if method.verbose && it % 100 == 0
            println("After iteration $it, error = $(round(err, digits=ceil(Int, 5 - log10(tol)))).")
        end
        Γ_old .= Γ_new
        Cold  .= C
        it += 1
    end
    if method.verbose
        println("Converged after $it iterations, error = $(round(err, digits=ceil(Int, 1 - log10(err)))).")
    end

    # Outputs: real-space c reported as SR piece (the true c has LR tails); k-space has the full ĉ.
    c_sr = C ./ r
    ĉ    = Ĉ ./ k
    g     = find_g_from_c_and_Γ(c_sr, Γ_new, r)     # your closure helper expects SR c(r)
    Sk    = find_S_from_ĉ_and_ρ(ĉ, ρ)

    return OZSolution(r, k, g, Sk, ĉ, c_sr)
end

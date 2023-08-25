using OrnsteinZernike, Plots
import Roots

potential = HardSpheres(1.0)
dims = 3 
ρ = 0.1*sqrt(2)
kBT = 1.0 
system = SimpleLiquid(dims, ρ, kBT, potential)
closure = HypernettedChain()
method = DensityRamp(NgIteration(), (0.1:0.2:0.9)*ρ)
sol = solve(system, closure, method)
p = compute_virial_pressure(sol, system)
χ = compute_compressibility(sol, system)



# using Plots
# plot(sol.r, sol.gr, xlims=(0,5), xlabel="r", ylabel="g(r)", lw=1, label="RY")
# @show maximum(sol.gr)
# sol = solve(system, PercusYevick(), method)
# plot!(sol.r, sol.gr, xlims=(0,5), xlabel="r", label="PY")
# sol = solve(system, HypernettedChain(), method)
# plot!(sol.r, sol.gr, xlims=(0,5), xlabel="r", label="HNC", lw=1)


function solve_self_consistent_RY(ρ)
    dims = 3 
    kBT = 1.0 
    β = 1/kBT
    function solve_RY(α, ρ)
        potential = HardSpheres(1.0)
        system = SimpleLiquid(dims, ρ, kBT, potential)
        closure = RogersYoung(α)
        method = DensityRamp(NgIteration(verbose=false, M=2^12, max_iterations=1000, tolerance=10^-8), (0.01:0.4:0.9)*ρ; verbose=false)
        sol = solve(system, closure, method)
        return system, sol
    end

    function find_thermodynamic_inconsistency(α, ρ)
        sys1, sol1 = solve_RY(α, ρ)
        χ = compute_compressibility(sol1, sys1)

        p1 = compute_virial_pressure(sol1, sys1)
        dρ = sqrt(eps(ρ))
        sys2, sol2 = solve_RY(α, ρ+dρ)
        p2 = compute_virial_pressure(sol2, sys2)
        βdpdρ = β*(p2-p1)/dρ
        inconsistency =  βdpdρ - β/(ρ*χ)

        @show inconsistency, α
        return inconsistency
    end

    fn = α -> find_thermodynamic_inconsistency(α, ρ)
    α =  Roots.find_zero(fn, (0.0000001,5), Roots.Bisection(), atol=0.01)
    sys, sol = solve_RY(α, ρ)
    return α, sys, sol
end

ρ = 0.6*sqrt(2)
α, sys, sol =  solve_self_consistent_RY(ρ)
p = compute_virial_pressure(sol, sys)

@show α, p/ρ-1
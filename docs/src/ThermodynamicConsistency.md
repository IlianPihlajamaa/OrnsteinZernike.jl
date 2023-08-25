# Thermodynamic Consistency

Some [Closures](@ref) have free parameters that can be tuned to obtain the most accurate bridge functions. For example, the famous Rogers-Young closure [`RogersYoung`](@ref) interpolates between the Percus-Yevick closure and the Hypernetted Chain one with a free parameter $\alpha$. While this parameter can be chosen freely, it is usually chosen such that the solutions satisfy thermodynamic consistency. 

## Free choice of $\alpha$

To warm up, let's solve a 3D hard-sphere system for $\alpha=0.5$
```@example 1
using OrnsteinZernike, Plots
M = 5000
ρ = 0.6 * sqrt(2)
dr = 20.0/M
kBT = 1.0
dims = 3 
α = 0.5
pot = HardSpheres(1.0)
system = SimpleLiquid(dims, ρ, kBT, pot)
method = NgIteration(M=M, dr=dr, verbose=false)
solRY = solve(system, RogersYoung(α), method)
solPY = solve(system, PercusYevick(), method)
solHNC = solve(system, HypernettedChain(), method)

plot(solRY.r, solRY.gr, label="Rogers-Young")
plot!(solPY.r, solPY.gr, label="Percus-Yevick")
plot!(solHNC.r, solHNC.gr, label="Hypernetted Chain")
plot!(xlims=(0.9,2.3), xlabel="r", ylabel="g(r)")
```
Here, we can see, that indeed, the Rogers-Young closure finds a middle ground between the two others. 

## Using thermodynamics

There are a number of thermodynamic relations that allow a consistent choice of $\alpha$. Here we opt for the requirement $\frac{1}{\rho\chi_T}=\left(\frac{\partial p}{\partial \rho}\right)_T$, which physically means that we require the virial and compressibility route of obtaining the pressure to give the same results. For simplicity, we use finite differences to obtain the derivative of the pressure. See the [Theory](@ref) section for more details. Subsequently, we use the bisection method from `Roots.jl` to find the optimal parameter. The results can be compared to Table 1 in Ref. 1. 



```@example 2
using OrnsteinZernike, Plots
import Roots
function find_self_consistent_solution(ρ, kBT, M, dr, dims, pot)

    function RY_inconsistency(ρ, α)
        system1 = SimpleLiquid(dims, ρ, kBT, pot)
        method = NgIteration(M=M, dr=dr, verbose=false)
        sol1 = solve(system1, RogersYoung(α), method)
        p1 = compute_virial_pressure(sol1, system1)

        dρ = sqrt(eps(ρ))
        system2 = SimpleLiquid(dims, ρ+dρ, kBT, pot)
        sol2 = solve(system2, RogersYoung(α), method)
        p2 = compute_virial_pressure(sol2, system2)
        dpdρ = (p2-p1)/dρ

        χ = compute_compressibility(sol1, system1)
        inconsistency = dpdρ - 1/(ρ*χ)
        return inconsistency
    end

    # we need to find α such that the inconsistency is zero.
    func = α ->  RY_inconsistency(ρ, α)
    α =  Roots.find_zero(func, (0.001,50.0), Roots.Bisection(), atol=0.0001)
    system = SimpleLiquid(dims, ρ, kBT, pot)
    method = NgIteration(M=M, dr=dr, verbose=false)
    sol = solve(system, RogersYoung(α), method)
    return system, sol, α
end

for ρstar = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.654]
    ρ = ρstar*sqrt(2)
    M = 1000
    dr = 10.0/M
    kBT = 1.0
    dims = 3 
    pot = HardSpheres(1.0)
    system, sol, α = find_self_consistent_solution(ρ, kBT, M, dr, dims, pot)
    P = compute_virial_pressure(sol, system)/ρ/kBT - 1
    gmax = maximum(sol.gr)
    println("At ρ/√2 = $(ρstar), we find α = $(round(α,digits=2)), and βp/ρ - 1 = $(round(P,digits=4)).")
end
```

References:

1. Rogers, Forrest J., and David A. Young. "New, thermodynamically consistent, integral equation for simple fluids." Physical Review A 30.2 (1984): 999.
import Pkg; Pkg.activate(".")
using Revise
using OrnsteinZernike,  Plots
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
        inconsistency = dpdρ/kBT - 1/(ρ*kBT*χ)

        return inconsistency
    end

    func = α ->  RY_inconsistency(ρ, α)
    α =  Roots.find_zero(func, (0.1,5.0), Roots.Bisection(), atol=0.0001)
    system = SimpleLiquid(dims, ρ, kBT, pot)
    method = NgIteration(M=M, dr=dr, verbose=false)
    sol = solve(system, RogersYoung(α), method)
    return system, sol, α
end
println("Hard Spheres")
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
    println("At ρ/√2 = $(ρstar), we find α = $(trunc(α,digits=2)), and βp/ρ - 1 = $(trunc(P,digits=4)).")
end

## 12 power fluid
println("1/r^12 fluid")

for z = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.813]
    ρ = z*sqrt(2)
    M = 1000
    dr = 20.0/M
    kBT = 1.0
    dims = 3 
    Γ = (4π*sqrt(2)*z/3)^(4)
    σ = 1.0
    ϵ = 1.0

    pot = PowerLaw(ϵ, σ, 12)
    system, sol, α = find_self_consistent_solution(ρ, kBT, M, dr, dims, pot)
    P = compute_virial_pressure(sol, system)/ρ/kBT - 1
    gmax = maximum(sol.gr)
    println("At z = $(z), we find Γ = $(trunc(Γ,digits=2)) α = $(trunc(α,digits=2)), and βp/ρ - 1 = $(trunc(P,digits=4)).")
end

println("1/r^9 fluid")

for z = [0.1, 0.25, 0.5, 0.943]
    ρ = z*sqrt(2)
    M = 1000
    dr = 20.0/M
    kBT = 1.0
    dims = 3 
    Γ = (4π*sqrt(2)*z/3)^(3)
    σ = 1.0
    ϵ = 1.0

    pot = PowerLaw(ϵ, σ, 9)
    system, sol, α = find_self_consistent_solution(ρ, kBT, M, dr, dims, pot)
    P = compute_virial_pressure(sol, system)/ρ/kBT - 1
    gmax = maximum(sol.gr)
    println("At z = $(z), we find Γ = $(trunc(Γ,digits=2)) α = $(trunc(α,digits=2)), and βp/ρ - 1 = $(trunc(P,digits=4)).")
end

println("1/r^6 fluid")

for z = [0.1, 0.25, 0.5, 1.0, 1.54]
    ρ = z*sqrt(2)
    M = 2000
    dr = 20/M
    kBT = 1.0
    dims = 3
    Γ = (4π*sqrt(2)*z/3)^(2)
    σ = 1.0
    ϵ = 1.0

    pot = PowerLaw(ϵ, σ, 6)
    system, sol, α = find_self_consistent_solution(ρ, kBT, M, dr, dims, pot)
    P = compute_virial_pressure(sol, system)/ρ/kBT - 1
    gmax = maximum(sol.gr)
    println("At z = $(z), we find Γ = $(trunc(Γ,digits=2)) α = $(trunc(α,digits=2)), and βp/ρ - 1 = $(trunc(P,digits=4)).")
end
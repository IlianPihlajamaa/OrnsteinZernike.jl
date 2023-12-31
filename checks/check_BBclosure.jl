import Pkg; Pkg.activate(".")
using Revise
using OrnsteinZernike,  Plots, Dierckx
import Roots, ForwardDiff

function find_self_consistent_solution(ρ, kBT, M, dr, dims, pot)

    function pressure(ρ, α)
        method = NgIteration(M=M, dr=dr, verbose=false)
        system = SimpleLiquid(dims, ρ, kBT, pot)
        sol = solve(system, BomontBretonnet(α), method)
        p = compute_virial_pressure(sol, system)
        return p
    end

    function find_inconsistency(ρ, α)
        system1 = SimpleLiquid(dims, ρ, kBT, pot)
        method = NgIteration(M=M, dr=dr, verbose=false)
        sol1 = solve(system1, BomontBretonnet(α), method)

        dpdρ = @time ForwardDiff.derivative(ρ -> pressure(ρ, α), ρ)

        @time begin
            p1 = compute_virial_pressure(sol1, system1)
            dρ = sqrt(eps(ρ))
            system2 = SimpleLiquid(dims, ρ+dρ, kBT, pot)
            sol2 = solve(system2, BomontBretonnet(α), method)
            p2 = compute_virial_pressure(sol2, system2)
            dpdρ2 = (p2-p1)/dρ
        end
        @time begin
            dρ = sqrt(eps(ρ))
            system2 = SimpleLiquid(dims, ρ+dρ, kBT, pot)
            sol2 = solve(system2, BomontBretonnet(α), method)
            p2 = compute_virial_pressure(sol2, system2)
            system3 = SimpleLiquid(dims, ρ-dρ, kBT, pot)
            sol3 = solve(system3, BomontBretonnet(α), method)
            p3 = compute_virial_pressure(sol3, system3)
            dpdρ3 = (p2-p3)/dρ/2
        end
        @show dpdρ, dpdρ2, dpdρ3

        χ = compute_compressibility(sol1, system1)
        inconsistency = dpdρ/kBT - 1/(ρ*kBT*χ)

        return inconsistency
    end

    func = α ->  find_inconsistency(ρ, α)
    α =  Roots.find_zero(func, (0.0,1.0), Roots.Bisection(), atol=0.0001)
    system = SimpleLiquid(dims, ρ, kBT, pot)
    method = NgIteration(M=M, dr=dr, verbose=false)
    sol = solve(system, BomontBretonnet(α), method)
    return system, sol, α
end
println("Hard Spheres")
for ρstar = [0.3, 0.5, 0.7, 0.8, 0.9]
    ρ = ρstar
    M = 1000
    dr = 10.0/M
    kBT = 1.0
    dims = 3
    pot = HardSpheres(1.0)
    system, sol, α = find_self_consistent_solution(ρ, kBT, M, dr, dims, pot)

    P = compute_virial_pressure(sol, system)/ρ/kBT - 1
    gmax = maximum(sol.gr)
    B = (OrnsteinZernike.bridge_function(BomontBretonnet(α), 0.0, 0.0, sol.gr .- 1 .- sol.cr))
    plot(sol.r, B) |> display
    println("At ρ = $(ρstar), we find α = $(trunc(α,digits=4)), and βp/ρ - 1 = $(trunc(P,digits=4)).")
    println("B(0) = $(B[1])")
end

for closure in [PercusYevick, MartynovSarkisov]
    for ρ in [0.3, 0.5, 0.7, 0.8, 0.9]
        M = 1000000
        dr = 10.0/M
        kBT = 1.0
        dims = 3
        pot = HardSpheres(1.0)
        system = SimpleLiquid(dims, ρ, kBT, pot)
        method = NgIteration(M=M, dr=dr, verbose=false)
        sol = solve(system, closure(), method)
        B = (OrnsteinZernike.bridge_function(closure(), 0.0, 0.0, sol.gr .- 1 .- sol.cr))

        println("for closure $(closure) at ρ = $(ρ), we have B(0) = $(Spline1D(sol.r,B)(0.0))")
    end
end




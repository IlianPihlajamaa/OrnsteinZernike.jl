function find_self_consistent_solution(ρ, kBT, M, dr, dims, pot)

    function RY_inconsistency(ρ, α)
        system1 = SimpleFluid(dims, ρ, kBT, pot)
        method = NgIteration(M=M, dr=dr, verbose=false)
        sol1, = solve(system1, RogersYoung(α), method)
        p1 = compute_virial_pressure(sol1, system1)

        dρ = sqrt(eps(ρ))
        system2 = SimpleFluid(dims, ρ+dρ, kBT, pot)
        sol2, = solve(system2, RogersYoung(α), method)
        p2 = compute_virial_pressure(sol2, system2)
        dpdρ = (p2-p1)/dρ

        χ = compute_compressibility(sol1, system1)
        inconsistency = dpdρ/kBT - 1/(ρ*kBT*χ)

        return inconsistency
    end

    func = α ->  RY_inconsistency(ρ, α)
    α =  Roots.find_zero(func, (0.1,5.0), Roots.Bisection(), atol=0.0001)
    system = SimpleFluid(dims, ρ, kBT, pot)
    method = NgIteration(M=M, dr=dr, verbose=false)
    sol, = solve(system, RogersYoung(α), method)
    return system, sol, α
end

for ρstar = [0.3]
    ρ = ρstar*sqrt(2)
    M = 1000
    dr = 10.0/M
    kBT = 1.0
    dims = 3
    pot = HardSpheres(1.0)
    system, sol, α = find_self_consistent_solution(ρ, kBT, M, dr, dims, pot)
    P = compute_virial_pressure(sol, system)/ρ/kBT - 1
    gmax = maximum(sol.gr)
    @test abs(P -  1.664) < 0.05
end


for dims = [1, 2, 3, 4, 5]
    M = 2^10
    ρ = 0.3
    kBT = 1.0

    pot = SingleComponentHardSpheres()
    system = SimpleLiquid(dims, ρ, kBT, pot)
    closure = PercusYevick()
    method = NgIteration(tolerance=10^-10, N_stages=5, M=M, verbose=false, max_iterations=10^3)
    sol = solve(system, closure, method)
    @test all(isfinite(sol.gr))
end
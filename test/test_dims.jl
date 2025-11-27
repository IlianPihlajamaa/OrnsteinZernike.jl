
function dimstest()
    for dims = [1, 2, 3, 4, 5]
        M = 2^10
        ρ = 0.8
        kBT = 1.0

        pot = HardSpheres(1.0)
        system = SimpleFluid(dims, ρ, kBT, pot)
        closure = PercusYevick()
        method = NgIteration(tolerance=10^-10, N_stages=5, M=M, verbose=false, max_iterations=10^3)
        sol, = solve(system, closure, method)
        @test all(isfinite.(sol.gr)) # should test with exact method
    end
end

dimstest()
## 1d M. E. Fisher and B. Widom, J. Chem. Phys. 50, 3756 (1969)
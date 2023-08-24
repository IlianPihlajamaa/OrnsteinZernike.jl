

for dims = [1, 2, 3, 4, 5] # dims = 1 does not converge
    local M = 2^10
    local ρ = 0.8
    local kBT = 1.0

    local pot = HardSpheres(1.0)
    local system = SimpleLiquid(dims, ρ, kBT, pot)
    local closure = PercusYevick()
    local method = DensityRamp(NgIteration(tolerance=10^-10, N_stages=5, M=M, verbose=false, max_iterations=10^3), ρ*(0.001:0.2:1.0); verbose=false)
    local sol = solve(system, closure, method)
    @test all(isfinite.(sol.gr)) # should test with exact method
end

## 1d M. E. Fisher and B. Widom, J. Chem. Phys. 50, 3756 (1969)
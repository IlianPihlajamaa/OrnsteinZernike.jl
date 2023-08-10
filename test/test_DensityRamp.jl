


M = 2^10
ρ = 1.0
kBT = 1.0
dims = 3

pot = SingleComponentHardSpheres()
system = SimpleLiquid(dims, ρ, kBT, pot)
closure = PercusYevick()
method0 = NgIteration(tolerance=10^-10, M=M, N_stages=5, max_iterations=1000, verbose=false)
method = DensityRamp(method0, ρ*(0.1:0.1:1.0), verbose=false)
sol = solve(system, closure, method);

method2 = Exact(M=M)
sol2 = solve(system, closure, method2);

atol = 0.1
@test all(abs.(sol.cr .- sol2.cr) .< 10atol) 
@test all(abs.(sol.gr .- sol2.gr) .< atol) 
@test all(abs.(sol.ck .- sol2.ck) .< 10atol) 
@test all(abs.(sol.Sk .- sol2.Sk) .< atol) 
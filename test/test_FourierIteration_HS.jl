
# Single Component
M = 2^10
ρ = 0.5
kBT = 1.1
dims = 3

pot = HardSpheres(1.0)
system = SimpleLiquid(dims, ρ, kBT, pot)
closure = PercusYevick()
method = FourierIteration(M = M, tolerance=10^-10, mixing_parameter=0.8, verbose=false, max_iterations=10^4)
sol = solve(system, closure, method)

sol2 = solve(system, closure, Exact(M=M))


atol = 0.1
@test all(abs.(sol.cr .- sol2.cr) .< atol) 
@test all(abs.(sol.gr .- sol2.gr) .< atol) 
@test all(abs.(sol.ck .- sol2.ck) .< 1.0) 
@test all(abs.(sol.Sk .- sol2.Sk) .< atol) 


M = 2^10
ρ = [0.1, 0.2, 0.12]
kBT = 1.0
dims = 3

pot = HardSpheres([1.0, 1.2, 0.8])
system = SimpleLiquid(dims, ρ, kBT, pot)

closure = PercusYevick()
method = FourierIteration(M = M, tolerance=10^-10, verbose=false, max_iterations=10^4)
sol =  solve(system, closure, method)

sol2 = solve(system, closure, Exact(M=M))

atol = 0.3
@test all((abs.(sol.cr .- sol2.cr)) .< atol) 
@test all((abs.(sol.gr .- sol2.gr)) .< atol) 
@test all((abs.(sol.ck .- sol2.ck)) .< atol) 
@test all((abs.(sol.Sk .- sol2.Sk)) .< atol) 
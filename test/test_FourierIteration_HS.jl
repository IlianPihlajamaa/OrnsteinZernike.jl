
# Single Component

M = 2^12
ρ = 0.5
kBT = 1.1
dims = 3

pot = SingleComponentHardSpheres()
system = SimpleLiquid(dims, ρ, kBT, pot)
closure = PercusYevick()
method = FourierIteration(M = M, tolerance=10^-10, mixing_parameter=0.8)
sol = solve(system, closure, method)

sol2 = solve(system, closure, Exact(M))


atol = 0.1
@test all(abs.(sol.cr .- sol2.cr) .< atol) 
@test all(abs.(sol.gr .- sol2.gr) .< atol) 
@test all(abs.(sol.ck .- sol2.ck) .< atol) 
@test all(abs.(sol.Sk .- sol2.Sk) .< atol) 




M = 2^14-1
ρ = [0.2, 0.3, 0.3]
kBT = 1.1
dims = 3

pot = MultiComponentHardSpheres([1.0, 1.2, 0.8])
system = SimpleLiquid(dims, ρ, kBT, pot)

closure = PercusYevick()
method = FourierIteration(M = M, tolerance=10^-10)
sol =  solve(system, closure, method)

sol2 = solve(system, closure, Exact(M))

atol = 0.1
@test all(maximum.(abs.(sol.cr .- sol2.cr)) .< atol) 
@test all(maximum.(abs.(sol.gr .- sol2.gr)) .< atol) 
@test all(maximum.(abs.(sol.ck .- sol2.ck)) .< atol) 
@test all(maximum.(abs.(sol.Sk .- sol2.Sk)) .< atol) 
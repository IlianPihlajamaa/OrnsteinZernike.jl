
dims = 3

M = 2^10
ρ = 0.3
kBT = 1.0

pot = SingleComponentHardSpheres()
system = SimpleLiquid(dims, ρ, kBT, pot)
closure = PercusYevick()
method = NgIteration(tolerance=10^-10, N_stages=5, M=M, verbose=false, max_iterations=10^3)
sol = solve(system, closure, method)
sol2 = solve(system, closure, Exact(M=M))

atol = 0.1

@test all((abs.(sol.cr .- sol2.cr)) .< atol) 
@test all((abs.(sol.gr .- sol2.gr)) .< atol) 
@test all((abs.(sol.ck .- sol2.ck)) .< atol) 
@test all((abs.(sol.Sk .- sol2.Sk)) .< atol) 


M = 1000
ρ = [0.08568282994118655,0.14359889826470624,0.056859386188983715,0.2203848611614371,0.26961494833839134]
kBT = 1.1
dims = 3

D = [0.40402349270703586,1.243480095653075,0.1173675810956587,1.327361311444319,0.3444445605976121]
pot = MultiComponentHardSpheres(D)
system = SimpleLiquid(dims, ρ, kBT, pot)
closure = PercusYevick()
method = NgIteration(tolerance=10^-10, N_stages=5, M=M, verbose=false, max_iterations=10^3)

sol = solve(system, closure, method)

method = Exact(M=M)
sol2 = solve(system, closure, method)

atol = 0.1
@test all(maximum.(abs.(sol.cr .- sol2.cr)) .< 10atol) 
@test all(maximum.(abs.(sol.gr .- sol2.gr)) .< atol) 
@test all(maximum.(abs.(sol.ck .- sol2.ck)) .< 10atol) 
@test all(maximum.(abs.(sol.Sk .- sol2.Sk)) .< atol) 


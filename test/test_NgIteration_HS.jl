function test2()
dims = 3

M = 1000
dr = 10/M
ρ = 0.3
kBT = 1.0

pot = HardSpheres(1.0)
system = SimpleFluid(dims, ρ, kBT, pot)
closure = PercusYevick()
method = NgIteration(tolerance=10^-10, N_stages=5, M=M, dr=dr, verbose=false, max_iterations=10^3)
sol = solve(system, closure, method)
sol2 = solve(system, closure, Exact(M=M, dr=dr))

atol = 0.1

@test all((abs.(sol.cr .- sol2.cr)) .< atol)
@test all((abs.(sol.gr .- sol2.gr)) .< atol)
@test all((abs.(sol.ck .- sol2.ck)) .< atol)
@test all((abs.(sol.Sk .- sol2.Sk)) .< atol)

# Check convergence info for NgIteration
@test sol.converged == true
@test sol.iterations > 0
@test sol.final_error < 10^-10
@test sol.termination_reason == :converged

# Check convergence info for Exact
@test sol2.converged == true
@test sol2.iterations == 0
@test sol2.final_error == 0.0
@test sol2.termination_reason == :exact

M = 2000
dr = 5.0/M
ρ = [0.08568282994118655,0.14359889826470624,0.056859386188983715,0.2203848611614371,0.26961494833839134]
kBT = 1.1
dims = 3

D = [0.40402349270703586,1.243480095653075,0.1173675810956587,1.327361311444319,0.3444445605976121]
pot = HardSpheres(D)
system = SimpleMixture(dims, ρ, kBT, pot)
closure = PercusYevick()
method = NgIteration(tolerance=10^-10,dr=dr, N_stages=5, M=M, verbose=false, max_iterations=10^3)

sol = solve(system, closure, method)

method = Exact(M=M, dr=dr,)
sol2 = solve(system, closure, method)


atol = 0.1


@test all(maximum.(abs.(sol.cr .- sol2.cr)) .< 10atol) 
@test all(maximum.(abs.(sol.gr .- sol2.gr)) .< atol) 
@test all(maximum.(abs.(sol.ck .- sol2.ck)) .< 10atol) 
@test all(maximum.(abs.(sol.Sk .- sol2.Sk)) .< atol) 

####
dims = 1

M = 2^8
ρ = 0.3
kBT = 1.0

pot = HardSpheres(1.0)
system = SimpleFluid(dims, ρ, kBT, pot)
closure = PercusYevick()
method = NgIteration(tolerance=10^-10, N_stages=5, M=M, verbose=false, max_iterations=10^3)
sol = solve(system, closure, method)
sol2 = solve(system, closure, Exact(M=M))

atol = 0.1

@test all((abs.(sol.cr .- sol2.cr)) .< atol) 
@test all((abs.(sol.gr .- sol2.gr)) .< atol) 
@test all((abs.(sol.ck .- sol2.ck)) .< atol) 
@test all((abs.(sol.Sk .- sol2.Sk)) .< atol) 


dims = 5

M = 2^8
ρ = 0.3
kBT = 1.0

pot = HardSpheres(1.0)
system = SimpleFluid(dims, ρ, kBT, pot)
closure = PercusYevick()
method = NgIteration(tolerance=10^-10, N_stages=5, M=M, verbose=false, max_iterations=10^3, dr=sqrt(π/(M+1))/(2π))
sol = solve(system, closure, method)
sol2 = solve(system, closure, Exact(M=M, dr=sqrt(π/(M+1))/(2π)))

atol = 0.1

@test all((abs.(sol.cr .- sol2.cr)) .< atol) 
@test all((abs.(sol.gr .- sol2.gr)) .< atol) 
@test all((abs.(sol.ck .- sol2.ck)) .< 5atol) 
@test all((abs.(sol.Sk .- sol2.Sk)) .< atol) 


end 

test2()

M = 1000
ρ = 1.0
kBT = 1.1
dims = 3
dr = 10/M
pot = HardSpheres(1.0)
system = SimpleFluid(dims, ρ, kBT, pot)
closure = PercusYevick()
method = FourierIteration(M = M,dr=dr, tolerance=10^-10, mixing_parameter=0.01, verbose=false, max_iterations=10^6)
sol = @time solve(system, closure, method)
sol = solve(system, closure, method)

sol2 = solve(system, closure, Exact(M=M,dr=dr,))


atol = 0.1
(abs.(sol.cr .- sol2.cr) ) |> maximum
(abs.(sol.gr .- sol2.gr) ) |> maximum
(abs.(sol.ck .- sol2.ck)) |> maximum
(abs.(sol.Sk .- sol2.Sk) ) |> maximum
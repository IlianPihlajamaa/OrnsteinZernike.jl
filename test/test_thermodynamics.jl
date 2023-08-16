

M = 1000
ρ = 0.5
kBT = 1.0
dims = 3

pot = HardSpheres(1.0)
system = SimpleLiquid(dims, ρ, kBT, pot)
closure = PercusYevick()
method = Exact(M=M)
sol = solve(system, closure, method)

method = NgIteration(M=M, tolerance=10^-10, verbose=false, max_iterations=10^3)
sol2 = solve(system, closure, method)


p = compute_virial_pressure(sol, system)
p2 = compute_virial_pressure(sol2, system)
χ = compute_compressibility(sol, system)
E = compute_excess_energy(sol, system)

η = ρ/6*π
pexact = ρ*kBT*(1+2η+3η^2)/(1-η)^2 
atol = 0.1
@test abs(pexact - p) < atol
@test abs(pexact - p2) < atol
@test E == 0.0

## add test for compressibility
## add test for multicomponent

import Pkg; Pkg.activate(".")

using Plots
import OrnsteinZernike
pl1 = plot()
pl2 = plot()

M = 1000
ρ = 0.5
kBT = 1.0
dims = 3

pot = OrnsteinZernike.SingleComponentHardSpheres()
system = OrnsteinZernike.SimpleLiquid(dims, ρ, kBT, pot)
closure = OrnsteinZernike.PercusYevick()
method = OrnsteinZernike.Exact(M=M)
sol = OrnsteinZernike.solve(system, closure, method)
method = OrnsteinZernike.FourierIteration(M=M, max_iterations=10000, tolerance=10^-10, init=zeros(M))
@time sol2 = OrnsteinZernike.solve(system, closure, method)

p = OrnsteinZernike.compute_virial_pressure(sol, system)
p2 = OrnsteinZernike.compute_virial_pressure(sol2, system)
χ = OrnsteinZernike.compute_compressibility(sol, system)
E = OrnsteinZernike.compute_excess_energy(sol, system)

η = ρ/6*π
pexact = ρ*kBT*(1+2η+3η^2)/(1-η)^2 
@show pexact, p, p2
@show findmin(sol.Cr),  findmin(sol2.Cr)
@show findmin(-sol.gr),  findmin(-sol2.gr)


# scatter!([log2(M)], [cr0], c="red", label=false)

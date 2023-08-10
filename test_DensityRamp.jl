import Pkg; Pkg.activate(".")

import OrnsteinZernike
using Plots



M = 2^12
ρ = 1.0
kBT = 1.0
dims = 3

pot = OrnsteinZernike.SingleComponentHardSpheres()
system = OrnsteinZernike.SimpleLiquid(dims, ρ, kBT, pot)
closure = OrnsteinZernike.PercusYevick()
method0 = OrnsteinZernike.NgIteration(tolerance=10^-10, M=M, N_stages=5, max_iterations=1000, verbose=false)
method = OrnsteinZernike.DensityRamp(method0, ρ*(0.1:0.01:1.0))
sol = OrnsteinZernike.solve(system, closure, method);

plot(sol.r, sol.gr , xlims=(0,6), ylims=(-1.1,10))
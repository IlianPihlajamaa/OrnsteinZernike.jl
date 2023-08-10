import Pkg; Pkg.activate(".")


import OrnsteinZernike


M = 2^12 - 1
ρ = 1.1
kBT = 1.0
dims = 3

pot = OrnsteinZernike.SingleComponentHardSpheres()
system = OrnsteinZernike.SimpleLiquid(dims, ρ, kBT, pot)
closure = OrnsteinZernike.PercusYevick()
method = OrnsteinZernike.Exact(M)
sol = OrnsteinZernike.solve(system, closure, method)
using Plots
plot(sol.k, sol.Ck) |> display
plot(sol.r, sol.Cr, xlims=(0,5)) |> display
plot(sol.k, sol.Sk) |> display
plot(sol.r, sol.gr, xlims=(0,10)) |> display

M = 2^14 - 1
ρ = [0.2, 0.3, 0.1]
kBT = 1.0
dims = 3

pot = OrnsteinZernike.MultiComponentHardSpheres([1.0, 0.6, 1.5]);
system = OrnsteinZernike.SimpleLiquid(dims, ρ, kBT, pot)
closure = OrnsteinZernike.PercusYevick();
method = OrnsteinZernike.Exact(M)
sol =  OrnsteinZernike.solve(system, closure, method);
using Plots
plot(sol.k, sol.Ck[:, 1, 1]) |> display
plot(sol.r, sol.Cr[:, 1, 1], xlims=(0,5)) |> display
plot(sol.k, sol.Sk[:, 1, 1], xlims=(0,20)) |> display
plot(sol.r, sol.gr[:, 1, 1], xlims=(0,5), ylims=(0,5)) |> display

# sol.gr[1200, 4, 6] == 1.463703395146971


M = 2^12 - 1

ρ = ones(10)*0.1
kBT = 1.0
dims = 3

pot = OrnsteinZernike.MultiComponentHardSpheres(ones(10));
system = OrnsteinZernike.SimpleLiquid(dims, ρ, kBT, pot)
closure = OrnsteinZernike.PercusYevick();
method = OrnsteinZernike.Exact(M)
@time sol = OrnsteinZernike.solve(system, closure, method);
using Plots
# plot(sol.k, sol.Ck[:, 1, 1]) |> display
# plot(sol.r, sol.Cr[:, 1, 1], xlims=(0,5)) |> display
plot(sol.k, sol.Sk[:, 1, 1], xlims=(0,20)) |> display
plot(sol.r, sol.gr[:, 1, 1], xlims=(0,10)) |> display


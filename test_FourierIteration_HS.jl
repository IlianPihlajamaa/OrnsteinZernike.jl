import Pkg; Pkg.activate(".")

import OrnsteinZernike
using Plots

# Single Component

M = 2^12
ρ = 0.9
kBT = 1.1
dims = 3

pot = OrnsteinZernike.SingleComponentHardSpheres()
system = OrnsteinZernike.SimpleLiquid(dims, ρ, kBT, pot)
closure = OrnsteinZernike.PercusYevick()
method = OrnsteinZernike.FourierIteration(M = M, tolerance=10^-10, mixing_parameter=0.8; dr=0.01)
@profview  sol = OrnsteinZernike.solve(system, closure, method)



p1 = plot(sol.k, sol.Ck, xlims=(0,30)) 
p2 = plot(sol.r, sol.Cr, xlims=(0,6)) 
p3 = plot(sol.k, sol.Sk, xlims=(0,30)) 
p4 = plot(sol.r, sol.gr, xlims=(0,6))  
sol = OrnsteinZernike.solve(system, closure, OrnsteinZernike.Exact(M))
plot!(p1, sol.k, sol.Ck) 
plot!(p2, sol.r, sol.Cr) 
plot!(p3, sol.k, sol.Sk) 
plot!(p4, sol.r, sol.gr) 
plot(p1, p2, p3, p4)

# multi Component
# M = 2^14 - 1
# ρ = [1.0, 0.8, 1.2]*0.2
# kBT = 1.1
# dims = 3

M = 2^14-1
ρ = [0.2, 0.3, 0.3]
kBT = 1.1
dims = 3

pot = OrnsteinZernike.MultiComponentHardSpheres([1.0, 1.2, 0.8])
system = OrnsteinZernike.SimpleLiquid(dims, ρ, kBT, pot)

closure = OrnsteinZernike.PercusYevick()
method = OrnsteinZernike.FourierIteration(M = M, tolerance=10^-10)
sol =  OrnsteinZernike.solve(system, closure, method)
p1 = plot(sol.k, sol.Ck[:, 2, 2], xlims=(0,30)) 
p2 = plot(sol.r, sol.Cr[:, 2, 3], xlims=(0,6)) 
p3 = plot(sol.k, sol.Sk[:, 2, 3], xlims=(0,90)) 
p4 = plot(sol.r, sol.gr[:, 2, 3], xlims=(0,6))  
sol = OrnsteinZernike.solve(system, closure, OrnsteinZernike.Exact(M))
plot!(p1, sol.k, sol.Ck[:, 2, 2]) 
plot!(p2, sol.r, sol.Cr[:, 2, 3]) 
plot!(p3, sol.k, sol.Sk[:, 2, 3]) 
plot!(p4, sol.r, sol.gr[:, 2, 3]) 
plot(p1, p2, p3, p4)
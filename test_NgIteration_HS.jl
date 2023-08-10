import Pkg; Pkg.activate(".")

import OrnsteinZernike


M = 2^18
ρ = 0.2
kBT = 1.0
dims = 3

pot = OrnsteinZernike.SingleComponentHardSpheres()
system = OrnsteinZernike.SimpleLiquid(dims, ρ, kBT, pot)
closure = OrnsteinZernike.PercusYevick()
method = OrnsteinZernike.NgIteration(tolerance=10^-10, N_stages=5, M=M)
sol1 = OrnsteinZernike.solve(system, closure, method)
    
using Plots
# p1 = plot(sol1.k, sol1.Ck, xlims=(0,30)) 
# p2 = plot(sol1.r, sol1.Cr, xlims=(0,6)) 
# p3 = plot(sol1.k, sol1.Sk, xlims=(0,30)) 
# p4 = plot(sol1.r, sol1.gr, xlims=(0,6)) 

closure = OrnsteinZernike.PercusYevick()
@time sol = OrnsteinZernike.solve(system, closure, OrnsteinZernike.Exact(M))
# plot!(p1, sol.k, sol.Ck) 
# plot!(p2, sol.r, sol.Cr) 
# plot!(p3, sol.k, sol.Sk) 
# plot!(p4, sol.r, sol.gr) 

# plot(p1, p2, p3, p4)

closure = OrnsteinZernike.PercusYevick()

rel_err = @. ( sol1.Cr - sol.Cr) / sol.Cr

plot(sol.r, rel_err, label = "gr_iter - gr_exact", xlims=(0.0, 2.0))



# MultiComponentHardSpheres
using Random
Random.seed!(3464)
M = 1000
ρ = rand(5)*0.6
kBT = 1.1
dims = 3

pot = OrnsteinZernike.MultiComponentHardSpheres(rand(5)*2)
system = OrnsteinZernike.SimpleLiquid(dims, ρ, kBT, pot)
closure = OrnsteinZernike.PercusYevick()
@time sol = OrnsteinZernike.solve(system, closure, OrnsteinZernike.Exact(M))

p1 = plot(sol.k, sol.Ck[:, 2, 2], xlims=(0,30)) 
p2 = plot(sol.r, sol.Cr[:, 2, 2], xlims=(0,6)) 
p3 = plot(sol.k, sol.Sk[:, 2, 2], xlims=(0,90)) 
p4 = plot(sol.r, sol.gr[:, 2, 3], xlims=(0,6))
method = OrnsteinZernike.NgIteration(M=M, tolerance=10^-10, N_stages=4)
sol = OrnsteinZernike.solve(system, closure, method)

plot!(p1, sol.k, sol.Ck[:, 2, 2]) 
plot!(p2, sol.r, sol.Cr[:, 2, 2]) 
plot!(p3, sol.k, sol.Sk[:, 2, 2]) 
plot!(p4, sol.r, sol.gr[:, 2, 3]) 
plot(p1, p2, p3, p4)


import Pkg; Pkg.activate(".")
using OrnsteinZernike,  Plots

# Make sure the discontinuity is a multiple of dr
Rmax = 10.0
Ms = 10 * round.(Int,  10 .^ (range(1,4,length=30)))

p1 = zeros(length(Ms))
p2 = zeros(length(Ms))
ρ = 0.189411447 * sqrt(2)
kBT = 1.0
dims = 3 

pot = HardSpheres(1.0)
system = SimpleLiquid(dims, ρ, kBT, pot)
for (i,M) in enumerate(Ms)
    @show i, M
    dr = Rmax/M
    method = NgIteration(M=M, dr=dr, verbose=false)
    sol1 = solve(system, PercusYevick(), method)
    p1[i] = compute_virial_pressure(sol1, system)/ρ/kBT-1.0
    @show p1[i]
    sol2 = solve(system, PercusYevick(), Exact(M=M, dr=dr))
    p2[i] = compute_virial_pressure(sol2, system)/ρ/kBT-1.0
end
η = ρ/6*π
@show pexact = ρ*kBT*(1+2η+3η^2)/(1-η)^2 /ρ/kBT-1.0
pl = scatter(Ms, abs.(p1.-pexact)./pexact, label="Iterative")
plot!(Ms, abs.(p2.-pexact)./pexact, label="from exact c(k)")
plot!(ylabel="log10(relative error)", xlabel="log10(M)", xscale=:log, yscale=:log)


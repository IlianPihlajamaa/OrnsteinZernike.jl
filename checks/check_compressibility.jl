import Pkg; Pkg.activate(".")
using OrnsteinZernike,  Plots

# Make sure the discontinuity is a multiple of dr
Rmax = 20.0
Ms = Rmax * round.(Int,  10 .^ (range(1,4,length=30)))

χ1 = zeros(length(Ms))
χ2 = zeros(length(Ms))
ρ = 1.0
kBT = 1.0
dims = 3 

pot = HardSpheres(1.0)
system = SimpleFluid(dims, ρ, kBT, pot)
for (i,M) in enumerate(Ms)
    @show i, M
    dr = Rmax/M
    method = NgIteration(M=M, dr=dr, verbose=false)
    sol1 = solve(system, PercusYevick(), method)
    χ1[i] = compute_compressibility(sol1, system)
    @show χ1[i]
    sol2 = solve(system, PercusYevick(), Exact(M=M, dr=dr))
    χ2[i] = compute_compressibility(sol2, system)
end
η = ρ/6*π
crint = ((-4 + η)* (2 + η^2))/(24(-1 + η)^4)
χ = 1/(1-ρ*4π*crint)/ρ/kBT

@show χexact = ((1+2η)^-2 * (1-η)^4)/ρ/kBT
pl = scatter(Ms, abs.(χ1.-χexact)./χexact, label="Iterative")
plot!(Ms, abs.(χ1.-χexact)./χexact, label="from exact c(k)")
plot!(ylabel="log10(relative error)", xlabel="log10(M)", xscale=:log, yscale=:log)


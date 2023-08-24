import Pkg; Pkg.activate(".")
using Revise
using OrnsteinZernike,  Plots
import Roots




N = 30
maxpower=5
# Make sure the discontinuity is a multiple of dr
Ms = 20*round.(Int,  10 .^ (range(1,maxpower,length=N)))

p1 = zeros(length(Ms))
p2 = zeros(length(Ms))
p3 = zeros(length(Ms))
ρ = 0.189411447 * sqrt(2)
kBT = 1.0
dims = 3 

pot = HardSpheres(1.0)
system = SimpleLiquid(dims, ρ, kBT, pot)
for (i,M) in enumerate(Ms)
    @show i, M
    dr = 20.0/M


    method = NgIteration(M=M, dr=dr, verbose=false)
    @time sol1 = solve(system, PercusYevick(), method)

    p1[i] = compute_virial_pressure(sol1, system)/ρ/kBT-1.0
    # p3[i] = 2/3 * pi * ρ * maximum(sol1.gr) 

    @time sol2 = solve(system, PercusYevick(), Exact(M=M, dr=dr))
    p2[i] = compute_virial_pressure(sol2, system)/ρ/kBT-1.0

end
η = ρ/6*π
@show pexact = ρ*kBT*(1+2η+3η^2)/(1-η)^2 /ρ/kBT-1.0
pl = scatter(log10.(Ms), log10.(abs.(p1.-pexact)./pexact), label="Iteration")
plot!(log10.(Ms), log10.(abs.(p2.-pexact)./pexact), ylims= (-10,0), label="Exact")
plot!(ylabel="log10(relative error)", xlabel="log10(M)")
pl |> display


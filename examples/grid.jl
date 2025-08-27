using OrnsteinZernike, Plots

system = SimpleFluid(3, 0.8, 1.0, LennardJones(1.0, 1.0))
closure = HypernettedChain()

sol_fine   = solve(system, closure; method=NgIteration(M=2000, dr=0.005))
sol_coarse = solve(system, closure; method=NgIteration(M=1000, dr=0.01))

plot(sol_fine.r, sol_fine.ck, label="c(r), fine grid")
plot!(sol_coarse.r, sol_coarse.ck, label="c(r), coarse grid", ls=:dash)
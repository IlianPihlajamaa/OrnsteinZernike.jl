# Accuracy

The Fourier tranforms used by the `FourierIteration`, and `NgIteration` solvers to represent the integrals using first order accuracy in $n$ dimensions and second order accuracy in three dimensions (with the trapezoidal rule). To obtain the latter accuracy, it is important that any discontinuities of the interaction potential lie on an exact multiple of `dr`. To test this, we can compute the pressure of a hard-sphere system, and compare to the exact value. Below, we compute the relative error for different values of the number of gridpoints `M`, and plot the result on a log-log-scale

```@example 
using OrnsteinZernike,  Plots

# Make sure the discontinuity is a multiple of dr
Rmax = 10.0
Ms = 10 * round.(Int,  10 .^ (range(1,4,length=30)))

p1 = zeros(length(Ms))
ρ = 0.189411447 * sqrt(2)
kBT = 1.0
dims = 3 

pot = HardSpheres(1.0)
system = SimpleLiquid(dims, ρ, kBT, pot)
for (i,M) in enumerate(Ms)
    dr = Rmax/M
    method = NgIteration(M=M, dr=dr, verbose=false)
    sol1 = solve(system, PercusYevick(), method)
    p1[i] = compute_virial_pressure(sol1, system)/ρ/kBT-1.0
end
η = ρ/6*π
pexact = ρ*kBT*(1+2η+3η^2)/(1-η)^2 /ρ/kBT-1.0
scatter(Ms, abs.(p1.-pexact)./pexact)
plot!(ylabel="relative error", xlabel="M", xscale=:log, yscale=:log)
```
We can see that the method has well-behaved second order convergence, and that with $M=10^4$, we get almost 6 digits of relative accuracy.
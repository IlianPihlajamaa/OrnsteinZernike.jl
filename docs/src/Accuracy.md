# Accuracy

The discrete Fourier tranforms used by the `FourierIteration`, and `NgIteration` solvers to represent their continuous counterparts using first order accuracy in $n$ dimensions and second order accuracy in three dimensions (with the trapezoidal rule). To obtain the latter accuracy, it is important that any discontinuities of the interaction potential lie on an exact multiple of `dr`. To test this, we can compute the pressure of a hard-sphere system, and compare to the exact value. Below, we compute the relative error for different values of the number of gridpoints `M`, and plot the result on a log-log-scale

```julia
using OrnsteinZernike,  Plots

# Make sure the discontinuity is a multiple of dr
Rmax = 10.0
M_array = 10 * round.(Int,  10 .^ (range(1,4,length=20)))
p = zeros(length(M_array))
ρ = 0.3
kBT = 1.0
dims = 3 
pot = HardSpheres(1.0)
system = SimpleFluid(dims, ρ, kBT, pot)

for (i,M) in enumerate(M_array)
    dr = Rmax/M
    method = NgIteration(M=M, dr=dr, verbose=false)
    sol = solve(system, PercusYevick(), method)
    pressure = compute_virial_pressure(sol, system)
    p[i] = pressure/ρ/kBT-1.0
    println("The pressure = ", round(pressure, digits=8), " with M = $(M) gridpoints.")
end
```
We can see that the method has well-behaved second order convergence, and that with $M=10^4$, we get almost 6 digits of relative accuracy.

```julia 5
η = ρ/6*π
p_exact = (1+2η+3η^2)/(1-η)^2-1.0
scatter(M_array, abs.(p.-p_exact)./p_exact)
plot!(ylabel="relative error", xlabel="M", xscale=:log, yscale=:log)
```

In principle, these results can be extrapolated to improve the accuracy further.
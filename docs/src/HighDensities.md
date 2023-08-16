# Solving at high densities

Sometimes, especially when dealing with high densities, the solvers do not converge out-of-the-box. Typically, the easiest way to solve this is to initialize the iterative procedures with solutions of the same equation at a lower density. While this can be done by hand using the `init` keyword argument of the `solve` function, a convenient method `DensityRamp` is implemented to do this automatically. For example, consider a hard sphere system at `ρ = 1.2`.

```julia
using OrnsteinZernike
ρ = 1.2
kBT = 1.0
dims = 3

pot = HardSpheres(1.0)
system = SimpleLiquid(dims, ρ, kBT, pot)
closure = PercusYevick()
method = NgIteration()
sol = solve(system, closure, method);
```
```
After iteration 0, the error is 114.96374587156.
After iteration 10, the error is 246.47358557791.
After iteration 20, the error is 48.60970316414.
[...]
After iteration 980, the error is 600.93302065646.
After iteration 990, the error is 346.00914196416.
After iteration 1000, the error is 78.17527233114.
ERROR: Recursive iteration did not converge within 1001 steps. Current error = 78.17527233113674.
 [1] error(s::String)
   @ Base .\error.jl:35
[...]
```

To solve an issue like this, we can define `DensityRamp` method. This takes two arguments. The first is the actual method by which to solve the equations, which we leave as `NgIteration()` and the second is a list of densities to be evaluated before the target density.

```
densities = [0.7, 0.9, 1.0, 1.1]
method2 = DensityRamp(method, densities)
sol = @time solve(system, closure, method2);
```

```
Solving the system at ρ = 0.7.

After iteration 0, the error is 3.7927462113.
After iteration 10, the error is 0.02748002125.
After iteration 20, the error is 2.459e-8.
Converged after 25 iterations, the error is 8.0e-11.

Solving the system at ρ = 0.9.

After iteration 0, the error is 1.46229405321.
After iteration 10, the error is 0.02624986205.
After iteration 20, the error is 1.78285e-6.
Converged after 29 iterations, the error is 3.0e-11.

Solving the system at ρ = 1.0.

After iteration 0, the error is 6.91022217054.
After iteration 10, the error is 0.25666355675.
After iteration 20, the error is 0.01233410847.
After iteration 30, the error is 1.888465e-5.
After iteration 40, the error is 1.94e-9.
Converged after 44 iterations, the error is 8.0e-11.

Solving the system at ρ = 1.1.

After iteration 0, the error is 42.29810538999.
After iteration 10, the error is 0.37566377791.
After iteration 20, the error is 0.27554017425.
[...].
After iteration 70, the error is 1.555e-8.
After iteration 80, the error is 5.2e-10.
Converged after 84 iterations, the error is 5.0e-11.

Solving the system at ρ = 1.2.

After iteration 0, the error is 57.02991026895.
After iteration 10, the error is 0.36390031646.
After iteration 20, the error is 0.38501430163.
[...]
After iteration 110, the error is 2.706384e-5.
After iteration 120, the error is 7.97e-7.
After iteration 130, the error is 4.98e-9.
Converged after 135 iterations, the error is 7.0e-11.
  0.048063 seconds (1.86 k allocations: 1.330 MiB)
```

Let's plot and compare with the exact solution:

```@example dens
using OrnsteinZernike # hide
ρ = 1.2 # hide
kBT = 1.0 # hide
dims = 3 # hide
pot = HardSpheres(1.0) # hide
system = SimpleLiquid(dims, ρ, kBT, pot) # hide
closure = PercusYevick() # hide
method = NgIteration() # hide
densities = [0.7, 0.9, 1.0, 1.1] # hide
method2 = DensityRamp(method, densities) # hide
sol = solve(system, closure, method2); # hide
using Plots
plot(sol.r, sol.gr, xlims=(0,5), xlabel="r", ylabel="g(r)", lw=4, label="iterative")

method3 = Exact()
sol = solve(system, closure, method3);
plot!(sol.r, sol.gr, lw=2, color=:black, label="exact")
```

We can see that while the density is so high that $g(r)$ shows a non-physical negative dip, it converged to the right mathematical solution. Note that this method will not always lead to convergence. For example, setting the density to `ρ = 2.0` in this example, probably won't be correctly solved no matter how many intermediate densities are considered.
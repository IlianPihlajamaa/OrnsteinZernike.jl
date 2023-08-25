# Solving at high densities

Sometimes, especially when dealing with high densities, the solvers do not converge out-of-the-box. It is sometimes necessary to play with the solver settings to get it to converge
```julia
using OrnsteinZernike
ρ = 1.5
kBT = 1.0
dims = 3

pot = HardSpheres(1.0)
system = SimpleLiquid(dims, ρ, kBT, pot)
closure = PercusYevick()
method = NgIteration()
sol = solve(system, closure, method);
```
```
After iteration 0, the error is 3.41302002652.
After iteration 10, the error is 1.98231818283.
After iteration 20, the error is 1.54425698984.
After iteration 30, the error is 1.31506989937.
[...]
After iteration 960, the error is 0.12249150013.
After iteration 970, the error is 0.1162956682.
After iteration 980, the error is 0.12306597032.
After iteration 990, the error is 0.10235109774.
After iteration 1000, the error is 0.09999300441.
ERROR: Recursive iteration did not converge within 1001 steps. Current error = 0.09999300441413587.
```

Instead, we can try to change some solver settings:
```@example pyhd
using OrnsteinZernike
ρ = 1.5 # hide
kBT = 1.0 # hide
dims = 3 # hide
pot = HardSpheres(1.0) # hide
system = SimpleLiquid(dims, ρ, kBT, pot) # hide
closure = PercusYevick() # hide

M = 10^4
method = NgIteration(M = M; dr = 10/M, max_iterations=10^4, N_stages=4)
sol = solve(system, closure, method);

using Plots
plot(sol.r, sol.gr, xlims=(0,5), xlabel="r", ylabel="g(r)", lw=4, label="iterative")
sol = solve(system, closure, Exact());
plot!(sol.r, sol.gr, lw=2, color=:black, label="exact")
```
At these unphysically high densities the result is not very accurate. 




Sometimes convergence is accelerated if the solution of the same equations at a slighlty lowe density is used as an initial condition. This is especially useful if the solution is needed at many different densities While this can be done by hand using the `init` keyword argument of the `solve` function, a convenient method `DensityRamp` is implemented to do this automatically. For example, consider a hard sphere system at `ρ = 1.2`.

```@example ramp
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

To solve an issue like this, we can define `DensityRamp` method. This takes two arguments. The first is the actual method by which to solve the equations, which we leave as `NgIteration()` and the second is a list of densities to be evaluated before the target density. It returns a `Vector` of solution objects.

```@example ramp
densities = [1.15, 1.18]
method = NgIteration()
method2 = DensityRamp(method, densities)
sol = @time solve(system, closure, method2)[end];
```

```@example ramp
using Plots
plot(sol.r, sol.gr, xlims=(0,5), xlabel="r", ylabel="g(r)", lw=4, label="iterative")
method3 = Exact()
sol = solve(system, closure, method3);
plot!(sol.r, sol.gr, lw=2, color=:black, label="exact")
```


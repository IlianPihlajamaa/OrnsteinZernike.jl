# Solving at high densities

Sometimes, especially when dealing with high densities, the solvers do not converge out-of-the-box. It may be necessary to play with the solver settings to get it to converge
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

Instead, we can try to change some solver settings (see [`NgIteration`](@ref) for detailed descriptions):

```@example pyhd
using OrnsteinZernike
ρ = 1.5 # hide
kBT = 1.0 # hide
dims = 3 # hide
pot = HardSpheres(1.0) # hide
system = SimpleLiquid(dims, ρ, kBT, pot) # hide
closure = PercusYevick() # hide

M = 10^4 # number of gridpoints
dr = 100.0/M # grid spacing
max_iterations = 10^4 # max number of iterations before convergence 
N_stages = 8 # number of previous iterations to use for the next guess

method = NgIteration(M=M; dr=dr, max_iterations=max_iterations, N_stages=N_stages)
sol = solve(system, closure, method);
```
```@example pyhd
using Plots
plot(sol.r, sol.gr, xlims=(0,5), xlabel="r", ylabel="g(r)", lw=4, label="iterative")
sol2 = solve(system, closure, Exact(M=M; dr=dr));
plot!(sol2.r, sol2.gr, lw=2, color=:black, label="exact")
```
Which converges even though the density is clearly nonphysically high.

Sometimes convergence is accelerated if the solution of the same equations at a slightly lower density is used as an initial condition. This is especially useful if the solution is needed at many different densities. While this can be done by hand using the `init` keyword argument of the `solve` function, a convenient method `DensityRamp` is implemented to do this automatically. For example, consider a hard sphere system at `ρ = 1.2`.

```julia
using OrnsteinZernike
ρ = 1.2
kBT = 1.0
dims = 3

pot = HardSpheres(1.0)
system = SimpleLiquid(dims, ρ, kBT, pot)
closure = PercusYevick()
method = NgIteration(M=5000, dr=0.01)
sol = solve(system, closure, method);
```
```
After iteration 0, the error is 2.51773674949.
After iteration 10, the error is 1.41253963056.
After iteration 20, the error is 1.2605856682.
After iteration 30, the error is 0.70489573233.
After iteration 40, the error is 0.44336537233.
After iteration 50, the error is 0.39947202862.
After iteration 60, the error is 0.25827952962.
After iteration 70, the error is 0.23273727229.
After iteration 80, the error is 0.13551597572.
After iteration 90, the error is 0.13460620172.
After iteration 100, the error is 0.29303789746.
After iteration 110, the error is 0.12428896865.
After iteration 120, the error is 0.13350846977.
After iteration 130, the error is 0.00388920336.
After iteration 140, the error is 0.00018271462.
After iteration 150, the error is 7.94778e-6.
After iteration 160, the error is 1.27727e-6.
After iteration 170, the error is 5.4383e-7.
After iteration 180, the error is 1.78e-9.
Converged after 185 iterations, the error is 7.0e-11.
```


To do this, we can use the `DensityRamp` solver. This takes two arguments. The first is the actual method by which to solve the equations, which we leave as `NgIteration()` and the second is a list of densities to be evaluated before the target density. It returns a `Vector` of solution objects.

```julia
densities = [1.15, 1.18]
method2 = DensityRamp(method, densities)
sol = @time solve(system, closure, method2);
```
```
Solving the system at ρ = 1.15.

After iteration 0, the error is 2.45107333815.
After iteration 10, the error is 1.14683100707.
After iteration 20, the error is 0.56643550932.
After iteration 30, the error is 0.49252243775.
After iteration 40, the error is 0.40727608832.
After iteration 50, the error is 0.03465798498.
After iteration 60, the error is 0.00212896547.
After iteration 70, the error is 2.296428e-5.
After iteration 80, the error is 4.4255e-7.
After iteration 90, the error is 1.634e-8.
Converged after 98 iterations, the error is 3.0e-11.

Solving the system at ρ = 1.18.

After iteration 0, the error is 234.73047273369.
After iteration 10, the error is 0.09622279909.
After iteration 20, the error is 0.11509610499.
After iteration 30, the error is 0.03246309423.
After iteration 40, the error is 0.01668981609.
After iteration 50, the error is 0.00030500371.
After iteration 60, the error is 1.38937e-6.
After iteration 70, the error is 1.1985e-7.
After iteration 80, the error is 1.516e-8.
After iteration 90, the error is 9.3e-10.
After iteration 100, the error is 3.0e-10.
Converged after 104 iterations, the error is 9.0e-11.

Solving the system at ρ = 1.2.

After iteration 0, the error is 193.67155095726.
After iteration 10, the error is 0.06503380948.
After iteration 20, the error is 0.11540431925.
After iteration 30, the error is 0.0836301059.
After iteration 40, the error is 0.00063923863.
After iteration 50, the error is 5.98769e-6.
After iteration 60, the error is 4.2188e-7.
After iteration 70, the error is 3.1e-9.
Converged after 80 iterations, the error is 9.0e-11.
  0.071976 seconds (17.20 k allocations: 4.204 MiB)
```

While, in total, this has increased the number of iterations done, we have obtained the solution at two different densities as well.
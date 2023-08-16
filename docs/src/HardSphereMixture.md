# Hard-sphere mixture

The previous example showed the case of a 1-component system. Let's instead look at mixtures here. Consider a 3:1 mixture of hard spheres with sizes 0.5 and 1.0. We solve the system with the same three steps as before.

First, we define the potential. Here, we must use `HardSpheres`, which takes a vector containing the diameters of each species.

```@example hs
using OrnsteinZernike
D = [0.5, 1.0]
potential = HardSpheres(D)
```

Secondly, we define the system. In this example, the total number density is $\rho = 1.6$. For mixtures, the system expects a vector of individual densities. Those are computed by multiplying the total density with the concentration fraction vector.

```@example hs
ρ_total = 1.6
ρ = ρ_total*[0.75, 0.25] ## it is a 3:1 system
dims = 3 # we consider a 3D system
kBT = 1.0 # thermal energy
system = SimpleLiquid(dims, ρ, kBT, potential)
```

Thirdly, we define the closure
```@example hs
closure = PercusYevick()
```

And now we solve the system. 
```@example hs
sol = solve(system, closure)
```

For mixtures, the fields `sol.gr`, `sol.cr`, `sol.ck`, and `sol.Sk` are now three dimensional arrays with shape `(Nr, Ns, Ns)`. For example, $g_{12}(r_6)$ is stored in `sol.gr[6,1,2]`.

We just solved the system using the default iterative solver introduced by Ng [`NgIteration`](@ref). However, in this specific case, an exact solution is implemented. To use this, we specify the method `Exact()`.

```@example hs
method = Exact()
sol_exact = solve(system, closure, method)
```

Let's plot the resulting $g(r)$. 

```@example hs
using Plots
plot(sol.r, sol.gr[:, 1, 1], xlims=(0,5), xlabel="r", ylabel="g(r)", lw=4, label="g11(r) iterative")
plot!(sol.r, sol.gr[:, 1, 2], xlabel="r", ylabel="g(r)", lw=4, label="g12(r) iterative")
plot!(sol.r, sol.gr[:, 2, 2], xlabel="r", ylabel="g(r)", lw=4, label="g22(r) iterative")
plot!(sol_exact.r, sol_exact.gr[:, 1, 1], lw=2, c=:black, label="exact")
plot!(sol_exact.r, sol_exact.gr[:, 1, 2], lw=2, c=:black, label=nothing)
plot!(sol_exact.r, sol_exact.gr[:, 2, 2], lw=2, c=:black, label=nothing)
```

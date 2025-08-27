# First steps

In order to solve an Ornstein-Zernike equation, there are a number of things that need to be specified. 
1. The interaction potential, 
2. The system parameters: density, temperature and dimensionality
3. The closure relation
4. The solver

Let's start with a very simple example: a three-dimensional 1-component system, where the particles interact according to an inverse power law potential $u(r)=\epsilon (\sigma/r)^n$. In this case, we can make use of the built-in potential 
```@example lj
using OrnsteinZernike
ϵ = 1.0
σ = 1.0
n = 8
potential = PowerLaw(ϵ, σ, n)
```

Now that we have the potential, we define the system

```@example lj
dims = 3 # we consider a 3D system
ρ = 0.6 # number density
kBT = 1.0 # thermal energy
system = SimpleFluid(dims, ρ, kBT, potential)
```

The `SimpleFluid` object is meant to be used when dealing with systems that have spherically symmetric interaction potentials and no external fields. 

The third step is to define a closure relation. For now, let's stick to the simple Hypernetted Chain closure
```@example lj
closure = HypernettedChain()
```
We can now solve the system. 
```@example lj
sol = solve(system, closure)
```

Note that we have not specified the method by which we do so. In this case, a simple default method is chosen that works well in most cases.

The `sol` object contains fields for the radial distribution function `gr`, direct correlation function (in real and Fourier space) `cr` and `ck`, the static structure factor `Sk`, and arrays for `r` and `k`. For example, we can plot the radial distribution function as follows:

```@example lj
using Plots
plot(sol.r, sol.gr, xlims=(0,5), xlabel="r", ylabel="g(r)")
```


Full code:
```julia
using OrnsteinZernike
ϵ = 1.0
σ = 1.0
n = 8
potential = PowerLaw(ϵ, σ, n)
dims = 3 # we consider a 3D system
ρ = 0.6 # number density
kBT = 1.0 # thermal energy
system = SimpleFluid(dims, ρ, kBT, potential)
closure = HypernettedChain()
sol = solve(system, closure)
using Plots
plot(sol.r, sol.gr, xlims=(0,5), xlabel="r", ylabel="g(r)")
```

# Defining your own potentials

Creating your own potentials type is very easy. We can do this by making use of the `CustomPotential`

### Example 

We want to implement the interaction potential 
$$u(r) = \epsilon \left(\frac{\sigma}{r}\right)^6.$$

To do this, first we define the function. This function should take two arguments, `r::Number` and `p`, the latter of which contains optional parameters that the function uses. For example:

```@example 1
function my_pot(r, p)
    return p.ϵ * (p.σ / r)^6
end
```

Now we can instantiate the `CustomPotential`, using e.g. a `NamedTuple` to pass in the parameters:

```@example 1
using OrnsteinZernike
p = (ϵ = 1.0, σ = 1.0)
potential = CustomPotential(u, p)
```

And use the potential as any other 

```@example 1
dims = 3 # we consider a 3D system
ρ = 0.6 # number density
kBT = 1.0 # thermal energy
system = SimpleLiquid(dims, ρ, kBT, potential)
closure = HypernettedChain()
sol = solve(system, closure)
using Plots
plot(sol.r, sol.gr, xlims=(0,5), xlabel="r", ylabel="g(r)")
```

## Mixtures

In the case of multicomponent systems, instead of a number the function should return a `StaticMatrix` from the `StaticArrays` package containing values for $u_{ij}$. 

### Example

Suppose we want to implement the same potential for the multicomponent case:
$$u_{ij}(r) = \epsilon_{ij} (\frac{\sigma_{ij}}{r_{ij}})^6.$$

While one could implement this in one line with broadcasting, here the function is written out fully for clarity:

```@example 2
using OrnsteinZernike, StaticArrays 

function mypotential(r, p)
    # we can construct a mutable sized matrix first
    Nspecies = size(p.ϵ, 1)
    out = MMatrix{Nspecies, Nspecies, Float64, Nspecies*Nspecies}(undef) 
    for species2 = 1:Nspecies
        for species1 = 1:Nspecies
            out[species1, species2] = p.ϵ[species1, species2] * (p.σ[species1, species2] / r) ^ 6
        end
    end
    # and convert it to an immutable variant
    return SMatrix(out) 
end
```

and now we can use it:

```@example 2
ϵ = SMatrix{2,2}([1.0 2.0; 0.4 0.9])
σ = SMatrix{2,2}([1.0 1.0; 1.0 0.8])
p = (ϵ = ϵ, σ = σ)
dims = 3 # we consider a 3D system
potential = CustomPotential(mypotential, p)
ρ = [0.25, 0.25] # number density
kBT = 1.0 # thermal energy
system = SimpleLiquid(dims, ρ, kBT, potential)
closure = HypernettedChain()
sol = solve(system, closure)
using Plots
plot(sol.r, sol.gr[:, 1, 1], xlims=(0,5), xlabel="r", ylabel="g(r)", label="g11(r)")
plot!(sol.r, sol.gr[:, 1, 2], xlims=(0,5), xlabel="r", ylabel="g(r)", label="g12(r)")
plot!(sol.r, sol.gr[:, 2, 2], xlims=(0,5), xlabel="r", ylabel="g(r)", label="g22(r)")
```

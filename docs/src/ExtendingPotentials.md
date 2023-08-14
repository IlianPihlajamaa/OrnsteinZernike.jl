# Defining your own potentials

Creating your own potentials type is very easy. It takes two steps. First, the type itself must be created, and secondly, function that evaluates the potentials must be overloaded. 

### Example 

We want to implement the interaction potential 
$$u(r) = \epsilon (\frac{\sigma}{r})^6.$$

To do this, first we define the type. Note that the new potential must be made a subtype of OrnsteinZernike.Potential
```@example 1
using OrnsteinZernike

import OrnsteinZernike.Potential
struct MyPot <: Potential 
    epsilon::Float64
    sigma::Float64
end
```

Now we can define how this potential should be evaluated. We must overload   `evaluate_potential(potential::Potential, r::Number)`


```@example 1
import OrnsteinZernike.evaluate_potential
function OrnsteinZernike.evaluate_potential(pot::MyPot, r::Number)
    return pot.epsilon * (pot.sigma / r) ^ 6
end
```

Now we can use the potential as any other 

```@example 1
ϵ = 1.0
σ = 1.0
potential = MyPot(ϵ, σ)
dims = 3 # we consider a 3D system
ρ = 0.5 # number density
kBT = 1.0 # thermal energy
system = SimpleLiquid(dims, ρ, kBT, potential)
closure = HypernettedChain()
sol = solve(system, closure)
using Plots
plot(sol.r, sol.gr, xlims=(0,5), xlabel="r", ylabel="g(r)")
```

## Mixtures

In the case of multicomponent systems, instead of a number the `evaluate_potential` should return a `StaticMatrix` from the `StaticArrays` package containing either values for $u_{ij}$. 

### Example

Suppose we want to implement the same potential for the multicomponent case:
$$u_{ij}(r) = \epsilon_{ij} (\frac{\sigma_{ij}}{r_{ij}})^6.$$


```@example 2
using OrnsteinZernike, StaticArrays

import OrnsteinZernike.Potential
struct MyPot2{Nspecies} <: Potential 
    epsilon::Matrix{Float64}
    sigma::Matrix{Float64}
end
```
Here we have made the type parametric with respect to the number of species

```@example 2
import OrnsteinZernike.evaluate_potential
function OrnsteinZernike.evaluate_potential(pot::MyPot2{Nspecies}, r::Number) where Nspecies
    # we can construct a mutable sized matrix first
    out = MMatrix{Nspecies, Nspecies, Float64, Nspecies*Nspecies}(undef) 
    for species2 = 1:Nspecies
        for species1 = 1:Nspecies
            out[species1, species2] = pot.epsilon[species1, species2] * (pot.sigma[species1, species2] / r) ^ 6
        end
    end
    # and convert it to an immutable variant
    return SMatrix(out) 
end
```

and now we can use it:

```@example 3
using OrnsteinZernike, StaticArrays # hide

import OrnsteinZernike.Potential # hide
struct MyPot2{Nspecies} <: Potential  # hide
    epsilon::Matrix{Float64} # hide
    sigma::Matrix{Float64} # hide
end # hide

import OrnsteinZernike.evaluate_potential # hide
function OrnsteinZernike.evaluate_potential(pot::MyPot2{Nspecies}, r::Number) where Nspecies # hide
    # we can construct a mutable sized matrix first # hide
    out = MMatrix{Nspecies, Nspecies, Float64, Nspecies*Nspecies}(undef)  # hide
    for species2 = 1:Nspecies # hide
        for species1 = 1:Nspecies # hide
            out[species1, species2] = pot.epsilon[species1, species2] * (pot.sigma[species1, species2] / r) ^ 6 # hide
        end # hide
    end # hide
    # and convert it to an immutable variant # hide
    return SMatrix(out)  # hide
end # hide
ϵ = [1.0 2.0; 0.4 0.9]
σ = [1.0 1.0; 1.0 0.8]
dims = 3 # we consider a 3D system
potential = MyPot2{2}(ϵ, σ)
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

# Defining your own closure

Creating your own closure type is very easy. It takes two steps. First, the type itself must be created, and secondly, one of the methods that evaluates the closure must be overloaded.

There are three options. One can overload either
1.  `bridge_function(closure, r, mayer_f, γ)`
2.  `closure_c_from_gamma(closure, r, mayer_f, γ, βu_LR)`
3.  `closure_cmulr_from_gammamulr(closure, r, mayer_f, γ, βu_LR)`

Here, `mayer_f` is the Mayer-f function $f(r) = \exp(-\beta u(r)) - 1$, if the closure needs that, and `βu_LR` is the long range part of the potential. In practise, which of the three methods is overloaded can be arbitrary and depends on what is most convenient.

### Example 

Assume that we have forgotten that the HypernettedChain closure is already implemented, and we wanted to reimplement it. The hypernetted chain closure approximates $c(r) \approx (f(r)+1)\exp(\gamma(r)) - \gamma(r) - 1$, or equivalently $b(r) \approx 0$.

First we define the type. Note that the new closure must be made a subtype of OrnsteinZernike.Closure
```@example 1
using OrnsteinZernike

import OrnsteinZernike.Closure
struct MyHNC <: Closure end
```

Now we can define how this closure should be evaluated. Here, we can for example either do

```@example 1
import OrnsteinZernike.bridge_function
function OrnsteinZernike.bridge_function(::MyHNC, _, _, _)
    return 0.0
end
```

or

```@example 1
import OrnsteinZernike.closure_c_from_gamma
function OrnsteinZernike.closure_c_from_gamma(::MyHNC, _, mayer_f, γ, _)
    return (mayer_f + 1) * exp(γ) - γ - 1
end
```

Now we can use the closure as any other 

```@example 1
ϵ = 1.0
σ = 1.0
n = 12
potential = PowerLaw(ϵ, σ, n)
dims = 3 
ρ = 0.6 
kBT = 1.0
system = SimpleFluid(dims, ρ, kBT, potential)
closure = MyHNC()
sol = solve(system, closure)
using Plots
plot(sol.r, sol.gr, xlims=(0,5), xlabel="r", ylabel="g(r)")
```

which can be compared to that of [First steps](@ref).

## Mixtures

In the case of multicomponent systems, instead of a number the function that is overloaded should return a `StaticMatrix` containing either values for $c_{ij}$ or $b_{ij}$. Since, in that case also the inputs are `StaticMatrix`, we can make use of `Julia`'s broadcasting syntax to perform an closure elementwise. 

```julia
import OrnsteinZernike.closure_c_from_gamma
function OrnsteinZernike.closure_c_from_gamma(::MyHNC, _, mayer_f, γ, _)
    return @. (mayer_f + 1) * exp(γ) - γ - 1 # note the @.
end
```
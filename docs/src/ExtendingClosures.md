# Defining your own closure

Creating a closure amounts to two small steps. First, define a subtype of
`OrnsteinZernike.Closure`. Second, implement how the bridge function should be
evaluated through

1. `bridge_function(closure, r, mayer_f, γ)` – the standard entry point.
2. (Optional) `closure_cmulr_point(closure, r, mayer_f, γ_SR, βu, βu_LR_disp,
   βu_LR_coul, q, uses_renorm)` – override this only if you need full control
   over the numerical formula.

The default evaluator calls `bridge_function`, takes care of Coulomb pieces, and
handles the dispersion tail when the closure reports that it expects the
renormalised quantity \(γ^* = γ_{SR} - βu_{LR}^{disp}\). This expectation is
communicated via the trait

```julia
uses_renormalized_gamma(::Closure) = false
```

and can be specialised for individual closures whenever needed.

```julia
import OrnsteinZernike: uses_renormalized_gamma
uses_renormalized_gamma(::MyClosure) = true
```

### Example 

Assume that we have forgotten that the HypernettedChain closure is already implemented, and we wanted to reimplement it. The hypernetted chain closure approximates $c(r) \approx (f(r)+1)\exp(\gamma(r)) - \gamma(r) - 1$, or equivalently $b(r) \approx 0$.

First we define the type. Note that the new closure must be made a subtype of OrnsteinZernike.Closure
```@example 1
using OrnsteinZernike

import OrnsteinZernike.Closure
struct MyHNC <: Closure end
```

The hypernetted-chain closure is recovered by saying that the bridge function is
zero everywhere:

```@example 1
import OrnsteinZernike.bridge_function
function OrnsteinZernike.bridge_function(::MyHNC, _, _, _)
    return 0.0
end
```

No explicit trait override is required because the closure does not rely on the
renormalised \(γ^*\).

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

In the case of multicomponent systems the bridge function should return a
`StaticMatrix` of the pair values. Because the arguments are themselves
`StaticMatrix` objects you can rely on broadcasting to write the expressions
elementwise.

```julia
import OrnsteinZernike.bridge_function
function OrnsteinZernike.bridge_function(::MyHNC, _, _, γ)
    return zero(γ) # or any other elementwise expression, e.g. @. ...
end
```

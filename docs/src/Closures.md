# Closures

This package defines several closure relations that can be used out of the box. It is straightforward to implement your own closure, see [Defining your own closure](@ref). 

Some closures, for example the Rogers-Young closure, include free parameters that may be fixed by the requirement of thermodynamic consistency. See the [Thermodynamic Consistency](@ref) page for an example. 

Some closures use a renormalized indirect correlation function $\gamma^*(r) = \gamma(r) - u_{LR}(r)$ instead of the standard one. Here, $u_{LR}(r)$ is the long-range tail of the interaction potential. The solver now makes this choice explicit: each closure advertises whether it wants $\gamma^*$ via the trait `uses_renormalized_gamma(closure)`, and each potential exposes its dispersion tail through `dispersion_tail(potential, kBT, r, βu)`. When a closure opts in, the appropriate tail is subtracted automatically before the bridge function is evaluated. If a renormalised-γ closure is paired with a potential that does not expose such a tail, the solver emits an error suggesting to wrap the potential in `WCADivision(...)` (Weeks–Chandler–Andersen split) or `AllShortRangeDivision(...)` (explicitly mark the interaction as fully short ranged). 


### Dispersion handling in practice

When you pick a closure for which `uses_renormalized_gamma` returns `true` (the
SMSA/HMSA family, Lee, Duh–Haymet, etc.) you must present the solver with a
potential that splits into short- and long-ranged pieces. You can do this in two
ways:

1. Wrap the base interaction in `WCADivision(potential, r_c)` to obtain the
   traditional Weeks–Chandler–Andersen short-range/dispersion tail separation.
2. Mark it as fully short ranged via `AllShortRangeDivision(potential)` so the solver knows that
   the dispersion tail is zero.

Failing to do so raises an error before the iteration starts. For example:

```julia
using OrnsteinZernike
closure = ZerahHansen()              # expects γ* = γ - βuₗᵣ
lj       = LennardJones(1.0, 1.0)
fluid    = SimpleFluid(3, 0.4, 1.0, lj)

# this throws: closure needs a dispersion split
# solve(fluid, closure, NgIteration())

# add a WCA split and the run succeeds
wca_fluid = SimpleFluid(3, 0.4, 1.0, WCADivision(lj, 2^(1/6)))
sol, = solve(wca_fluid, closure, NgIteration())
```

For charged systems, the charges are always dealt with analytically. They are not included in the dispersion tails.

## Implemented Closures
Below is an alphabetical list of implemented closures. We use the notation shown in the [Theory](@ref) section.

```@index
Modules = [OrnsteinZernike]
Pages = ["Closures.md"]
Order   = [:type]
Private = false
```

```@autodocs
Modules = [OrnsteinZernike]
Filter = t -> ((typeof(t) === DataType || typeof(t) === UnionAll) && t <: OrnsteinZernike.Closure)
Private = false
```

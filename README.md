# OrnsteinZernike.jl
A generic solver for Ornstein-Zernike equations from liquid state theory

[![Build status (Github Actions)](https://github.com/IlianPihlajamaa/OrnsteinZernike.jl/workflows/CI/badge.svg)](https://github.com/IlianPihlajamaa/OrnsteinZernike.jl/actions)
[![codecov.io](http://codecov.io/github/IlianPihlajamaa/OrnsteinZernike.jl/coverage.svg?branch=main)](http://codecov.io/github/IlianPihlajamaa/OrnsteinZernike.jl?branch=main)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://IlianPihlajamaa.github.io/OrnsteinZernike.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://IlianPihlajamaa.github.io/OrnsteinZernike.jl/dev)

One of the big triumphs in liquid state theory is the ability to approximate the structure of a liquid from the way its constituent particles interact. 
This can be done using the exact Ornstein-Zernike equation: 
$$g(r) - 1 = c(r) + \rho \int d\textbf{r}' c(\textbf{r}')(g(\textbf{r}- \textbf{r}') - 1)$$

Here $g(r)$ is the radial distribution function, which describes the structure of a liquid, $\rho$ is the number density, and $c(r)$ is the direct correlation function. Together with an approximate closure relation, which links $c(r)$ to the interaction potential of the particles, this integral equation can be solved. 

This package implements common solution methods for the Ornstein Zernike equation, for single component systems as well as mixtures. It implements many predefined closure relations and many interaction potentials. Moreover, extending the package to include other closure relations or interaction potentials takes minimal effort.

## Example use

Let's solve the Ornstein-Zernike equation for a single component three-dimensional system of hard spheres at number density $ρ = 0.5$. 

```julia
using OrnsteinZernike
dims = 3; kBT = 1.0; ρ = 0.5;
potential = SingleComponentHardSpheres()
system = SimpleLiquid(dims, ρ, kBT, potential)
closure = PercusYevick()
sol = @time solve(system, closure);
```
which prints
```
After iteration 0, the error is 10.4396417.
After iteration 10, the error is 1.2794593.
After iteration 20, the error is 1.1e-6.
Converged after 22 iterations, the error is 2.0e-7.
0.006574 seconds (176 allocations: 254.781 KiB)
```
Now that we have solved the equation, we can plot the solution:
```julia
using Plots
plot(sol.r, sol.gr, xlims=(0,5))
```
![image](docs/src/Figs/example.png)

See the <a href="https://ilianpihlajamaa.github.io/OrnsteinZernike.jl/dev/">Documentation</a> for more details.

Please open an issue if anything is unclear in the documentation, if any unexpected errors arise or for feature requests. PRs are of course also welcome.

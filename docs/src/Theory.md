# Theory

In this section, we describe our conventions and notation.

## Single component systems

The Ornstein zernike equation 
$h(r) = c(r) + \rho \int d\textbf{r}' c(\textbf{r}')h(|\textbf{r} - \textbf{r}'|) $
links the pair correlation function $g(r) = h(r)+1$ to the number density $\rho$ and the direct correlation function $c(r)$. In order to solve it, a second relation between $h(r)$ and $c(r)$ must be specified. This is called the closure relation. 

In practise, the closure relation typically takes the form $c(r) = f(\gamma(r), r, u(r))$, where $u(r)$ is the interaciton potential and $\gamma(r) = h(r) - c(r)$ is the indirect correlation function. Sometimes, instead the closure is defined for the bridge function $b(r)$, and the closure relation is then given by $h(r) - 1 = \exp\left(-\beta u(r) + h(r) - c(r) + b(r) \right)$, where $\beta = 1/k_BT$. 

## Mixtures
Everything above generalizes to the mixture case:

The Ornstein zernike equation 

$h_{ij}(r) = c_{ij}(r) + \sum_l \rho_l \int d\textbf{r}' c_{il}(\textbf{r}')h_{lj}(|\textbf{r} - \textbf{r}'|)$

links the pair correlation function $g_{ij}(r) = h_{ij}(r)+1$ to the species specific number density $\rho_{i}$ and the direct correlation function $c_{ij}(r)$. In order to solve it, a second relation between $h_{ij}(r)$ and $c_{ij}(r)$ must be specified. This is called the closure relation. 

In practise, the closure relation typically takes the form $c_{ij}(r) = f(\gamma_{ij}(r), r, u_{ij}(r))$, where $u_{ij}(r)$ is the interaciton potential and $\gamma_{ij}(r) = h_{ij}(r) - c_{ij}(r)$ is the indirect correlation function. Sometimes, instead the closure is defined for the bridge function $b(r)$, and the closure relation is then given by $h_{ij}(r) - 1 = \exp\left(-\beta u_{ij}(r) + h_{ij}(r) - c_{ij}(r) + b_{ij}(r) \right)$, where $\beta = 1/k_BT$. 

## Fourier Transforms

The ability to numerically solve the Ornstein-Zernike equation relies heavily on doing repeated Fourier Transforms. In arbitrary dimensions, these Fourier transforms can be written as Hankel transforms in the case that the argument is a radially symmetric function. In particular, in $d$ dimensions, the radial Fourier transform and its inverse are given by

$\hat{F}(k) = (2\pi)^{d/2} k ^{1-d/2}\int_0^\infty dr r^{d/2}J_{d/2-1}(kr)F(r)$

$F(r) = (2\pi)^{-d/2} r ^{1-d/2}\int_0^\infty dk k^{d/2}J_{d/2-1}(kr)\hat{F}(k),$

in which $J_{d/2-1}(x)$ is the bessel function of order $d/2-1$. In the special cases of 1 and 3 dimensions, the transform simplifies into:

$\hat{F}(k) = 2\int_0^\infty dr \cos(kr)F(r)$

$F(r) = \frac{1}{\pi}\int_0^\infty dk \cos(kr)\hat{F}(k),$

in 1$d$, and 

$\hat{F}(k) = \frac{4\pi}{k}\int_0^\infty dr r \sin(kr)F(r)$

$F(r) = \frac{1}{2\pi^2r}\int_0^\infty dk k\sin(kr)\hat{F}(k),$

in 3$d$. This package uses discrete versions of all of the above. 

## Thermodynamic properties

Using the structure as determined by this package, several thermodynamic properies can be computed. In particular, this package contains methods to compute the (virial) pressure $p$, the isothermal compressibility $\chi$, and the excess internal energy per particle $E_x$.

For mixtures, they are computed respectively from the following definitions

$$p =  k_BT \rho_0\sum_i x_i - 1/6 \rho_0^2 \sum_{ij} x_i x_j \int d\textbf{r} r g_{ij}(r) u'_{ij}(r)$$

$$1/(\rho_0 k_BT \chi)  = 1 - ρ_0 \sum_{ij} x_i x_j \hat{c}_{ij}(k\to0)$$,

$$E_x =   1/2 \rho_0 \sum_{ij} x_i x_j  \int d\textbf{r}  g_{ij}(r) u_{ij}(r)$$

in which $\rho_0=N/V$. The functions to use are [`compute_virial_pressure`](@ref), [`compute_compressibility`](@ref), and, [`compute_excess_energy`](@ref).


# Thermodynamic Routes

The Ornstein–Zernike framework gives several thermodynamic quantities:

- **Virial route (pressure)**: from \(g(r)\) and \(u'(r)\).
- **Compressibility route**:
  - `compute_compressibility` returns **isothermal compressibility** \(\chi\).
  - Pressure via compressibility route requires **integration over density** using
    \[
      \frac{\partial (\beta p)}{\partial \rho} \;=\; \frac{1}{S(0)}, 
      \qquad S(0) = \rho\,k_BT\,\kappa_T,
    \]
    so
    \[
      \beta p(\rho) \;=\; \int_0^\rho \frac{d\rho'}{S(0;\rho')}.
    \]
- **Energy route**: excess internal energy from \(g(r)\).

Because closures are approximate, routes generally **disagree**; the gap is a measure of closure error.

---

## Example: Virial Pressure vs. Compressibility-Route Pressure

Below we:
1) compute **virial pressure** directly,  
2) compute \(\chi\) at a sequence of densities,  
3) integrate \(1/S(0)\) to obtain the **compressibility-route pressure**.

```@example 9
using OrnsteinZernike

# --- System definition at target temperature and potential ---
kBT     = 1.0
σ, ϵ    = 1.0, 1.0
closure = Verlet()   # any closure works; V shown here
potential = InversePowerLaw(ϵ, σ, 8)

# Choose a target density and also a ramp to integrate from 0 → ρ_target
ρ_target   = 0.8
delta_ρ = 0.01
ρ_grid     = collect((delta_ρ/2):delta_ρ:(ρ_target-delta_ρ/2))  # avoid 0 to keep numerics stable
system_at = ρ -> SimpleFluid(3, ρ, kBT, potential)

# (1) Virial-route pressure at ρ_target
sol_target = solve(system_at(ρ_target), closure, NgIteration(M=5000, dr=0.01))
p_virial   = compute_virial_pressure(sol_target, system_at(ρ_target))   # this is "virial pressure"

# (2) Compressibility κ_T along a density grid
chi_values = similar(ρ_grid)
sols = solve(system_at(ρ_grid[end]), closure, DensityRamp(NgIteration(M=5000, dr=0.01), ρ_grid))

for (i, ρ) in enumerate(ρ_grid)
    chi_values[i] = compute_compressibility(sols[i], system_at(ρ))  # Isothermal compressibility κ_T
end

# (3) Integrate to get compressibility-route pressure:
#     d(βp)/dρ = 1 / S(0) with S(0) = ρ * kBT * κ_T
S0 = ρ_grid .* kBT .* chi_values

# Midpoint integration of 1/S0 over ρ to get β p(ρ)
βp  = delta_ρ * sum(1 ./ S0) 
p_comp = kBT * βp    # convert βp → p

println("Virial-route pressure at          ρ = $(ρ_target): p_virial = $(p_virial)")
println("Compressibility-route pressure at ρ = $(ρ_target): p_comp   = $(p_comp)")
```

**Notes**
- `compute_compressibility(sol, system)` returns \(\chi\) (units of inverse pressure).  
- We used a **density ramp** inside the loop to help convergence and reuse previous solutions.  

---

## Why Routes Differ

- **Virial route** is sensitive to short-range structure (contact region).  
- **Compressibility route** probes long-wavelength fluctuations via \(S(0)\).  
- Discrepancies reflect **closure approximation error**.

---

## Making Routes Agree

Closures like **Rogers–Young** include a tunable parameter (e.g., \(\alpha\)) that can be chosen so that
\( p_{\mathrm{vir}} \approx p_{\mathrm{comp}} \). See the *Thermodynamic Consistency* tutorial for a bisection example.
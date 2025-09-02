# How `OrnsteinZernike.jl` Solves Charged Systems

Charged fluids have long-ranged Coulomb interactions that require special handling in the Ornstein–Zernike (OZ) solver. This page explains the decomposition used in the code, how the *short-ranged* (SR) and *long-ranged* (LR) pieces are propagated through the OZ relations, and how the solver iterates to convergence.

---

## Field Splitting: Definitions and Notation

Potentials and correlation functions are consistently decomposed into **short-ranged** and **long-ranged** parts:

- **Potential**
  \[
  \beta u_{\text{coul}} = \underbrace{\beta u_{\text{SR,coul}}}_{\text{Coulomb SR split}}
  + \underbrace{\beta u_{\text{LR,coul}}}_{\text{Coulomb LR split}}.
  \]

- **Direct correlation**
  \[
  c = c_{\text{short\_range}} + \Phi,
  \qquad \Phi \equiv -\,\beta u_{\text{LR,coul}}.
  \]

- **Indirect correlation**
  \[
  \gamma = \gamma_{\text{short\_range}} + \gamma_{\text{long\_range}}.
  \]

- **Total correlation**
  \[
  h = h_{\text{short\_range}} + q,
  \]
  where the LR part \(q\) is defined by the **LR OZ relation** in Fourier space:
  \[
  \widehat{q}(k) = \widehat{\Phi} + \widehat{\Phi} \rho \widehat{q}(k).
  \]

---

## Coulomb Splitting in Practice

The Coulomb potential \(z_i z_j \,\ell_B / r\) is split into SR and LR pieces by a user-selectable strategy. Other strategies can be easily implemented by a user.

```@docs
NoCoulombSplitting
EwaldSplitting
```

The SR/LR split is produced by `split_coulomb_potential(r, system, coulombsplitting)`.

Given \(\Phi\) and its transform, the LR OZ part is solved analytically in \(k\)-space. This gives a fixed LR reference \(q\) used throughout the SR iteration.

---


## Minimal Usage Example: 1-3 electrolyte

```@example charges
using OrnsteinZernike, StaticArrays, Plots

# Define a (neutral) base mixture first
dims = 3
ρ = [0.6, 0.2]              # number densities (reduced)
Z = [1, -3]                   # charges
ℓB = 7.0                      # Bjerrum length (reduced)

pot = HardSpheres([1.0, 0.5])
sys = SimpleMixture(dims, ρ, 1, pot)

# Turn it into a charged mixture and choose a Coulomb splitting
charged = SimpleChargedMixture(sys, Z, ℓB)

closure = HypernettedChain()
method  = FourierIteration(M=4096, dr=0.01, tolerance=1e-8, mixing_parameter=0.5)

# Solve with an explicit splitting choice (e.g., NoCoulombSplitting or EwaldSplitting(α))
sol = solve(charged, closure, method; coulombsplitting=NoCoulombSplitting())
# e.g.: sol = solve(charged, closure, method; coulombsplitting=EwaldSplitting(α=3.0))
# plot the unlike charges:
plot(sol.r, sol.gr[:, 1, 2], xlims=(0,3), label="g_{12}(r)", xlabel="r", ylabel="g(r)")
plot!(sol.r, sol.gr[:, 1, 1], label="g_{11}(r)")
plot!(sol.r, sol.gr[:, 2, 2], label="g_{22}(r)")
```

We see that the unlike charges have a strong peak in g(r).

**Tip:** If convergence is delicate, try:
- increasing `M` (grid size) and/or decreasing `dr`,
- reducing the `mixing_parameter` or changing the number of stages if using `NgIteration`.


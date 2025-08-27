# Mixtures in OrnsteinZernike.jl

For a mixture with \(N_s\) species, the Ornstein–Zernike equation becomes matrix-valued:  

\[
h_{ij}(r) = c_{ij}(r) + \sum_{k=1}^{N_s} \rho_k \, (c_{ik} * h_{kj})(r)
\]

Outputs such as \(g_{ij}(r)\), \(S_{ij}(k)\), \(c_{ij}(r)\) are stored as 3D arrays:  
- `gr[r, i, j]` = pair distribution between species *i* and *j*.  
- `Sk[k, i, j]` = partial structure factor between species *i* and *j*.

---

## Example: Binary Lennard–Jones Mixture

```@example 10
using OrnsteinZernike, Plots, StaticArrays

ρ = [0.3, 0.2]  # densities of species A and B

function lj_matrix(r, params)
    ϵ_AA, σ_AA, ϵ_BB, σ_BB, ϵ_AB, σ_AB = params
    return SMatrix{2,2, Float64, 4}(
        4*ϵ_AA*((σ_AA/r)^12 - (σ_AA/r)^6),  4*ϵ_AB*((σ_AB/r)^12 - (σ_AB/r)^6),
        4*ϵ_AB*((σ_AB/r)^12 - (σ_AB/r)^6),  4*ϵ_BB*((σ_BB/r)^12 - (σ_BB/r)^6)
    )
end

pot = CustomPotential(lj_matrix, (1.0,1.0, 1.2,0.8, 1.1,0.9))
system = SimpleMixture(3, ρ, 1.0, pot)

sol = solve(system, HypernettedChain())

plot(sol.r, sol.gr[:,1,1], label="g_AA(r)")
plot!(sol.r, sol.gr[:,1,2], label="g_AB(r)")
plot!(sol.r, sol.gr[:,2,2], label="g_BB(r)")
```
---

## Interpretation

- \(g_{AA}(r)\): correlations within species A  
- \(g_{BB}(r)\): correlations within species B  
- \(g_{AB}(r)\): cross-correlations  

These can be combined into total structure factors, or into Bhatia–Thornton components.

# Choosing Grid Parameters (M and dr)

The Ornstein–Zernike solvers use **radial Fourier/Hankel transforms** to go between real and Fourier space.  
The grid is controlled by two parameters:

- `M` : number of grid points  
- `dr`: approximate radial step size (may not be truly uniform, since Hankel transforms use non-uniform grids in arbitrary dimensions)  

The **maximum radius** is

\[
r_{\max} = M \cdot dr
\]

which sets the size of the box in real space.

---

## What Must Be Resolved?

- In **real space**, it is the **direct correlation function \(c(r)\)** that is Fourier transformed.  
  - Therefore, `dr` and `rmax` must be chosen such that oscillations and decay of \(c(r)\) are well resolved.  
  - Poor resolution of \(c(r)\) will produce inaccurate structure factors \(S(k)\).
  - In a three-dimensional system, if the potential has discontinuities, they must lie exactly between two gridpoints to ensure second order convergence.

- In **Fourier space**, the discretization also needs to be accurate.  
  - If the \(k\)-grid is too coarse, long-wavelength properties (e.g. \(S(k\to 0)\)) will be unreliable. Ensure that the direct correlation function \(\gamma(k)\) is well discretized.
  - If \(r_{\max}\) is too small, aliasing and finite-box effects may distort the Fourier transform.

---

## Practical Guidelines

- **Choose dr small enough** to resolve short-range features of the potential and resulting \(c(r)\).  
  - For hard spheres or Lennard–Jones: `dr ≲ 0.05σ`.  
- **Choose rmax large enough** so that \(c(r)\) has decayed essentially to zero before the cutoff.  
  - A good starting point: `rmax ≈ 10σ`.  
- **Increase M** together with decreasing dr if necessary — this improves both real-space and Fourier-space resolution.  
- **Check convergence**: repeat a calculation with finer grid and verify \(g(r)\), \(S(k)\), and thermodynamic routes are stable.

---

## Example: Grid Convergence with dr

```@example 8
using OrnsteinZernike, Plots

system = SimpleLiquid(3, 0.8, 1.0, LennardJones(1.0, 1.0))
closure = HypernettedChain()

Rmax = 10.0
p1 = plot(xlims=(0,2), xlabel="r", ylabel="c(r)")
p2 = plot(xlims=(0,20), xlabel="k", ylabel="c(k)")
for M in [1000, 100, 50]
    dr = Rmax/M
    sol  = @time solve(system, closure, NgIteration(M=M, dr=dr))
    kmax = round(maximum(sol.k), digits=2)
    plot!(p1, sol.r, sol.cr, label="M=$M, kmax=$(kmax)", ls=:dash, markers=:o, legend=false)
    plot!(p2, sol.k, sol.ck, label="M=$M, kmax=$(kmax)", ls=:dash, markers=:o, legend=:bottomright)
end
display(plot(p1,p2))
```

Here the coarse grid under-resolves oscillations in \(c(r)\) and doesn't capture it's decay in Fourier space, which in turn leads to distorted structure factors \(S(k)\).

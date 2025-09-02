# Tutorial: Reproducing Electrolyte HNC Results from Literature

In this tutorial, we demonstrate how to use `OrnsteinZernike.jl` to reproduce classical results for symmetric electrolytes obtained with the **Hypernetted Chain (HNC)** closure. Specifically, this page will reproduce some of the results presented in.

> D.‐M. Duh; A. D. J. Haymet, *Integral equation theory for charged liquids: Model 2–2 electrolytes and the bridge function*, J. Chem. Phys. **97**, 7716–7729 (1992).

This paper studies **2:2 electrolytes** (e.g., divalent cations and anions) using integral equation theory. Our goal here is to set up the same model, run the HNC closure, and reproduce the radial distribution functions (RDFs) and excess energies across different concentrations,.  

---

## Step 1. Setting up the environment

```@example elec
using OrnsteinZernike, StaticArrays
using Plots
```

---

## Step 2. Define physical constants

We work with a set of standard physical constants: particle size `σ`, dielectric constant `ϵr`, elementary charge `e`, vacuum permittivity `ϵ0`, and Boltzmann constant `kB`. The parameter `B` sets the strength of the short-ranged repulsive potential. The values are taken from the paper above.

```@example elec
σ = 2.8428 * 10^-10      # particle size (m)
ϵr = 78.358              # relative permittivity of water at 25 °C
e = 1.602176634e-19      # elementary charge (C)
ϵ0 = 8.8541878128e-12    # vacuum permittivity (F/m)
ϵ = ϵ0 * ϵr
kB = 1.380649e-23        # Boltzmann constant (J/K)
B = 5377.75 * 4 * 10^-10 # repulsion parameter (m K)
```
The potential used is given by

$$u_{ij}(r) = \frac{k_B B}{\sigma} * \left(\frac{\sigma}{r}\right)^9 +  \frac{Z_i Z_j e^2}{4\pi\epsilon r}$$

For convenience, we can also compute the approximate position of the potential minimum (not strictly needed, but useful for diagnostics):

```@example elec
Rm = (3^(1/4) * B^(1/8) * kB^(1/8)*(4*π*ϵ)^(1/8)*σ)/(e^(1/4)*2^(1/4))
```

---

## Step 3. Prepare plotting

We will plot two RDFs: the like–like distribution function `g₁₁(r)` and the unlike distribution function `g₁₂(r)`.

```@example elec
p1 = plot(xlabel="r (Å)", ylab="g₁₁(r)")
p2 = plot(xlabel="r (Å)", ylab="g₁₂(r)")
nothing # hide
```

---

## Step 4. Ion charges and thermodynamics

We consider a symmetric 2:2 electrolyte, so the charges are `+2` and `–2`. The temperature is set to room temperature (298 K). We also prepare a list of mixing parameters to ensure numerical convergence of the iterative solver.

```@example elec
Z = [2, -2] 
T0 = 273.15
kBT = kB * (T0 + 25)   # thermal energy at 298 K

mixing_parameters = [0.05, 0.005, 0.05, 0.05, 0.05, 0.05]

# Bjerrum length in meters
bjerrum_length = e^2/(4*π*ϵ*kBT) 
```

---

## Step 5. Loop over concentrations

We now study several salt concentrations (in mol/L). For each concentration:

1. Convert to number density (per m³).
2. Build the mixture model.
3. Solve the Ornstein–Zernike equation with the HNC closure.
4. Extract the RDFs and excess energy.
5. Add the results to the plots.

```@example elec
for (i,c) in enumerate([0.001, 0.005, 0.02, 0.0625, 0.2, 0.5625]) # mol/L
    # Convert molar concentration to number density
    ρ0 = c * 6.02214076e23 * 1e3  
    ρ = [ρ0, ρ0] # equal density of cations and anions

    # Compute Debye screening parameter
    κD = sqrt(1/(ϵr * ϵ0 * kBT)* e^2 * sum(ρ .* Z.^2))
    @show c, Rm * κD

    # Numerical grid settings
    # From here we work in units of kBT and σ
    M = 4096
    dims = 3
    dr = 0.105 * 10^-10 / σ  
    Sones = ones(SMatrix{2,2, Float64, 4})

    # Define short-range repulsive potential
    pot = PowerLaw(kB*B/σ*Sones/kBT, Sones, 9)

    # Build system and add Coulomb interactions
    system = SimpleMixture(dims, ρ*σ^3, 1, pot)
    system = SimpleChargedMixture(system, Z, bjerrum_length / σ)

    # Choose closure and numerical method
    closure = HypernettedChain()
    method = FourierIteration(M = M,
                              dr = dr,
                              tolerance = 1e-8,
                              mixing_parameter = mixing_parameters[i],
                              verbose = false,
                              max_iterations = 10^6)

    # Solve OZ equation with HNC closure
    sol = solve(system, closure, method,
                coulombsplitting = OrnsteinZernike.NoCoulombSplitting())

    # Extract RDFs
    r = sol.r 
    g = sol.gr

    # Compute excess energy (per particle, in reduced units)
    Ex = OrnsteinZernike.compute_excess_energy(sol, system)
    @show -Ex

    # Plot g₁₁(r) and g₁₂(r)
    plot!(p1, r*1e10*σ, g[:, 2, 2], label="c=$(c)M")
    plot!(p2, r*1e10*σ, g[:, 1, 2], label="c=$(c)M")
end
```

---

## Step 6. Display results

Finally, we set axis limits and show the combined plots.

```@example elec
plot!(p2, xlims=(0.0, 16.0))
plot!(p1, xlims=(0.0, 16.0))

display(plot(p1, p2, layout=(2,1)))
```

---

## Results

Running this script produces the ion–ion radial distribution functions at different concentrations. These curves reproduce the exact behavior reported in Duh & Haymet (1992). The excess energies match exactly the data from table III (except for one data point for some reason).  


# Systems

A `System` describes *what* you want to solve: dimensionality, density, temperature, and the interaction potential.  
The package exports two **neutral** system types and two **charged** wrappers:

- `SimpleFluid` — single component
- `SimpleMixture` — multi-component (mixtures)
- `SimpleChargedFluid` — **one-component plasma (OCP)**: a single mobile charged species in a *uniform neutralizing background*
- `SimpleChargedMixture` — electrolyte mixture: multiple mobile charged species (no background; must be electroneutral)

!!! tip "Choosing a system"
    | Type                     | Components | ρ (density) input           | `evaluate_potential` return |
    |--------------------------|------------|-----------------------------|-----------------------------|
    | `SimpleFluid`            | 1          | `Number`                    | `Number`                    |
    | `SimpleMixture`          | ≥ 2        | `AbstractVector` (length = Ns) | `Ns×Ns` matrix (e.g. `SMatrix`) |
    | `SimpleChargedFluid`     | 1 + background | (in `base`)             | (from `base`)               |
    | `SimpleChargedMixture`   | ≥ 2        | (in `base`)                 | (from `base`)               |

Internally, for mixtures the density vector is stored as a diagonal matrix with `StaticArrays.SVector` storage for performance.

---

## `SimpleFluid`

`SimpleFluid` assumes rotationally symmetric pair interactions (depend only on the center-to-center distance) and no external field.

```@docs
SimpleFluid
```

### Example: 2D Lennard–Jones (single component)

```@example 1
using OrnsteinZernike
potential = LennardJones(1.0, 1.0)  # ϵ=1, σ=1
kBT = 1.0
ρ   = 0.5
dims = 2
system = SimpleFluid(dims, ρ, kBT, potential)
```

---

## `SimpleMixture`

`SimpleMixture` is for multi-component systems. Pass a vector of species densities; the potential must evaluate to an `Ns×Ns` matrix giving all pair interactions.

```@docs
SimpleMixture
```

!!! note
    The density vector `ρ::AbstractVector` is converted to `Diagonal(SVector{Ns}(ρ))` internally.

### Example: 3D 10-component hard-sphere mixture

```@example 1
using OrnsteinZernike # hide
D        = collect(0.1:0.1:1.0)        # species diameters
potential = HardSpheres(D)
kBT = 1.0
ρ   = fill(0.1, 10)                     # number densities per species
dims = 3
system = SimpleMixture(dims, ρ, kBT, potential)
```

---

## Charged systems

### `SimpleChargedFluid` — One-Component Plasma (OCP)

`SimpleChargedFluid` models a **single mobile charged species** embedded in a **uniform neutralizing background** (“jellium”).  
This means you do **not** need a counter-ion species: electroneutrality is enforced by the background.

- Background charge density is implicitly `ρ_bg = − ρ · z` (not mobile, not part of the potential; it only neutralizes and fixes the small-k behavior).
- Debye screening uses *mobile* charges only: `κ_D^2 = 4π ℓ_B ρ z^2` (in 3D reduced units).
- The solver handles the `k→0` Coulomb singularity via an Ewald-like split; background terms cancel the divergent constants so structure is well-defined.
- Use `SimpleChargedMixture` instead if you want explicit co-/counter-ions (no background, mixture must be electroneutral).

```@docs
SimpleChargedFluid
```

#### Example: 3D classical OCP (hard-sphere core + Coulomb + neutralizing background)

```@example 1
using OrnsteinZernike # hide
dims = 3
ρ    = 0.5
kBT  = 1.0
core = HardSpheres(1.0)                 # short-range core (optional)
base = SimpleFluid(dims, ρ, kBT, core)

# OCP: single species of charge z in a uniform neutralizing background
z    = 1.0
lB = 7.0                                # Bjerrum length
sys  = SimpleChargedFluid(base, z, lB)  # κ chosen automatically (Debye)

# Now use closures/solvers as usual:
# sol = solve(sys, HypernettedChain(); method=NgIteration(M=2000, dr=0.01))
```

### `SimpleChargedMixture` — Electrolyte mixtures

`SimpleChargedMixture` wraps a neutral `SimpleMixture` and adds per-species charges `z`.  
No background is added; the system must satisfy electroneutrality `∑ᵢ ρᵢ zᵢ ≈ 0`.

```@docs
SimpleChargedMixture
```

#### Example: 3D 1:1 restricted primitive model (RPM)

```@example 1
using OrnsteinZernike # hide
dims = 3
ρ    = [0.3, 0.3]                      # mobile species densities
kBT  = 1.0
hs   = HardSpheres([1.0, 1.0])         # per-species diameters (short-range)
base = SimpleMixture(dims, ρ, kBT, hs)

z    = [ 1.0, -1.0 ]                   # charges (must be electroneutral with ρ)
lB = 7.0                               # Bjerrum length
sys  = SimpleChargedMixture(base, z, lB) 
```

---

## Backward compatibility

The old `SimpleLiquid` constructors are still available **as deprecated aliases**:

- `SimpleLiquid(dims, ρ::Number, kBT, potential) → SimpleFluid`
- `SimpleLiquid(dims, ρ::Vector, kBT, potential) → SimpleMixture`

You’ll see a deprecation warning; please migrate to the new names.

---

## Common pitfalls

- **Density/potential mismatch**:  
  If `ρ` is a vector of length `Ns`, then `evaluate_potential(potential, r)` must return an `Ns×Ns` matrix (e.g., `SMatrix{Ns,Ns}`).
- **Units**:  
  Ensure the potential and `kBT` are in consistent reduced units.
- **For OCP (`SimpleChargedFluid`)**:  
  Remember the neutralizing background is implicit and non-mobile; use *mixtures* if you need explicit counter-ions.
- **For `SimpleChargedMixture`**:  
  The mixture must be electroneutral: `∑ᵢ ρᵢ zᵢ ≈ 0`.


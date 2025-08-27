# Systems

A `System` describes the physical system you want to solve: dimensionality, density, temperature, and the interaction potential.  
The package exports two **neutral** system types and two **charged** wrappers:

- `SimpleFluid` — single component
- `SimpleMixture` — multi-component (mixtures)
- `SimpleChargedFluid` — single component + electrostatics
- `SimpleChargedMixture` — mixture + electrostatics

!!! tip "Choosing a system"
    | Type                   | Components | ρ (density) input           | `evaluate_potential` return |
    |------------------------|------------|-----------------------------|-----------------------------|
    | `SimpleFluid`          | 1          | `Number`                    | `Number`                    |
    | `SimpleMixture`        | ≥ 2        | `AbstractVector` (length = Ns) | `Ns×Ns` matrix (e.g. `SMatrix`) |
    | `SimpleChargedFluid`   | 1          | (in `base`)                 | (from `base`)               |
    | `SimpleChargedMixture` | ≥ 2        | (in `base`)                 | (from `base`)               |

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

## Charged systems (electrolytes)

Electrostatics are modeled at the **system level** via wrappers around neutral systems:

- `SimpleChargedFluid(base; z, εr=78.4, κ=:auto)`
- `SimpleChargedMixture(base; z, εr=78.4, κ=:auto)`

Here `z` are charges (units of *e*), `εr` is the relative dielectric, and `κ` is the Gaussian/Ewald split (defaults to Debye).

```@docs
SimpleChargedFluid
SimpleChargedMixture
```

!!! warning "Electroneutrality"
    `SimpleChargedMixture` checks ∑ᵢ ρᵢ zᵢ ≈ 0 and will throw if violated.

### Example: 3D 1:1 restricted primitive model (RPM)

```@example 1
using OrnsteinZernike # hide
# Base mixture (short-range): equal-diameter hard spheres
dims = 3
ρ    = [0.3, 0.3]                      # species densities
kBT  = 1.0
hs   = HardSpheres([1.0, 1.0])         # per-species diameters
base = SimpleMixture(dims, ρ, kBT, hs)

# Add electrostatics on the system level
z    = [1.0, -1.0]
sys  = SimpleChargedMixture(base; z=z, εr=78.4)  # κ auto-set to Debye

# Now you can call `solve(sys, closure)` just like for neutral systems.
nothing # hide
```

### Example: 3D charged single component

```@example 1
using OrnsteinZernike # hide
dims = 3
ρ    = 0.5
kBT  = 1.0
pot  = HardSpheres(1.0)
base = SimpleFluid(dims, ρ, kBT, pot)
sys  = SimpleChargedFluid(base; z=1.0, εr=78.4)
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
- **Charged mixtures must be neutral**:  
  `SimpleChargedMixture` enforces electroneutrality via ∑ᵢ ρᵢ zᵢ ≈ 0.


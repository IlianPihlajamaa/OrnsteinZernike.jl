abstract type CoulombSplitting end

"""
    EwaldSplitting(α) <: CoulombSplitting

Splits the Coulomb potential into short-range and long-range parts using the Ewald splitting parameter `α`.
Fields:
- α::Float64 : Ewald splitting  (inverse length scale)

The short-range part is given by:
    u_short_range(r) = (z_i * z_j * ℓB / r) * erfc(α * r)
The long-range part is given by:
    u_long_range(r) = (z_i * z_j * ℓB / r) * erf(α * r)
"""
struct EwaldSplitting <: CoulombSplitting
    α::Float64
end

"""
    NoSplitting <: CoulombSplitting

No splitting of the Coulomb potential; the entire potential is treated as long-range.
"""
struct NoCoulombSplitting <: CoulombSplitting end



# ---- Dimension constant for the Coulomb kernel's FT:  v̂(k) = (C_d / k^2) / (ε_r β) ----
# (Matches the surface area of the unit sphere: C_3=4π, C_2=2π, C_1=2.)
C_d(d::Integer) = d == 3 ? 4π : d == 2 ? 2π : d == 1 ? 2 : error("Unsupported dimension d=$d")

# ---- Debye screening parameter κ_D  (scalar, single-species) ----
# κ_D^2 = C_d * ℓB * ρ * z^2  ⇒  κ_D = sqrt(C_d * ℓB * ρ * z^2)
debye_kappa(ρ::Number, z::Number, ℓB; dims::Integer=3) =
    sqrt(C_d(dims) * ℓB * ρ * z^2)

# ---- Debye κ_D for mixtures (vector densities) ----
# Accept either a Diagonal of densities or a plain vector; z is the vector of valences.
debye_kappa(ρdiag::Diagonal, z::AbstractVector, ℓB; dims::Integer=3) =
    sqrt(C_d(dims) * ℓB * sum(diag(ρdiag) .* (z .^ 2)))

debye_kappa(ρ::AbstractVector, z::AbstractVector, ℓB; dims::Integer=3) =
    sqrt(C_d(dims) * ℓB * sum(ρ .* (z .^ 2)))

function evaluate_coulomb_potential(r::Real, system::SimpleChargedSystem)
    @assert dims_of(system) == 3 "Coulomb potential only implemented in 3D."
    z = system.z
    βu = z * z' * system.ℓB / r
    return βu
end

function evaluate_coulomb_potential(r::AbstractVector, system::SimpleChargedSystem)
    return [evaluate_coulomb_potential(ri, system) for ri in r]
end

function split_coulomb_potential(r::AbstractVector, system, CS::CoulombSplitting)
    usr1, ulr1 = split_coulomb_potential(first(r), system, CS)
    usr = typeof(usr1)[usr1]
    ulr = typeof(ulr1)[ulr1]
    for i in (firstindex(r)+1):lastindex(r)
        usri, ulri = split_coulomb_potential(r[i], system, CS)
        push!(usr, usri)
        push!(ulr, ulri)
    end

    return usr, ulr 
end

function split_coulomb_potential(r::Real, system::SimpleChargedSystem, ::NoCoulombSplitting)
    βu = evaluate_coulomb_potential(r, system)
    return zero(βu), βu
end

function split_coulomb_potential(r::Real, system::SimpleChargedSystem, CS::EwaldSplitting)
    @assert dims_of(system) == 3 "EwaldSplitting of Coulomb potential only implemented in 3D."
    α = CS.α
    z = system.z
    ℓB = system.ℓB
    erfpart = erf(α * r)
    βu_long_range = z * z' * ℓB * erfpart / r
    βu_short_range = z * z' * ℓB * (1 - erfpart) / r
    return βu_short_range, βu_long_range
end
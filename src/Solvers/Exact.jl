function find_Ck_3dPYSC_exact(η, k)
    A = -(1 - η)^-4 * (1 + 2η)^2
    B = (1 - η)^-4 * 6η * (1 + η / 2)^2
    DD = -(1 - η)^-4 * 1 / 2 * η * (1 + 2η)^2
    Ck = @. 4π / k^6 *
            (
        24 * DD - 2 * B * k^2 - (24 * DD - 2 * (B + 6 * DD) * k^2 + (A + B + DD) * k^4) * cos(k)
        +
        k * (-24 * DD + (A + 2 * B + 4 * DD) * k^2) * sin(k)
    )
    return Ck
end

function find_Cr_3dPYSC_exact(η, r)
    Cr = @. -(1 - η)^-4 * ((1 + 2η)^2 - 6η * (1 + η / 2)^2 * r + 0.5 * η * (1 + 2η)^2 * r^3)
    Cr[r.>1.0] .= 0.0
    return Cr
end

function solve(system::SimpleLiquid{3, 1, T1, T2, HardSpheres{T3}}, ::PercusYevick, method::Exact) where {T1,T2,T3}
    # D = 1 by definition
    @assert system.potential.D == 1.0 "This method assumes that the hard sphere diameter D = 1.0"
    r, k = construct_r_and_k_grid(system, method)

    ρ = system.ρ
    η = ρ * π / 6

    # baxter
    Cr = find_Cr_3dPYSC_exact(η, r)

    # analytical tranform of the above:
    Ck = find_Ck_3dPYSC_exact(η, k)

    Sk = @. 1.0 + ρ * Ck / (1 - ρ * Ck)

    # γ is a smooth function, therefore it has a nicer transform
    γmulk = @. ρ * (Ck * k)^2 / (k - ρ * Ck * k)
    γmulr = inverse_radial_fourier_transform_3d(γmulk, r, k)

    gr = @. Cr + γmulr / r + 1
    gr[r.<1.0] .= 0.0
    return OZSolution(r, k, gr, Sk, Ck, Cr)
end


function solve_3dPYMC_exact(ρ, diameters, kᵢ::Number)
    Ns = length(ρ)
    p = length(diameters)
    @assert Ns == p
    d = (diameters .+ diameters') / 2
    s = (diameters .- diameters') / 2
    ξ = [π / 6 * sum(ρ[j] * d[j, j]^ν for j = 1:p) for ν in 1:3]
    a = [(1 - ξ[3])^(-2) * (1 - ξ[3] + 3 * ξ[2] * d[i, i]) for i = 1:p]
    b = [-3 / 2 * d[i, i]^2 * (1 - ξ[3])^(-2) * ξ[2] for i = 1:p]
    Q̃ = zeros(ComplexF64, Ns, Ns)
    
    #analytical solution of integrals
    for μ = 1:Ns, ν = 1:Ns
        I0 = -1im / kᵢ * (cis(kᵢ * d[μ, ν]) - cis(kᵢ * s[μ, ν]))
        I1 = -1im / kᵢ * (d[μ, ν] * cis(kᵢ * d[μ, ν]) - s[μ, ν] * cis(kᵢ * s[μ, ν])) + 1im / kᵢ * I0
        I2 = -1im / kᵢ * (d[μ, ν]^2 * cis(kᵢ * d[μ, ν]) - s[μ, ν]^2 * cis(kᵢ * s[μ, ν])) + 2im / kᵢ * I1
        I3 = 0.5 * a[μ] * I2 + b[μ] * I1 - (0.5 * a[μ] * d[μ, ν]^2 + b[μ] * d[μ, ν]) * I0
        Q̃[μ, ν] = I[μ, ν] - 2π * sqrt(ρ[μ] * ρ[ν]) * I3
    end

    C̃kᵢ = I - real.(Q̃' * Q̃)
    H̃kᵢ = (I - C̃kᵢ) \ C̃kᵢ

    Hkᵢ = H̃kᵢ ./ sqrt.(ρ .* ρ')
    Ckᵢ = C̃kᵢ ./ sqrt.(ρ .* ρ')
    return Ckᵢ, Hkᵢ
end


# """
# ref: Baxter, R.J. Ornstein–Zernike Relation and Percus–Yevick Approximation for Fluid Mixtures, J. Chem. Phys. 52, 4559 (1970)
# """
function solve(system::SimpleLiquid{3, species, T1, T2, HardSpheres{T3}}, ::PercusYevick, method::Exact) where {species, T1, T2, T3<:AbstractMatrix}
    r, k = construct_r_and_k_grid(system, method)
    ρ = (system.ρ).diag
    T = eltype(ρ)
    Ns = species
    x = ρ / sum(ρ)
    Ck = zeros(length(k), length(ρ), length(ρ))
    Hk = zeros(length(k), length(ρ), length(ρ))
    Sk = zeros(length(k), length(ρ), length(ρ))
    gr = zeros(length(k), length(ρ), length(ρ))
    Cr = zeros(length(k), length(ρ), length(ρ))
    diameters = diag(system.potential.D)

    @assert all((diameters .+ diameters')/2 .≈ system.potential.D) "This exact method does not work for non-additive diameters"

    p = length(diameters)
    @assert Ns == p
    d = (diameters .+ diameters') / 2

    for (iₖ, kᵢ) in enumerate(k)
        Ckᵢ, Hkᵢ = solve_3dPYMC_exact(ρ, diameters, kᵢ)
        for μ = 1:Ns, ν = 1:Ns
            Ck[iₖ, μ, ν] = Ckᵢ[μ, ν]
            Hk[iₖ, μ, ν] = Hkᵢ[μ, ν]
            Sk[iₖ, μ, ν] = x[μ] * I[μ, ν] + x[μ] * x[μ] * sum(ρ) * Hkᵢ[μ, ν]
        end
    end

    for μ = 1:Ns, ν = 1:Ns
        expnegβu = (r .> d[μ, ν])
        γμν = inverse_radial_fourier_transform_3d((Hk[:, μ, ν] - Ck[:, μ, ν]) .* k, r, k) ./ r
        @. gr[:, μ, ν] = (one(T) + γμν) * expnegβu
        @. Cr[:, μ, ν] = (expnegβu - one(T)) * (1 + γμν)
    end
    return OZSolution(r, k, gr, Sk, Ck, Cr)
end
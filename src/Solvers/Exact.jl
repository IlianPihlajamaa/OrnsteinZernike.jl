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

function find_Ck_1dPYSC_exact(η, k)
    Q0 = -1/(1-η)
    c0 = -Q0^2
    c1 = η * Q0^2
    Ck = @. (2*(-c1 + c1*cos(k) + (c0 + c1) * k* sin(k)))/k^2
    return Ck
end

function find_Cr_1dPYSC_exact(η, r)
    Q0 = -1/(1-η)
    c0 = -Q0^2
    c1 = η * Q0^2
    Cr = @. c0 + c1*r
    Cr[r.>1.0] .= 0.0
    return Cr
end

function find_Ck_5dPYSC_exact(η, k)
    T = 1+18η+6η^2
    Q0 = 1/(120η*(1-η)^3)*(1-33η-87η^2-6η^3-T^(3/2))
    Q1 = -1/(12*(1-η)^3)*((3+2η)*T^(1/2) + 3+19η+3η^2)
    Q2 =-T^(1/2)/(24(1-η)^3)*(2+3η+T^(1/2))

    Q̃0 = -8Q2
    c0 = -(Q̃0)^2
    c1 = 120η*Q0^2
    c3 = 20η*(8Q0*Q2-3Q1^2)
    c5 = -3/8*η*c0
    Ck = @. -(1/(k^10))*
        8*π^2 * (8 * (720c5 - 18c3 * k^2 + c1 * k^4) + (-5760c5 +
        144*(c3 + 20c5) * k^2 - 8*(c1 + 9c3 + 30c5) * k^4
        + (3c0 + 4c1 + 6c3 + 8c5) * k^6)*cos(k) +
        k * (-5760c5 + 48*(3c3 + 20c5) * k^2 -
        (3c0 + 8*(c1 + 3c3 + 6c5)) * k^4 +
        (c0 + c1 + c3 + c5) * k^6) * sin(k))

    return Ck
end

function find_Cr_5dPYSC_exact(η, r)
    T = 1+18η+6η^2
    Q0 = 1/(120η*(1-η)^3)*(1-33η-87η^2-6η^3-T^(3/2))
    Q1 = -1/(12*(1-η)^3)*((3+2η)*T^(1/2) + 3+19η+3η^2)
    Q2 =-T^(1/2)/(24(1-η)^3)*(2+3η+T^(1/2))

    Q̃0 = -8Q2
    c0 = -(Q̃0)^2
    c1 = 120η*Q0^2
    c3 = 20η*(8Q0*Q2-3Q1^2)
    c5 = -3/8*η*c0
    Cr = @. c0 + c1*r + c3*r^3 + c5*r^5
    Cr[r.>1.0] .= 0.0
    return Cr
end



function solve(system::SimpleFluid{3, T1, T2, HardSpheres{T3}}, ::PercusYevick, method::Exact) where {T1,T2,T3}
    # D = 1 by definition
    @assert system.potential.D == 1.0 "This method assumes that the hard sphere diameter D = 1.0. Try using a SimpleMixture instead of a SimpleFluid."
    dk = π  ./ (method.dr * (method.M))
    r = [i*method.dr for i = 0.5:(method.M-0.5)]
    k = [j*dk for j = 0.5:(method.M-0.5)]

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

    return OZSolution(r, k, gr, Sk, Cr, Ck, γmulr ./ r, γmulk ./ k, true, 0, 0.0, :exact)
end

function solve(system::SimpleFluid{1, T1, T2, HardSpheres{T3}}, ::PercusYevick, method::Exact) where {T1,T2,T3}
    # D = 1 by definition
    @assert system.potential.D == 1.0 "This method assumes that the hard sphere diameter D = 1.0"
    r = method.dr * (1:method.M) |> collect
    dk = π  ./ (method.dr * (method.M))
    k = dk * (1:method.M) |> collect
    mayer_f = find_mayer_f_function(system, r)
    elementtype = typeof(r[1] .* (system.kBT) .* system.ρ * (mayer_f[1]))
    mayer_f = elementtype.(mayer_f)
    fourierplan = get_fourier_plan(system, method, mayer_f)

    r, k = fourierplan.r, fourierplan.k

    ρ = system.ρ
    η = ρ

    # baxter
    Cr = find_Cr_1dPYSC_exact(η, r)

    # analytical tranform of the above:
    Ck = find_Ck_1dPYSC_exact(η, k)

    Sk = @. 1.0 + ρ * Ck / (1 - ρ * Ck)

    # γ is a smooth function, therefore it has a nicer transform
    γmulk = @. ρ * (Ck * k)^2 / (k - ρ * Ck * k)
    γmulr = similar(γmulk)
    inverse_fourier!(γmulr, γmulk, fourierplan)

    gr = @. Cr + γmulr / r + 1
    gr[r.<1.0] .= 0.0
    return OZSolution(r, k, gr, Sk, Cr, Ck, γmulr ./ r, γmulk ./ k, true, 0, 0.0, :exact)
end

function solve(system::SimpleFluid{5, T1, T2, HardSpheres{T3}}, ::PercusYevick, method::Exact) where {T1,T2,T3}
    # D = 1 by definition
    @assert system.potential.D == 1.0 "This method assumes that the hard sphere diameter D = 1.0"
    r = rand(method.M)
    mayer_f = find_mayer_f_function(system, r)
    elementtype = typeof(r[1] .* (system.kBT) .* system.ρ * (mayer_f[1]))
    mayer_f = elementtype.(mayer_f)
    fourierplan = get_fourier_plan(system, method, mayer_f)

    r, k = fourierplan.r, fourierplan.k

    ρ = system.ρ
    η = ρ * 8π^2 / (15 * 2^5)

    # baxter
    Cr = find_Cr_5dPYSC_exact(η, r)

    # analytical tranform of the above:
    Ck = find_Ck_5dPYSC_exact(η, k)

    Sk = @. 1.0 + ρ * Ck / (1 - ρ * Ck)

    # γ is a smooth function, therefore it has a nicer transform
    γmulk = @. ρ * (Ck * k)^2 / (k - ρ * Ck * k)
    γmulr = similar(γmulk)
    inverse_fourier!(γmulr, γmulk, fourierplan)

    gr = @. Cr + γmulr / r + 1
    gr[r.<1.0] .= 0.0
    return OZSolution(r, k, gr, Sk, Cr, Ck, γmulr ./ r, γmulk ./ k, true, 0, 0.0, :exact)
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
function solve(system::SimpleMixture{3, species, T1, T2, HardSpheres{T3}}, ::PercusYevick, method::Exact) where {species, T1, T2, T3<:AbstractMatrix}
    dk = π  ./ (method.dr * (method.M))
    r = [i*method.dr for i = 0.5:(method.M-0.5)]
    k = [j*dk for j = 0.5:(method.M-0.5)]
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
        @. gr[:, μ, ν] = (1.0 + γμν) * expnegβu
        @. Cr[:, μ, ν] = (expnegβu - 1.0) * (1 + γμν)
    end
    γr = gr .- Cr .- 1.0
    γk = Hk .- Ck
    return OZSolution(r, k, gr, Sk, Cr, Ck, γr, γk, true, 0, 0.0, :exact)
end

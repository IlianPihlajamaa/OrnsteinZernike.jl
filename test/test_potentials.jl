using Test
using StaticArrays

# ========= Helpers =========
# Central difference away from discontinuities
numdiff(f, r; h = sqrt(eps(float(r)))) = f(r + h) - f(r - h) |> x -> x / (2h)
isapprox0(x; atol=1e-12) = isapprox(x, 0.0; atol=atol)

# ========= Yukawa =========
@testset "Yukawa" begin
    A, κ = 1.7, 0.8
    p = Yukawa(A, κ)

    @test evaluate_potential(p, 2.0) ≈ A*exp(-κ*2.0)/2.0
    @test evaluate_potential_derivative(p, 2.0) ≈ (A*exp(-κ*2.0)) * (-(κ*2.0 + 1.0)) / (2.0^2)

    # Compare analytic derivative vs numerical
    for r in (0.7, 1.0, 3.0, 10.0)
        f = r -> evaluate_potential(p, r)
        @test evaluate_potential_derivative(p, r) ≈ numdiff(f, r) rtol=1e-6 atol=1e-10
    end

    # Broadcasting over vectors
    rs = collect(range(0.5, 5.0; length=7))
    us = evaluate_potential.((p,), rs)
    @test length(us) == length(rs)
    @test all(isfinite, us)

    # Mixture from charges (A_ij = q_i q_j)
    q = [1.0, -2.0, 0.5]
    pmix = Yukawa(q, κ)
    @test size(pmix.A) == (3,3)
    @test pmix.A ≈ q * q'

    # Mixture explicit A_ij
    Aij = @SMatrix [1.0 0.2; 0.2 0.5]
    pA = Yukawa(Matrix(Aij), κ)
    @test pA.A ≈ Aij
end

# ========= Gaussian Core Model =========
@testset "GaussianCore" begin
    ϵ, σ = 2.3, 1.4
    p = GaussianCore(ϵ, σ)

    @test evaluate_potential(p, 0.0) ≈ ϵ          # u(0) = ϵ
    @test evaluate_potential_derivative(p, 0.0) |> isapprox0  # du/dr at 0 should be 0

    # Analytic vs numerical derivative
    for r in (0.1, 0.7, 2.0, 5.0)
        f = r -> evaluate_potential(p, r)
        @test evaluate_potential_derivative(p, r) ≈ numdiff(f, r) rtol=1e-6 atol=1e-10
    end

    # Mixing rules: σ_ij additive, ϵ_ij geometric mean
    eps = [1.0, 0.5, 2.0]
    sig = [1.0, 1.2, 0.8]
    pmix = GaussianCore(eps, sig)
    ϵij_expect = sqrt.(eps * eps')
    σij_expect = (sig .+ sig') ./ 2
    @test pmix.ϵ ≈ ϵij_expect
    @test pmix.σ ≈ σij_expect

    # Broadcast sanity
    rs = range(0.0, 3.0; length=6) |> collect
    @test length(evaluate_potential.((p,), rs)) == length(rs)

    # No discontinuities
    @test isempty(discontinuities(p))
end

# ========= TabulatedPotential =========
@testset "TabulatedPotential" begin
    # Build a strictly linear function for exact interpolation/derivative checks
    rgrid = collect(range(0.5, 5.0; length=50))
    ugrid = @. 3.0 * rgrid + 2.0
    p_lin  = TabulatedPotential(rgrid, ugrid, :error)

    # Exact hits on nodes
    for i in (1, 10, 25, 50)
        @test evaluate_potential(p_lin, rgrid[i]) == ugrid[i]
    end

    # Interpolation within cells (should be exact for linear)
    @test evaluate_potential(p_lin, 1.23) ≈ 3.0*1.23 + 2.0 atol=1e-12

    # Derivative equals slope everywhere inside range (use right-slope at nodes by design)
    for r in (0.6, 1.0, 2.5, 4.9)
        @test evaluate_potential_derivative(p_lin, r) ≈ 3.0 atol=1e-12
    end

    # Extrapolation modes
    p_flat   = TabulatedPotential(rgrid, ugrid, :flat)
    p_linear = TabulatedPotential(rgrid, ugrid, :linear)
    @test evaluate_potential(p_flat, 0.4) == ugrid[1]
    @test evaluate_potential(p_flat, 6.0) == ugrid[end]
    @test evaluate_potential(p_linear, 0.4) ≈ ugrid[1] + (ugrid[2]-ugrid[1])/(rgrid[2]-rgrid[1])*(0.4 - rgrid[1])
    @test evaluate_potential(p_linear, 6.0) ≈ ugrid[end] + (ugrid[end]-ugrid[end-1])/(rgrid[end]-rgrid[end-1])*(6.0 - rgrid[end])

    # :error throws outside range
    @test_throws ErrorException evaluate_potential(p_lin, 0.49)
    @test_throws ErrorException evaluate_potential_derivative(p_lin, 5.01)

    # No discontinuities
    @test isempty(discontinuities(p_lin))
end

# ========= SquareWell =========
@testset "SquareWell" begin
    σ, ϵ, λ = 1.0, 2.0, 1.5
    p = SquareWell(σ, ϵ, λ)

    # Regions
    @test isinf(evaluate_potential(p, 0.99))
    @test evaluate_potential(p, 1.0) == -ϵ
    @test evaluate_potential(p, λ*σ - 1e-12) == -ϵ
    @test evaluate_potential(p, λ*σ + 1e-12) == 0.0

    # Derivative is zero away from discontinuities
    for r in (0.2, 1.1, 1.49, 2.0)
        @test evaluate_potential_derivative(p, r) == 0.0
    end

    # Discontinuities reported correctly
    discs = discontinuities(p)
    @test discs ≈ [σ, λ*σ]
end

# ========= Morse =========
@testset "Morse" begin
    ϵ, σ, α = 1.3, 1.1, 2.2
    p = Morse(ϵ, σ, α)

    # u has minimum at r = σ; du/dr = 0 there
    @test evaluate_potential_derivative(p, σ) |> isapprox0

    # Compare analytic derivative vs numerical at several r
    for r in (0.6, 1.1, 2.0, 3.5)
        f = r -> evaluate_potential(p, r)
        @test evaluate_potential_derivative(p, r) ≈ numdiff(f, r) rtol=1e-6 atol=1e-10
    end

    # Smooth -> no discontinuities
    @test isempty(discontinuities(p))
end

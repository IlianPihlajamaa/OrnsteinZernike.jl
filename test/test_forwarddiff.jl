using Test
using OrnsteinZernike
using ForwardDiff

@testset "ForwardDiff BomontBretonnet" begin
    dims = 3
    kBT = 1.0
    pot = HardSpheres(1.0)
    M = 1000
    method = NgIteration(M=M, dr=10.0 / M, verbose=false)
    closure = BomontBretonnet(f=0.5)

    pressure(ρ) = begin
        system = SimpleFluid(dims, ρ, kBT, pot)
        sol = solve(system, closure, method)
        compute_virial_pressure(sol, system)
    end

    ρ0 = 0.35
    autodiff = ForwardDiff.derivative(pressure, ρ0)

    δ = 5e-4
    finite_diff = (pressure(ρ0 + δ) - pressure(ρ0 - δ)) / (2δ)

    @test isapprox(autodiff, finite_diff; rtol=1e-5, atol=1e-6)
end

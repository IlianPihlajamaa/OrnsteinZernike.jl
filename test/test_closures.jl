function test_closure(closure::OrnsteinZernike.Closure)
    M = 1000
    dr = 10/M
    ρ = 0.6
    kBT = 1.0
    dims = 3

    base_pot = HardSpheres(1.0)
    pot = OrnsteinZernike.uses_renormalized_gamma(closure) ? AllShortRangeDivision(base_pot) : base_pot
    system = SimpleFluid(dims, ρ, kBT, pot)
    method0 = NgIteration(tolerance=10^-8, M=M, dr=dr, N_stages=2, max_iterations=1000, verbose=false)
    sol, = solve(system, closure, method0)
    @test !(any(isnan, sol.cr))  # Check for NaN in cr
    @test !(any(isnan, sol.gr))  # Check for NaN in gr
    @test !(any(isnan, sol.ck))  # Check for NaN in ck
    @test !(any(isnan, sol.Sk))  # Check for NaN in Sk
end


function test_closures()
    @testset "Closure Tests" begin
        for closure_type in subtypes(OrnsteinZernike.Closure)
            parentmodule(closure_type) === OrnsteinZernike || continue
            test_closure(closure_type())
        end
    end
end


test_closures()

@testset "Renormalized gamma support" begin
    closure = ZerahHansen()
    lj = LennardJones(1.0, 1.0)
    sys = SimpleFluid(3, 0.6, 1.0, lj)
    method = NgIteration(M=1000, dr=0.01, verbose=false)
    @test_throws ErrorException  sol = solve(sys, closure, method) 
end

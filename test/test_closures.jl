function test_closure(closure)
    M = 1000
    dr = 10/M
    ρ = 0.6
    kBT = 1.0
    dims = 3

    pot = HardSpheres(1.0)
    system = SimpleFluid(dims, ρ, kBT, pot)
    method0 = NgIteration(tolerance=10^-8, M=M, dr=dr, N_stages=2, max_iterations=1000, verbose=false)
    sol = solve(system, closure, method0)
    @test !(any(isnan, sol.cr))  # Check for NaN in cr
    @test !(any(isnan, sol.gr))  # Check for NaN in gr
    @test !(any(isnan, sol.ck))  # Check for NaN in ck
    @test !(any(isnan, sol.Sk))  # Check for NaN in Sk
end


function test_closures()
    @testset "Closure Tests" begin
        for closure in subtypes(OrnsteinZernike.Closure)
            test_closure(closure())
        end
    end
end


test_closures()
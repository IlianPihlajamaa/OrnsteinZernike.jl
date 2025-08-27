import Pkg; Pkg.activate(".")
using Revise
using OrnsteinZernike,  Plots, Dierckx
import Roots, ForwardDiff

function find_pressure_derivative(ρ, kBT, dims, pot, closure, method)

    function pressure(ρ)
        system = SimpleFluid(dims, ρ, kBT, pot)
        sol = solve(system, closure, method)
        p = compute_virial_pressure(sol, system)
        return p
    end

    dρ = sqrt(eps(ρ))

    # @time dpdρ = ForwardDiff.derivative(pressure, ρ)
    @time dpdρ2 =  (pressure(ρ+dρ) - pressure(ρ))/(dρ)
    @time dpdρ3 =  (pressure(ρ+dρ) - pressure(ρ-dρ))/(2dρ)
    @time dpdρ4 = (3kBT* (-72+π*ρ* (-60+π *ρ* (-18+π *ρ))))/(-6+π* ρ)^3
    @show dpdρ2, dpdρ3, dpdρ4

    return
end



println("Hard Spheres")
@profview for ρstar = [0.3]
    ρ = ρstar
    M = 100000
    dr = 10.0/M
    kBT = 1.0
    dims = 3
    pot = HardSpheres(1.0)
    closure = PercusYevick()
    method = NgIteration(M=M, dr=dr, verbose=false)
    find_pressure_derivative(ρ, kBT, dims, pot, closure, method)
end

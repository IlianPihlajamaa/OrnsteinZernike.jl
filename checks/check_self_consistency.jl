import Pkg; Pkg.activate(".")
using Revise
using OrnsteinZernike,  Plots
import Roots 

# function find_self_consistent_solution(ρ, kBT, method, dims, pot; lims=(0.1, 2.0))

#     function RY_inconsistency(ρ, α)
#         system1 = SimpleFluid(dims, ρ, kBT, pot)
#         sol1 = solve(system1, RogersYoung(α), method)
#         p1 = compute_virial_pressure(sol1, system1)

#         dρ = sqrt(eps(ρ))
#         system2 = SimpleFluid(dims, ρ+dρ, kBT, pot)
#         sol2 = solve(system2, RogersYoung(α), method)
#         p2 = compute_virial_pressure(sol2, system2)
#         dpdρ = (p2-p1)/dρ

#         χ = compute_compressibility(sol1, system1)
#         inconsistency = dpdρ/kBT - 1/(ρ*kBT*χ)
#         return inconsistency
#     end

#     func = α ->  RY_inconsistency(ρ, α)
#     α =  Roots.find_zero(func, lims, Roots.Bisection(), atol=0.0001)
#     system = SimpleFluid(dims, ρ, kBT, pot)
#     sol = solve(system, RogersYoung(α), method)
#     return system, sol, α
# end
# println("Hard Spheres")
# for ρstar = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.654]
#     ρ = ρstar*sqrt(2)
#     M = 1000
#     dr = 10.0/M
#     kBT = 1.0
#     dims = 3
#     pot = HardSpheres(1.0)
#     method = NgIteration(dr=dr, M=M, verbose=false)
#     system, sol, α = find_self_consistent_solution(ρ, kBT, method, dims, pot)
#     P = compute_virial_pressure(sol, system)/ρ/kBT - 1
#     gmax = maximum(sol.gr)
#     println("At ρ/√2 = $(ρstar), we find α = $(trunc(α,digits=2)), and βp/ρ - 1 = $(trunc(P,digits=4)).")
# end

# ## 12 power fluid
# println("1/r^12 fluid")

# for z = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.813]
#     ρ = z*sqrt(2)
#     M = 1000
#     dr = 20.0/M
#     kBT = 1.0
#     dims = 3 
#     Γ = (4π*sqrt(2)*z/3)^(4)
#     σ = 1.0
#     ϵ = 1.0

#     pot = PowerLaw(ϵ, σ, 12)
#     method = NgIteration(dr=dr, M=M, verbose=false)
#     system, sol, α = find_self_consistent_solution(ρ, kBT, method, dims, pot)
#     P = compute_virial_pressure(sol, system)/ρ/kBT - 1
#     gmax = maximum(sol.gr)
#     println("At z = $(z), we find Γ = $(trunc(Γ,digits=2)) α = $(trunc(α,digits=2)), and βp/ρ - 1 = $(trunc(P,digits=4)).")
# end

# println("1/r^9 fluid")

# for z = [0.1, 0.25, 0.5, 0.943]
#     ρ = z*sqrt(2)
#     M = 1000
#     dr = 20.0/M
#     kBT = 1.0
#     dims = 3 
#     Γ = (4π*sqrt(2)*z/3)^(3)
#     σ = 1.0
#     ϵ = 1.0

#     pot = PowerLaw(ϵ, σ, 9)
#     method = NgIteration(dr=dr, M=M, verbose=false)
#     system, sol, α = find_self_consistent_solution(ρ, kBT, method, dims, pot)
#     P = compute_virial_pressure(sol, system)/ρ/kBT - 1
#     gmax = maximum(sol.gr)
#     println("At z = $(z), we find Γ = $(trunc(Γ,digits=2)) α = $(trunc(α,digits=2)), and βp/ρ - 1 = $(trunc(P,digits=4)).")
# end

# println("1/r^6 fluid")

# for z = [0.1, 0.25, 0.5, 1.0, 1.54]
#     ρ = z*sqrt(2)
#     M = 2000
#     dr = 20/M
#     kBT = 1.0
#     dims = 3
#     Γ = (4π*sqrt(2)*z/3)^(2)
#     σ = 1.0
#     ϵ = 1.0

#     pot = PowerLaw(ϵ, σ, 6)
#     method = NgIteration(dr=dr, M=M, verbose=false)
#     system, sol, α = find_self_consistent_solution(ρ, kBT, method, dims, pot)
#     P = compute_virial_pressure(sol, system)/ρ/kBT - 1
#     gmax = maximum(sol.gr)
#     println("At z = $(z), we find Γ = $(trunc(Γ,digits=2)) α = $(trunc(α,digits=2)), and βp/ρ - 1 = $(trunc(P,digits=4)).")
# end


# println("Yukawa")
# #J. Chem. Phys. 128, 184507 共2008兲
# for ϕ = [0.2, 0.3, 0.35, 0.4, 0.45, 0.5]
#     ρ = ϕ*6/pi
#     M = 10000
#     dr = 10.0/M
#     kBT = 1.0
#     dims = 3
#     method = NgIteration(dr = dr, M = M, verbose=false, max_iterations=10^6, N_stages=2)
#     pot = CustomPotential((r,p)->500*exp(-4r)/r, nothing)
#     system, sol, α = find_self_consistent_solution(ρ, kBT, method, dims, pot; lims=(0.9, 2.0))
#     P = compute_virial_pressure(sol, system)/ρ/kBT - 1
#     println("At ϕ = $(ϕ), α = $(trunc(α,digits=4)), and βp/ρ = $(trunc(P,digits=4)).")
# end


using Optim
function find_self_consistent_solution_ERY(ρ, kBT, method, dims, pot; x0=[0.2, 0.2])
    function ERY_inconsistency(ρ, α2, a2)
        α = sqrt(α2)
        a = sqrt(a2)
        closure = ExtendedRogersYoung(α, a)
        dρ = ρ*0.001
        dkBT = sqrt(eps(kBT))

        system0 = SimpleFluid(dims, ρ-dρ, kBT, pot)
        sol0 = solve(system0, closure, method)
        system1 = SimpleFluid(dims, ρ, kBT, pot)
        sol1 = solve(system1, closure, method)
        system2 = SimpleFluid(dims, ρ+dρ, kBT, pot)
        sol2 = solve(system2, closure, method)
        system3 = SimpleFluid(dims, ρ, kBT+dkBT, pot)
        sol3 = solve(system3, closure, method)

        p0 = compute_virial_pressure(sol0, system0)
        p2 = compute_virial_pressure(sol2, system2)

        ρU0 = (ρ-dρ)*compute_excess_energy(sol0, system0)
        ρU1 = ρ*compute_excess_energy(sol1, system1)
        ρU2 = (ρ+dρ)*compute_excess_energy(sol2, system2)


        dpdρ = (p2-p0)/(2dρ)
        d2ρUdρ2 = (ρU2 + ρU0 - 2ρU1)/(dρ^2)

        ĉ0_1 = OrnsteinZernike.get_ĉ0(sol1, system1)
        ĉ0_3 = OrnsteinZernike.get_ĉ0(sol3, system3)
        dĉ0dkBT = (ĉ0_3-ĉ0_1)/dkBT
        dĉ0dβ = - kBT^2 * dĉ0dkBT

        χ = compute_compressibility(sol1, system1)
        inconsistency1 = dpdρ/kBT - 1/(ρ*kBT*χ)

        inconsistency2 = d2ρUdρ2 + dĉ0dβ

        inconsistency = sum([inconsistency1, inconsistency2].^2)
        @show inconsistency, α, a, d2ρUdρ2
        return inconsistency
    end

    func = x ->  ERY_inconsistency(ρ, x[1]^2, x[2]^2)
    α =  optimize(func, [x0...]+rand(2)./10, NelderMead(), Optim.Options(x_tol=0.0001, f_tol=10^-6)).minimizer
    system = SimpleFluid(dims, ρ, kBT, pot)
    sol = solve(system, ExtendedRogersYoung(α[1], α[2]), method)
    return system, sol, α
end


x0s = [(0.323, 0.177), (0.375, 0.213), (0.398, 0.231), (0.423, 0.243), (0.446, 0.255), (0.467, 0.270)]
for (i,ϕ) = enumerate([0.2, 0.3, 0.35, 0.4, 0.45, 0.5])
    ρ = ϕ*6/pi
    M = 10000
    dr = 10.0/M
    kBT = 1.0
    dims = 3
    method = NgIteration(dr = dr, M = M, verbose=false, max_iterations=10^6, N_stages=2)
    pot = CustomPotential((r,p)->500*exp(-4r)/r, nothing)
    system, sol, α = find_self_consistent_solution_ERY(ρ, kBT, method, dims, pot; x0=x0s[i])
    P = compute_virial_pressure(sol, system)/ρ/kBT - 1
    println("At ϕ = $(ϕ), (α,a) = $(trunc(α[1],digits=4)), $(trunc(α[2],digits=4)), and βp/ρ = $(trunc(P,digits=4)).")
end
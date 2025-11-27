using NLsolve, OrnsteinZernike

function find_self_consistent_solution_Lee(ρ, kBT, M, dr, dims, pot)
    function find_inconsistency(ρ, params)
        ζ, ϕ, α = params
        @show ζ, ϕ, α
        closure = Lee(ζ, ϕ, α, ρ)
        system1 = SimpleFluid(dims, ρ, kBT, pot)
        method = NgIteration(M=M, dr=dr, verbose=false)
        sol1, = solve(system1, closure, method)
        p1 = compute_virial_pressure(sol1, system1)

        dρ = sqrt(eps(ρ))
        system2 = SimpleFluid(dims, ρ+dρ, kBT, pot)
        sol2, = solve(system2, closure, method)
        p2 = compute_virial_pressure(sol2, system2)

        dpdρ = (p2-p1)/dρ

        χ = compute_compressibility(sol1, system1)
        inconsistency1 = dpdρ - 1/(ρ*χ)

        βu = OrnsteinZernike.evaluate_potential(system1.potential, sol1.r)
        mayer_f = OrnsteinZernike.find_mayer_f_function.((system1,), βu)

        gammastar0 = sol1.gr[1] .- sol1.cr[1] .- 1 .+ ρ/2*mayer_f[1]
        gamma0 = sol1.gr[1] .- sol1.cr[1] .- 1 
        
        η = ρ/6*π
        grmax = (1-0.5η)/(1-η)^3

        dBdgamma0 = -gammastar0*ζ - (α^2*gammastar0^3*ϕ*ζ)/(2*(1*+ α*gammastar0)^2) + (3*α*gammastar0^2*ϕ*ζ)/(2*(1 + α*gammastar0))
        inconsistency2 = dBdgamma0 - 1 + 1/grmax

        b = OrnsteinZernike.bridge_function(closure, sol1.r, mayer_f, sol1.gr .- sol1.cr .- 1.0)
        mu = (8η-9η^2+3η^3)/(1-η)^3
        inconsistency3 = b[1] + gamma0 - mu
        return [inconsistency1,inconsistency2, inconsistency3]
    end
    @show find_inconsistency(ρ, [1.2041, 0.9962, 1.0])
    func = α ->  find_inconsistency(ρ, α)
    s = nlsolve(func,  [1.2041, 0.9962, 1.0], xtol=0.001, factor=0.1)
    @show s
    params = s.zero
    ζ, ϕ, α = params
    system = SimpleFluid(dims, ρ, kBT, pot)
    method = NgIteration(M=M, dr=dr, verbose=false)
    closure = Lee(ζ, ϕ, α, ρ)
    sol, = solve(system, closure, method)
    return system, sol, params
end

ρ = 0.9
sys, sol, params = find_self_consistent_solution_Lee(ρ, 1.0, 10^4, 0.001, 3,  HardSpheres(1.0))
ζ, ϕ, α= params
closure = Lee(ζ, ϕ, α, ρ)

βu = OrnsteinZernike.evaluate_potential(sys.potential, sol.r)
mayer_f = OrnsteinZernike.find_mayer_f_function.((sys,), βu)

br = OrnsteinZernike.bridge_function(closure, sol.r, mayer_f, sol.gr - sol.cr .- 1)
@show br[1]
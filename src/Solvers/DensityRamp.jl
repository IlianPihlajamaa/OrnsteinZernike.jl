function solve(system::SimpleLiquid, closure::Closure, method::DensityRamp)
    densities = method.densities
    ρtarget = system.ρ
    system.ρ = densities[1]
    if method.verbose
        println("\nSolving the system at ρ = $(densities[1]).\n")
    end
    sol0 = solve(system, closure, method.method)
    γ_old = @. sol0.gr - one(eltype(sol0.gr)) - sol0.cr
    for i = (firstindex(densities)+1):lastindex(densities)
        if method.verbose
            println("\nSolving the system at ρ = $(densities[i]).\n")
        end
        system.ρ = densities[i]
        sol = solve(system, closure, method.method, init=γ_old)
        @. γ_old =  sol.gr - one(eltype(sol.gr)) - sol.cr
    end
    if method.verbose
        println("\nSolving the system at ρ = $(ρtarget).\n")
    end
    system.ρ = ρtarget
    sol = solve(system, closure, method.method, init=γ_old)
    return sol
end 

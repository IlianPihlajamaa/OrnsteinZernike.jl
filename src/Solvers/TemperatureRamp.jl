function solve(system::SimpleUnchargedSystem, closure::Closure, method::TemperatureRamp)
    temperatures = method.temperatures
    kBTtarget = system.kBT
    sols = []
    system.kBT = temperatures[1]
    if method.verbose
        println("\nSolving the system at kBT = $(system.kBT).\n")
    end
    sol0 = solve(system, closure, method.method)
    push!(sols, sol0)
    γ_old = @. sol0.gr - one(eltype(sol0.gr)) - sol0.cr
    for i = (firstindex(temperatures)+1):lastindex(temperatures)
        system.kBT = temperatures[i]
        if method.verbose
            println("\nSolving the system at kBT = $(system.kBT).\n")
        end
        sol = solve(system, closure, method.method, init=γ_old)
        push!(sols, sol)
        @. γ_old =  sol.gr - one(eltype(sol.gr)) - sol.cr
    end
    if method.verbose
        println("\nSolving the system at kBT = $(kBTtarget).\n")
    end
    system.kBT = kBTtarget
    sol = solve(system, closure, method.method, init=γ_old)
    push!(sols, sol)
    return sols
end 

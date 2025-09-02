function solve(system::System, closure::Closure, method::TemperatureRamp; kwargs...)
    temperatures = method.temperatures
    kBTtarget = kBT_of(system)
    sols = []
    base_of(system).kBT = temperatures[1]
    if method.verbose
        println("\nSolving the system at kBT = $(kBT_of(system)).\n")
    end
    sol0 = solve(system, closure, method.method)
    push!(sols, sol0)
    γ_old = @. sol0.gr - one(eltype(sol0.gr)) - sol0.cr
    for i = (firstindex(temperatures)+1):lastindex(temperatures)
        base_of(system).kBT = temperatures[i]
        if method.verbose
            println("\nSolving the system at kBT = $(kBT_of(system)).\n")
        end
        sol = solve(system, closure, method.method, gamma_0=γ_old; kwargs...)
        push!(sols, sol)
        @. γ_old =  sol.gr - one(eltype(sol.gr)) - sol.cr
    end
    if method.verbose
        println("\nSolving the system at kBT = $(kBT_of(system)).\n")
    end
    base_of(system).kBT = kBTtarget
    sol = solve(system, closure, method.method, gamma_0=γ_old; kwargs...)
    push!(sols, sol)
    return sols
end 

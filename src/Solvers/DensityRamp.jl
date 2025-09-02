function recast_γ(γ::Array{T, 3}) where {T}
    Ns = size(γ, 2)
    Nk = size(γ, 1)
    γ2 = zeros(SMatrix{Ns, Ns, T, Ns^2}, Nk)
    for i = 1:Nk
        γ2[i] = γ[i, :, :]
    end
    return γ2
end

recast_γ(γ::Array{T, 1}) where {T} = γ


function solve(system::System, closure::Closure, method::DensityRamp; kwargs...)
    densities = method.densities
    ρtarget = ρ_of(system)
    sols = []
    base_of(system).ρ = densities[1]
    if method.verbose
        println("\nSolving the system at ρ = $(densities[1]).\n")
    end
    sol0 = solve(system, closure, method.method)
    push!(sols, sol0)
    γ_old = @. sol0.gr - one(eltype(sol0.gr)) - sol0.cr
    for i = (firstindex(densities)+1):lastindex(densities)
        if method.verbose
            println("\nSolving the system at ρ = $(densities[i]).\n")
        end
        base_of(system).ρ = densities[i]
        γ_old2 = recast_γ(γ_old)
        sol = solve(system, closure, method.method, gamma_0=γ_old2; kwargs...)
        push!(sols, sol)
        @. γ_old =  sol.gr - one(eltype(sol.gr)) - sol.cr
    end
    if method.verbose
        println("\nSolving the system at ρ = $(ρtarget).\n")
    end
    base_of(system).ρ = ρtarget
    γ_old2 = recast_γ(γ_old)
    sol = solve(system, closure, method.method, gamma_0=γ_old2; kwargs...)
    push!(sols, sol)
    return sols
end 

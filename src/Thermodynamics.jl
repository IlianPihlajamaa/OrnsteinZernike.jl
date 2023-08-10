## Single component

function compute_excess_energy(sol::OZSolution, system::SimpleLiquid{3, 1, T1, T2, P}) where {T1, T2, P}
    r = sol.r
    dr = r[1] - r[2]
    gr = sol.gr
    ρ = system.ρ
    u = evaluate_potential.(Ref(system.potential), r)
    E = zero(eltype(gr))
    for i in eachindex(r, gr, u)
        if isinf(u[i]) 
            if iszero(gr[i])
                continue
            else
                error("Inconsistency: g(r) is finite where u(r) is infinite.")
            end
        end
        E += gr[i]*u[i]*r[i]^2
    end
    E *= 2π*ρ*dr
    return E
end

function compute_virial_pressure(sol::OZSolution, system::SimpleLiquid{3, 1, T1, T2, P}) where {T1, T2, P}
    r = sol.r
    ρ = system.ρ
    kBT = system.kBT
    β = 1/kBT
    dr = r[1] - r[2]
    gr = sol.gr
    potential = system.potential
    dudr = evaluate_potential_derivative.(Ref(potential), r)
    p = zero(eltype(gr))
    for i in eachindex(r, gr, dudr)
        if isinf(dudr[i]) 
            if iszero(gr[i])
                continue
            else
                error("Inconsistency: g(r) is finite where du(r)/dr is infinite.")
            end
        end
        p += gr[i]*dudr[i]*r[i]^3
    end
    p = kBT*ρ - 2/3 * π * ρ^2 * dr * p

    ## now add the terms for the discontinuities
    discs = discontinuities(system.potential)
    
    for discontinuity in discs
        eleft = exp(-β*evaluate_potential(potential, prevfloat(discontinuity)))
        eright = exp(-β*evaluate_potential(potential, nextfloat(discontinuity)))
        index = searchsortedfirst(r, discontinuity) # first index >= disc

        rmin1 = r[index-1]
        rmin2 = r[index-2]
        ymin1 = exp(β*evaluate_potential(potential, rmin1))*gr[index-1]
        ymin2 = exp(β*evaluate_potential(potential, rmin2))*gr[index-2]
        yleft = ymin1 + (discontinuity - rmin1) * (ymin2 - ymin1) / (rmin2 - rmin1)
        
        r1 = r[index+1]
        r2 = r[index+2]
        y1 = exp(β*evaluate_potential(potential, r1))*gr[index+1]
        y2 = exp(β*evaluate_potential(potential, r2))*gr[index+2]
        yright = y1 + (discontinuity - r1) * (y2 - y1) / (r2 - r1)
        if isfinite(yleft) && isfinite(yright)
            ydisc = (yleft + yright)/2
            if abs(yleft - yright) > 0.1
                error("This is weird, the cavity disctribution function looks discontinuous")
            end
            p += (2π*ρ^2)/(3β)*discontinuity^3*(eright-eleft)*ydisc
        elseif isfinite(yleft)
            ydisc = yleft
            p += (2π*ρ^2)/(3β)*discontinuity^3*(eright-eleft)*ydisc
        elseif isfinite(yright)
            ydisc = yright
            p += (2π*ρ^2)/(3β)*discontinuity^3*(eright-eleft)*ydisc
        else
            @show yleft, yright
            error("This cannot happen.")
        end
    end
    return p
end

function compute_compressibility(sol::OZSolution, system::SimpleLiquid{dims, 1, T1, T2, P}) where {dims, T1, T2, P}
    ρ = system.ρ
    ĉ = sol.ck
    k = sol.k
    dĉdk = (ĉ[2]-ĉ[1])/(k[2]-k[1])
    ĉ0 = ĉ[1] - k[1] * dĉdk
    invχ = one(ρ) - ρ * ĉ0
    return one(ρ)/invχ
end

## Mixtures

function compute_excess_energy(sol::OZSolution, system::SimpleLiquid{3, species, T1, T2, P}) where {species, T1, T2, P}
    error("Not implemented yet")
    r = sol.r
    ρ = system.ρ
    ρ0 = sum(ρ.diag)
    x = ρ.diag/ρ0
    dr = r[1] - r[2]
    gr = sol.gr
    u = evaluate_potential.(Ref(system.potential), r)
    E = zero(eltype(gr))
    for i in eachindex(r, gr, u)
        if isinf(u[i]) 
            if iszero(gr[i])
                continue
            else
                error("Inconsistency: g(r) is finite where u(r) is infinite.")
            end
        end
        E += sum((x*x') .* gr[i] .* u[i])*r[i]^2
    end
    E *= 2π*ρ*dr
    return E
end

function compute_virial_pressure(sol::OZSolution, system::SimpleLiquid{3, species, T1, T2, P}) where {species, T1, T2, P}
    r = sol.r
    ρ = system.ρ
    ρ0 = sum(ρ.diag)
    x = ρ.diag/ρ0
    kBT = system.kBT
    dr = r[1] - r[2]
    gr = sol.gr
    dudr = evaluate_potential_derivative.(Ref(system.potential), r)
    p = zero(eltype(gr))
    for i in eachindex(r, gr, dudr)
        if isinf(dudr[i]) 
            if iszero(gr[i])
                continue
            else
                error("Inconsistency: g(r) is finite where du(r)/dr is infinite.")
            end
        end
        p += sum((x*x') .* gr[i] .* dudr[i])*r[i]^3 
    end
    p = kBT*ρ - 2/3 * π * ρ^2 * dr * p
    return p
end

function compute_compressibility(sol::OZSolution, system::SimpleLiquid{dims, species, T1, T2, P}) where {dims, species, T1, T2, P}
    ρ = system.ρ
    ρ0 = sum(ρ.diag)
    x = ρ.diag/ρ0
    ĉ = sol.ck
    T = typeof(ρ)
    k = sol.k
    dcdk = (c[2]-c[1])/(k[2]-k[1])
    ĉ0 = ĉ[1] - k[1] * dcdk
    invχ = one(T) - ρ * sum((x*x') .* ĉ0)
    return one(T)/invχ
end
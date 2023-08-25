## Single component

"""
    compute_excess_energy(sol::OZSolution, system::SimpleLiquid)

Computes the excess energy per particle Eₓ, such that E = (dims/2*kBT + Eₓ)*N.

uses the formula Eₓ = 1/2 ρ ∫dr g(r) u(r) for single component systems
and Eₓ = 1/2 ρ Σᵢⱼ xᵢxⱼ ∫dr gᵢⱼ(r) uᵢⱼ(r) for mixtures. Here x is the concentration fraction xᵢ=ρᵢ/sum(ρ).

"""
function compute_excess_energy(sol::OZSolution, system::SimpleLiquid{dims, species, T1, T2, P}) where {dims,species, T1, T2, P}
    r = sol.r
    ρ = system.ρ
    ρ0 = sum(ρ)
    x = get_concentration_fraction(system)
    dr = r[1] - r[2]
    gr = sol.gr
    u = evaluate_potential.(Ref(system.potential), r)
    E = zero(eltype(eltype(gr)))

    rpow = dims-1
    sphere_surface = surface_N_sphere(dims)

    for i in eachindex(r)
        if any( isinf.(u[i]) .& .!(iszero.(gr[i])))
            error("Inconsistency: g(r) is finite where u(r) is infinite.")
        end
        if any(isinf.(u[i]))
            continue
        end
        E += sum((x*x') .* gr[i] .* u[i])*r[i]^rpow
    end
    E *= sphere_surface*ρ0*dr/2
    return E
end


function find_left_and_right_lim_y(potential, β, gr, r::AbstractArray, r0::Number)
    index = searchsortedfirst(r, r0) # first index >= disc
    rmin1 = r[index-1]
    rmin2 = r[index-2]
    ymin1 = exp.(β*evaluate_potential(potential, rmin1)).*gr[index-1]
    ymin2 = exp.(β*evaluate_potential(potential, rmin2)).*gr[index-2]
    yleft = ymin1 + (r0 - rmin1) * (ymin2 - ymin1) / (rmin2 - rmin1)
    
    r1 = r[index+1]
    r2 = r[index+2]
    y1 = exp.(β*evaluate_potential(potential, r1)).*gr[index+1]
    y2 = exp.(β*evaluate_potential(potential, r2)).*gr[index+2]
    yright = y1 + (r0 - r1) * (y2 - y1) / (r2 - r1)
    return yleft, yright
end


function find_de_mul_y0(discontinuity, β, potential, gr, r)
    eleft = exp.(-β*evaluate_potential(potential, prevfloat(discontinuity)))
    eright = exp.(-β*evaluate_potential(potential, nextfloat(discontinuity)))
    de = eright - eleft
    yleft, yright = find_left_and_right_lim_y(potential, β, gr, r, discontinuity)

    if all(isfinite.(yleft)) && all(isfinite.(yright))
        ydisc = (yleft + yright)/2
        if abs(yleft - yright) > 0.1
            error("This is weird, the cavity distribution function looks discontinuous")
        end
        return de.*ydisc
    elseif all(isfinite.(yleft))
        ydisc = yleft
        return de.*ydisc
    elseif all(isfinite.(yright))
        ydisc = yright
        return de.*ydisc
    else
        @show yleft, yright
        error("This cannot happen. File an issue")
    end
end


# surface of sphere embedded in n dimensions d=3->4pi
function surface_N_sphere(n)
    return 2π^(n/2)/gamma(n/2)
end


"""
    compute_virial_pressure(sol::OZSolution, system::SimpleLiquid)

Computes the pressure via the virial route

uses the formula p = kBTρ - 1/6 ρ^2 ∫dr r g(r) u'(r) for single component systems
and p =  kBT Σᵢρᵢ - 1/6 Σᵢⱼ ρᵢρⱼ ∫dr r gᵢⱼ(r) u'ᵢⱼ(r) for mixtures.

It handles discontinuities in the interaction potential analytically if `discontinuities(potential)` is defined.
For additional speed/accuracy define a method of `evaluate_potential_derivative(potential, r::Number)` that analytically computes du/dr. 
By default this is done using finite differences.
"""
function compute_virial_pressure(sol::OZSolution, system::SimpleLiquid{dims, Nspecies, T1, T2, P}) where {dims, Nspecies, T1, T2, P}
    r = sol.r
    ρ = system.ρ
    ρ0 = sum(ρ)
    x = get_concentration_fraction(system)
    kBT = system.kBT
    β = 1/kBT
    dr = r[1] - r[2]
    gr = sol.gr
    potential = system.potential
    dudr = evaluate_potential_derivative(potential, r)
    rpow = dims
    p = zero(eltype(gr))
    for i in eachindex(r)
        p += sum((x*x') .* gr[i] .* dudr[i])*r[i]^rpow 
    end

    sphere_surface = surface_N_sphere(dims)
    p = kBT*ρ0 - sphere_surface/6 * ρ0^2 * dr * p
    ## now add the terms for the discontinuities

    discs = unique(discontinuities(system.potential))
    for discontinuity in discs
        dey0 = find_de_mul_y0(discontinuity, β, potential, gr, r)
        dp = (sphere_surface*ρ0^2)/(6β)*discontinuity^rpow*sum((x*x') .* dey0)
        p += dp
    end

    return p
end

function get_concentration_fraction(system)
    ρ = system.ρ
    if ρ isa AbstractArray
        return ρ.diag / sum(ρ.diag)
    elseif ρ isa Number
        return one(ρ)
    end
    error("Unreachable code: file an issue!")
end



"""
    compute_compressibility(sol::OZSolution, system::SimpleLiquid)

Computes the isothermal compressibility χ of the system

uses the formula 1/χ = 1 - ρ ĉ(k=0) for single component systems and
1/χ = 1 - ρ Σᵢⱼ xᵢxⱼ ĉᵢⱼ(k=0) for mixtures. 
"""
function compute_compressibility(sol::OZSolution, system::SimpleLiquid{dims, species, T1, T2, P}) where {dims, species, T1, T2, P}
    ρ = system.ρ
    kBT = system.kBT

    ρ0 = sum(ρ)
    x = get_concentration_fraction(system)
    ĉ = sol.ck
    T = typeof(ρ0)
    k = sol.k
    dcdk = (ĉ[2]-ĉ[1])/(k[2]-k[1])
    ĉ0 = ĉ[1] - k[1] * dcdk
    invρkBTχ = one(T) - ρ0 * sum((x*x') .* ĉ0)
    χ = (one(T)/invρkBTχ)/(kBT*ρ0)
    return χ
end
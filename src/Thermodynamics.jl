## Single component


"""
    integrate(x::AbstractVector, y::AbstractVector, ::SimpsonFast)

Use Simpson's rule on an irregularly spaced grid x.
"""
function simpsons_rule(x::AbstractVector, y::AbstractVector)
    length(x) == length(y) || error("x and y vectors must be of the same length!")
    length(x) ≥ 2 || error("vectors must contain at least 3 elements")
    N = length(x)
    retval = zero(eltype(y))*x[1]
    yjp2 = y[1]
    for i in 0:floor(Int64, (N-1)/2)-1
        j = 2i+1
        jp1 = j+1
        jp2 = jp1+1

        dxj = x[jp1] - x[j]
        dxjp1 = x[jp2] - x[jp1]

        yj = yjp2
        yjp1 = y[jp1]
        yjp2 = y[jp2]

        term1 = (2 - dxjp1/dxj) * yj
        term2 = (dxjp1 + dxj)^2 / (dxjp1*dxj) * yjp1
        term3 = (2 - dxj/dxjp1) * yjp2

        retval += (dxj+dxjp1)*(term1+term2+term3)/6
    end
    if iseven(N)
        dxNm1 = x[N] - x[N-1]
        dxNm2 = x[N-1] - x[N-2]
        retval += (2dxNm1^2+3dxNm1*dxNm2) / (6(dxNm1+dxNm2)) * y[N]
        retval += (dxNm1^2+3dxNm1*dxNm2) / (6dxNm2) * y[N-1]
        retval -= (dxNm1^3) / (6dxNm2*(dxNm2+dxNm1)) * y[N-2]
    end
    return retval
end


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
    gr = sol.gr
    u = evaluate_potential.(Ref(system.potential), r)
    E = zero(eltype(eltype(gr)))

    rpow = dims-1
    sphere_surface = surface_N_sphere(dims)
    fraction_matrix = (x*x')

    for s1 = axes(gr,2)
        for s2 = axes(gr,3)
            integrand = gr[:, s1, s2] .* getindex.(u, s1, s2) .* r[:] .^ rpow 
            for i in eachindex(integrand) # if u is inf or gr is very small, set contribution to zero 
                if isnan(integrand[i]) || gr[i, s1, s2] < 10^-6
                    integrand[i] = zero(eltype(integrand))
                end
            end
            integral = simpsons_rule(r, integrand)
            E += fraction_matrix[s1,s2] * integral
        end
    end

    E *= sphere_surface*ρ0/2
    return E
end


function find_left_and_right_lim_y(potential, β, gr, r::AbstractArray, r0::Number)
    index = searchsortedfirst(r, r0) # first index >= disc
    rmin1 = r[index-1]
    rmin2 = r[index-2]
    ymin1 = exp.(β*evaluate_potential(potential, rmin1)).*gr[index-1, :, :]
    ymin2 = exp.(β*evaluate_potential(potential, rmin2)).*gr[index-2, :, :]
    yleft = ymin1 + (r0 - rmin1) * (ymin2 - ymin1) / (rmin2 - rmin1)
    
    r1 = r[index+1]
    r2 = r[index+2]
    y1 = exp.(β*evaluate_potential(potential, r1)).*gr[index+1, :, :]
    y2 = exp.(β*evaluate_potential(potential, r2)).*gr[index+2, :, :]
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

uses the formula p = kBTρ - 1/(2d) ρ^2 ∫dr r g(r) u'(r) for single component systems
and p =  kBT Σᵢρᵢ - 1/(2d) Σᵢⱼ ρᵢρⱼ ∫dr r gᵢⱼ(r) u'ᵢⱼ(r) for mixtures.

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
    gr = sol.gr
    potential = system.potential
    dudr = evaluate_potential_derivative(potential, r)
    rpow = dims
    p1 = zero(eltype(gr))
    fraction_matrix = (x*x')
    for s1 = axes(gr,2)
        for s2 = axes(gr,3)
            integrand = gr[:, s1, s2] .* getindex.(dudr, s1, s2) .* r[:] .^ rpow 
            for i in eachindex(integrand) # if u is inf or gr is very small, set contribution to zero 
                if isnan(integrand[i]) || gr[i, s1, s2] < 10^-6
                    integrand[i] = zero(eltype(integrand))
                end
            end
            integral = simpsons_rule(r, integrand)
            p1 += fraction_matrix[s1,s2]*integral
        end
    end
    sphere_surface = surface_N_sphere(dims)
    p = kBT*ρ0 - sphere_surface/(2*dims) * ρ0^2 * p1

    ## now add the terms for the discontinuities

    discs = unique(discontinuities(system.potential))
    for discontinuity in discs
        dey0 = find_de_mul_y0(discontinuity, β, potential, gr, r)
        dp = (sphere_surface*ρ0^2)/((2*dims)*β)*discontinuity^rpow*sum((x*x') .* dey0)
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

_eachslice(a;dims=1) = eachslice(a,dims=dims)
_eachslice(a::Vector;dims=1) = a


"""
    compute_compressibility(sol::OZSolution, system::SimpleLiquid)

Computes the isothermal compressibility χ of the system

uses the formula 1/ρkBTχ = 1 - ρ ĉ(k=0) for single component systems and
1/ρkBTχ = 1 - ρ Σᵢⱼ ĉᵢⱼ(k=0) for mixtures. 
Eq. (3.6.16) in Hansen and McDonald
"""
function compute_compressibility(sol::OZSolution, system::SimpleLiquid{dims, species, T1, T2, P}) where {dims, species, T1, T2, P}
    ρ = system.ρ
    kBT = system.kBT
    x = get_concentration_fraction(system)
    ρ0 = sum(ρ)
    T = typeof(ρ0)
    ĉ0 = get_ĉ0(sol, system)
    invρkBTχ = one(T) - ρ0 * sum((x*x') .* ĉ0)
    χ = (one(T)/invρkBTχ)/(kBT*ρ0)
    return χ
end

function get_ĉ0(sol::OZSolution, ::SimpleLiquid{dims, 1, T1, T2, P}) where {dims, T1, T2, P}
    spl = Spline1D(sol.r, sol.cr[:, 1, 1].*sol.r.^(dims-1))
    ĉ0 = surface_N_sphere(dims)*integrate(spl, zero(eltype(sol.r)), maximum(sol.r))
    return ĉ0
end

function get_ĉ0(sol::OZSolution, ::SimpleLiquid{dims, species, T1, T2, P}) where {dims, species, T1, T2, P}
    ĉ0 = zeros(eltype(eltype(sol.ck)), species, species)
    for i = 1:species
        for j = 1:species
            spl = Spline1D(sol.r, sol.cr[:, i, j].*sol.r.^(dims-1))
            ĉ0[i,j] = surface_N_sphere(dims)*integrate(spl, zero(eltype(sol.r)), maximum(sol.r))
        end
    end
    return ĉ0
end



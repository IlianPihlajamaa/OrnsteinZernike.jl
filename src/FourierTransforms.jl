struct My3DPlan{T1,T2,T3,T4, T5, T6}
    plan::T1
    r::T5
    k::T6
    dr::T2
    dk::T3
    M::T4
end

struct My1DPlan{T2,T3,T4, T5, T6}
    r::T5
    k::T6
    dr::T2
    dk::T3
    M::T4
end

function get_fourier_plan(::SimpleLiquid{1, species, T1, T2, P}, method, F) where {species, T1, T2, P}
    M = length(F)
    dr = method.dr
    dk =  π/((M+0.5)*dr)
    r = [i*dr for i = 0.5:(M-0.5)]
    k = [j*dk for j = 0.5:(M-0.5)]
    return My1DPlan(r, k, dr, dk, M)
end


function get_fourier_plan(::SimpleLiquid{3, species, T1, T2, P}, method, F) where {species, T1, T2, P}
    plan =  find_fourier_plan_3d(F)
    M = length(F)
    dr = method.dr
    dk =  π/(M+1)/dr
    r = collect((1:M)*dr)
    k = collect((1:M)*dk)
    return My3DPlan(plan, r, k, dr, dk, M)
end

function get_fourier_plan(system::SimpleLiquid{ndims, species, T1, T2, P}, method, F) where {ndims, species, T1, T2, P}
    dr = method.dr
    return find_fourier_plan_nd(system, F, dr)
end

function find_fourier_plan_3d(F::Vector{T}) where T
    if T <: Number
        plan = FFTW.plan_r2r!(copy(F), FFTW.RODFT00; flags=FFTW.ESTIMATE)
        return plan
    elseif T <: AbstractArray
        F2 = copy(F)
        Nspecies = size(F2[1], 1)
        elT = eltype(T)
        F2 = reinterpret(reshape, elT, F2)
        if Nspecies == 1
            plan = FFTW.plan_r2r!(F2, FFTW.RODFT00; flags=FFTW.ESTIMATE)
            return plan 
        else
            plan = FFTW.plan_r2r!(F2, FFTW.RODFT00, 2, flags=FFTW.ESTIMATE)
            return plan
        end
    end
end



function find_fourier_plan_nd(::SimpleLiquid{ndims, species, T1, T2, P}, F::Vector{T}, dr::Number) where {ndims, species, T1, T2, P, T}
    M = length(F)
    plan = QDHT(Ndims/2-1, 1, M*dr, M)
    return plan
end

function inverse_radial_fourier_transform_3d(F̂, r, k)
    M = length(r)
    @assert length(k) == length(F̂) ==  M
    dk = k[2] - k[1]
    dr = r[2] - r[1]
    # @assert dk*dr ≈ π/(M+1)
    F = FFTW.r2r(F̂, FFTW.RODFT00)*dk/(4π^2)
    # same as F = dk/(2π^2)*sum(F̂ .* sin.((1:M).*(1:M)'*pi/(M+1)), dims=1)'
    return F
end

function radial_fourier_transform_3d(F, r, k)
    M = length(r)
    @assert length(k) == length(F) ==  M
    dk = k[2] - k[1]
    dr = r[2] - r[1]
    # @assert dk*dr ≈ π/(M+1)
    F̂ = FFTW.r2r(F, FFTW.RODFT00)*2π*dr
    # same as F̂ = 4π*dr*sum(F .* sin.((1:M).*(1:M)'*pi/(M+1)), dims=1)'
    return F̂ 
end


"""
    fourier!(F̂::Vector{T}, F::Vector{T}, plan::My3DPlan, dr) where T

computes the three dimensional radial fourier transform of F = r*f(r), returning F̂ = k*f̂(k), where f̂ is the fourier transform of f.
it uses the discrete sine tranform provided by the r2r function of FFTW internally.
"""
function fourier!(F̂::Vector{T}, F::Vector{T}, myplan::My3DPlan) where T
    dr = myplan.dr
    plan = myplan.plan
    @. F̂ = F*2π*dr
    if T <: Number
        plan*F̂
    elseif T<:AbstractMatrix
        F̂2 = reinterpret(reshape, eltype(T), F̂)
        plan*F̂2
    end
end

"""
    inverse_fourier!(F::AbstractVector{T}, F̂::AbstractVector{T}, plan::My3DPlan, dk) where T

computes the three dimensional radial fourier transform of F̂ = k*f̂(k), returning F = r*f(r), where f̂ is the fourier transform of f.
it uses the discrete sine tranform provided by the r2r function of FFTW internally.
"""
function inverse_fourier!(F::AbstractVector{T}, F̂::AbstractVector{T}, myplan::My3DPlan) where T
    dk = myplan.dk
    plan = myplan.plan
    @. F = F̂ * dk/(4π^2)
    if T <: Number
        plan*F
    elseif T<:AbstractMatrix
        F2 = reinterpret(reshape, eltype(T), F)
        plan*F2
    end
end

"""
    fourier!(F̂::Vector{T}, F::Vector{T}, plan::My1DPlan, dr) where T

computes the one dimensional fourier transform of F = r*f(r), returning F̂ = k*f̂(k), where f̂ is the fourier transform of f.

"""
function fourier!(F̂::Vector{T}, F::Vector{T}, myplan::My1DPlan) where T
    dr = myplan.dr
    r = myplan.r
    k = myplan.k
    for j in eachindex(F̂)
        F̂[j] = zero(eltype(F̂))
        for i in eachindex(F)
            F̂[j] += 2 * dr * F[i] / r[i] * cos(k[j]*r[i])
        end
        F̂[j] = F̂[j] * k[j]
    end
end

"""
    inverse_fourier!(F::AbstractVector{T}, F̂::AbstractVector{T}, plan::My1DPlan, dk) where T

computes the one dimensional fourier transform of F̂ = k*f̂(k), returning F = r*f(r), where f̂ is the fourier transform of f.
"""
function inverse_fourier!(F::AbstractVector{T}, F̂::AbstractVector{T}, myplan::My1DPlan) where T
    dk = myplan.dk
    r = myplan.r
    k = myplan.k
    for i in eachindex(F)
        F[i] = zero(eltype(F))
        for j in eachindex(F̂)
            F[i] += dk/π * F̂[j] / k[j] * cos(k[j]*r[i])
        end
        F[i] = F[i] * r[i]
    end
end


"""
    fourier!(F̂, F, plan, equation::OZEquation{3, T1, T2, T3}) where {T1, T2, T3}

computes the N-dimensional radial fourier transform of F = r*f(r), returning F̂ = k*f̂(k), where f̂ is the fourier transform of f.
it uses the hankel tranform of Hankel.jl internally.
"""
function fourier!(F̂::AbstractVector{T}, F::AbstractVector{T}, Q::Hankel.QDHT{p, 1, T2}) where {p, T, T2}
    # here p = dims/2 - 1
    halfdims = p+1
    k = Q.k
    r = Q.r
    @. F = F * r^(halfdims-2)
    Hankel.mul!(F̂, Q, F)
    @. F̂ = F̂ * (2π)^halfdims / k^(halfdims-2)
end

"""
inverse_fourier!(F, F̂, plan, equation::OZEquation{3, T1, T2, T3}) where {T1, T2, T3}

computes the N-dimensional radial fourier transform of F̂ = k*f̂(k), returning F = r*f(r), where f̂ is the fourier transform of f.
it uses the hankel tranform of Hankel.jl internally.
"""
function inverse_fourier!(F::AbstractVector{T}, F̂::AbstractVector{T}, Q::Hankel.QDHT{p, 1, T2}) where {p, T, T2}
    # here p = dims/2 - 1
    halfdims = p+1
    k = Q.k
    r = Q.r
    @. F̂ = F̂ * k^(halfdims-2)
    Hankel.mul!(F, Q, F̂)
    @. F ./= Q.scaleRK^2
    @. F = F * (2π)^(-halfdims) / r^(halfdims-2)
end



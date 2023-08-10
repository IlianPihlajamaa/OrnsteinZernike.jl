function get_forward_and_backward_plan(::SimpleLiquid{3, species, T1, T2, P}, F) where {species, T1, T2, P}
    return find_fourier_plans_3d(F)
end

function get_forward_and_backward_plan(::SimpleLiquid{2, species, T1, T2, P}, F) where {species, T1, T2, P}
    return find_fourier_plans_2d(F)
end

function find_fourier_plans_3d(F::Vector{T}) where T
    if T <: Number
        forward = FFTW.plan_r2r!(copy(F), FFTW.RODFT00; flags=FFTW.ESTIMATE)
        backward = FFTW.plan_r2r!(copy(F), FFTW.RODFT00; flags=FFTW.ESTIMATE)
        return forward, backward
    elseif T <: AbstractArray
        F2 = copy(F)
        Nspecies = size(F2[1], 1)
        elT = eltype(T)
        F2 = reinterpret(reshape, elT, F2)
        if Nspecies == 1
            forward = FFTW.plan_r2r!(F2, FFTW.RODFT00; flags=FFTW.ESTIMATE)
            backward = FFTW.plan_r2r!(F2, FFTW.RODFT00; flags=FFTW.ESTIMATE)
            return forward, backward
        else
            forward = FFTW.plan_r2r!(F2, FFTW.RODFT00, 2, flags=FFTW.ESTIMATE)
            backward = FFTW.plan_r2r!(F2, FFTW.RODFT00, 2, flags=FFTW.ESTIMATE)
            return forward, backward
        end
    end
end

function inverse_radial_fourier_transform_3d(F̂, r, k)
    M = length(r)
    @assert length(k) == length(F̂) ==  M
    dk = k[2] - k[1]
    dr = r[2] - r[1]
    # @assert dk*dr ≈ π/(M+1)
    F = FFTW.r2r(F̂, FFTW.RODFT00)*dk/(4π^2)
    return F
end

function radial_fourier_transform_3d(F, r, k)
    M = length(r)
    @assert length(k) == length(F) ==  M
    dk = k[2] - k[1]
    dr = r[2] - r[1]
    # @assert dk*dr ≈ π/(M+1)
    F̂ = FFTW.r2r(F, FFTW.RODFT00)*2π*dr
    return F̂ 
end


"""
    fourier!(F̂, F, plan, equation::OZEquation{3, T1, T2, T3}) where {T1, T2, T3}

computes the three dimensional radial fourier transform of F = r*f(r), returning F̂ = k*f̂(k), where f̂ is the fourier transform of f.
it uses the discrete sine tranform provided by the r2r function of FFTW internally.
"""
function fourier!(F̂::Vector{T}, F::Vector{T}, plan::FFTW.r2rFFTWPlan, dr) where T
    @. F̂ = F*2π*dr
    if T <: Number
        plan*F̂
    elseif T<:AbstractMatrix
        F̂2 = reinterpret(reshape, eltype(T), F̂)
        plan*F̂2
    end
end

"""
inverse_fourier!(F, F̂, plan, equation::OZEquation{3, T1, T2, T3}) where {T1, T2, T3}

computes the three dimensional radial fourier transform of F̂ = k*f̂(k), returning F = r*f(r), where f̂ is the fourier transform of f.
it uses the discrete sine tranform provided by the r2r function of FFTW internally.
"""
function inverse_fourier!(F::AbstractVector{T}, F̂::AbstractVector{T}, plan::FFTW.r2rFFTWPlan, dk) where T
    @. F = F̂ * dk/(4π^2)
    if T <: Number
        plan*F
    elseif T<:AbstractMatrix
        F2 = reinterpret(reshape, eltype(T), F)
        plan*F2
    end
end

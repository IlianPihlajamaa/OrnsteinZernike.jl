abstract type System end

mutable struct SimpleLiquid{dims, species, T1, T2, P} <: System
    ρ::T1
    kBT::T2
    potential::P
end

function SimpleLiquid(dims, ρ::Number, kBT, potential)
    Tρ = typeof(ρ)
    TkT = typeof(kBT)
    Tpot = typeof(potential)
    return SimpleLiquid{dims, 1, Tρ, TkT, Tpot}(ρ, kBT, potential)
end

function SimpleLiquid(dims, ρ::AbstractVector, kBT, potential)
    TkT = typeof(kBT)
    Tpot = typeof(potential)
    Ns = length(ρ)
    ρ = Diagonal(SVector{Ns}(ρ))
    Tρ = typeof(ρ)
    return SimpleLiquid{dims, Ns, Tρ, TkT, Tpot}(ρ, kBT, potential)
end

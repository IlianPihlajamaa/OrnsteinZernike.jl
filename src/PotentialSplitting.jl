abstract type DividedPotential <: Potential end


"""
    WCADivision{P<:Potential, T, UC}

fields
 - `potential::Potential` 
 - `cutoff` : (\$r_{c}\$)
 - `U_c` : \$U(r=r_c)\$

Splits the `potential` at the cutoff point using the WCA splitting rule. This means 

\$u(r) = u_{SR}(r) + U_{LR}(r),\$

where \$U_{LR}(r) = u(r)\$ for \$r > r_{c}\$, and \$U(r_{c})\$        for \$r < r_{c}\$,
and   \$USR(r) = 0\$    for \$r > r_{c}\$, and \$U(r) - U(r_{c})\$ for \$r < r_{c}\$.


"""
struct WCADivision{P<:Potential, T, UC} <: DividedPotential
    potential::P
    cutoff::T
    U_c::UC
    function WCADivision(pot::Potential, rc) 
        Uc = evaluate_potential(pot, rc)
        return new{typeof(pot), typeof(rc), typeof(Uc)}(pot, rc, Uc)
    end
end


"""
    dispersion_tail(potential::Potential, kBT, r, βu)

Return the dispersion long-range contribution associated with `potential`. The
default implementation returns zero, signalling that no tail is provided.
"""
dispersion_tail(::Potential, kBT, r::Number, βu) = zero(βu)

function dispersion_tail(potential::Potential, kBT, r::AbstractArray, βu::AbstractArray)
    return dispersion_tail.((potential,), kBT, r, βu)
end

"""
    βu, βu_LR = evaluate_long_range_potential(potential, kBT, r)

Return the total reduced potential βu together with its dispersion tail βuₗᵣ.
"""
function evaluate_long_range_potential(potential::Potential, kBT, r::Number)
    βu = evaluate_potential(potential, r) / kBT
    βu_LR = dispersion_tail(potential, kBT, r, βu)
    return βu, βu_LR
end

function evaluate_long_range_potential(potential::Potential, kBT, r::AbstractArray)
    βu_βuLR = evaluate_long_range_potential.((potential,), kBT, r)
    return first.(βu_βuLR), last.(βu_βuLR)
end

function dispersion_tail(potdiv::WCADivision, kBT, r::Number, βu)
    if r > potdiv.cutoff
        return βu
    else
        return potdiv.U_c / kBT
    end
end

function evaluate_potential(potential::WCADivision, r)
    evaluate_potential(potential.potential, r)
end

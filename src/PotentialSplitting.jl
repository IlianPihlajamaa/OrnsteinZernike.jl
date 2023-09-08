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
    βU, βULR = evalutate_long_range_potential(pot::DividedPotential, kBT, r)

return the potential and the long range part of the potential, at distance r in units of kBT. 
"""
function evalutate_long_range_potential(pot::Potential, kBT, r::Number)
    U = evaluate_potential(pot, r) / kBT
    return U, zero(U)
end
function evalutate_long_range_potential(potdiv::Potential, kBT, r::AbstractArray)
    Utuples = evalutate_long_range_potential.((potdiv, ), kBT, r)
    return first.(Utuples), last.(Utuples)
end

function evalutate_long_range_potential(potdiv::WCADivision, kBT, r::Number)
    pot = potdiv.potential
    U = evaluate_potential(pot, r) / kBT
    if r > potdiv.cutoff
        ULR = U
    else
        ULR = potdiv.U_c / kBT
    end
    return U, ULR
end

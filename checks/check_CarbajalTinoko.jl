using OrnsteinZernike
M = 10000
ρ = 0.9
kBT = 1.0
dims = 3
dr = 10.0/M
pot = HardSpheres(1.0)
system = SimpleFluid(dims, ρ, kBT, pot)
closure = CarbajalTinoko(0.3)
method = NgIteration()
sol = solve(system, closure, method)

using Roots
λ = 0.3
r, γ = 2.180762339605223, 0.0007099178871052923
e = ifelse(λ > 0, 3 + λ, 3exp(λ*r))

f = function (b)
    ω = γ + b
    @show ω
    if  abs(ω) < 0.0001
        bfunc = e*(-(ω^2/6)+ω^4/360)
    else
        y = exp(ω)
        bfunc = e*((2-ω)*y - 2 - ω)/(y - 1)
    end
    obj = b - bfunc
    @show r, γ, obj
    return obj
end
b = find_zero(f, γ)

using Plots
plot(-10:0.1:10.0, f.(-10:0.1:10.0))
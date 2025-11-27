using Test
using OrnsteinZernike, StaticArrays
#Integral equation theory for charged liquids: Model 2–2 electrolytes and the bridge function  D.‐M. Duh; A. D. J. Haymet J. Chem. Phys. 97, 7716–7729 (1992)
σ = 2.8428 * 10^-10
ϵr = 78.358
e = 1.602176634e-19
ϵ0 = 8.8541878128e-12
ϵ = ϵ0 * ϵr
kB = 1.380649e-23
B = 5377.75 * 4 * 10^-10 # m K

#  u(r,i,j) = kB*B/σ * (σ/r)^9 + Z[i]*Z[j]* (e^2)/(4*π*ϵ*r)  # J
Rm = (3^(1/4) * B^(1/8) * kB^(1/8)*(4 * π*ϵ)^(1/8)*σ)/(e^(1/4)*2^(1/4)) # minimum u

Z = [2, -2] 
T0 = 273.15
kBT = kB * (T0 + 25)
mixing_parameters = [0.05, 0.05, 0.05, 0.05, 0.05]
bjerrum_length = e^2/(4*π*ϵ*kBT) 

Es = [0.4328, 1.279, 1.713, 2.197, 2.645] # from table III


for (i,c) = enumerate([0.001, 0.02, 0.0625, 0.2, 0.5625]) # mol/L
    ρ0 = c * 6.02214076e23 * 1e3  # number density in m^-3, times 2 for 1:1 electrolyte
    ρ = [ρ0, ρ0] # number density of each species in m^-3
    κD = sqrt(1/(ϵr * ϵ0 * kBT)* e^2 * sum(ρ  .* Z.^2))

    # from here we work in units of sigma and kBT
    M = 4096
    dims = 3
    dr = 0.105 * 10^-10  / σ  
    Sones = ones(SMatrix{2,2, Float64, 4})
    pot = InversePowerLaw(kB*B/σ*Sones/kBT, Sones, 9)
    system = SimpleMixture(dims, ρ*σ^3, 1, pot)
    system = SimpleChargedMixture(system, Z, bjerrum_length / σ)
    closure = HypernettedChain()
    method = FourierIteration(M = M,dr=dr, tolerance=10^-8, mixing_parameter=mixing_parameters[i], verbose=false, max_iterations=10^6)
    sol, = solve(system, closure, method, coulombsplitting=OrnsteinZernike.NoCoulombSplitting())
    r = sol.r 
    g = sol.gr

    Ex = OrnsteinZernike.compute_excess_energy(sol, system)
    @test -Ex ≈ Es[i] atol=0.01 
end

@testset "Charged Fourier vs Ng" begin
    dims = 3
    ρ = 0.05
    kBT = 1.0
    pot = HardSpheres(1.0)

    base = SimpleFluid(dims, ρ, kBT, pot)
    system = SimpleChargedFluid(base, 1.0, 5.0)

    closure = HypernettedChain()
    M = 256
    dr = 10.0 / M
    ng_method = NgIteration(M=M, dr=dr, tolerance=1e-8, verbose=false, max_iterations=10^4)
    fourier_method = FourierIteration(M=M, dr=dr, tolerance=1e-8, verbose=false, max_iterations=10^5, mixing_parameter=0.5)

    sol_ng, = solve(system, closure, ng_method; coulombsplitting=OrnsteinZernike.NoCoulombSplitting())
    sol_fourier, = solve(system, closure, fourier_method; coulombsplitting=OrnsteinZernike.NoCoulombSplitting())

    @test maximum(abs.(sol_ng.cr .- sol_fourier.cr)) < 5e-3
    @test maximum(abs.(sol_ng.gr .- sol_fourier.gr)) < 5e-3
    @test maximum(abs.(sol_ng.ck .- sol_fourier.ck)) < 5e-3
    @test maximum(abs.(sol_ng.Sk .- sol_fourier.Sk)) < 5e-3
end

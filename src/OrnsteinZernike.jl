
"""
A generic solver package for Ornstein-Zernike equations from liquid state theory

"""
module OrnsteinZernike
using FFTW, StaticArrays, LinearAlgebra, Hankel, SpecialFunctions, Dierckx 
using Bessels: besselj
using FunctionZeros
using Roots: find_zero

export solve
export SimpleFluid, SimpleMixture, SimpleChargedFluid, SimpleChargedMixture
export OZSolution
export Exact, FourierIteration, NgIteration, DensityRamp, TemperatureRamp
export PercusYevick,  HypernettedChain, MeanSpherical, ModifiedHypernettedChain, Verlet, MartynovSarkisov
export SoftCoreMeanSpherical, RogersYoung, ZerahHansen, DuhHaymet, Lee, ChoudhuryGhosh, BallonePastoreGalliGazzillo
export VompeMartynov, CharpentierJackse, BomontBretonnet, Khanpour, ModifiedVerlet, CarbajalTinoko, ExtendedRogersYoung
export CustomPotential, PowerLaw, HardSpheres, LennardJones, SquareWell, Morse, TabulatedPotential, Yukawa, GaussianCore
export compute_compressibility, compute_excess_energy, compute_virial_pressure
export evaluate_potential, evaluate_potential_derivative, discontinuities
export WCADivision
export NoCoulombSplitting, EwaldSplitting


include("Systems.jl")
include("Potentials.jl")
include("Coulomb.jl")
include("PotentialSplitting.jl")
include("Closures.jl")
include("FourierTransforms.jl")
include("Solutions.jl")
include("Solvers.jl")
include("Thermodynamics.jl")

end # module OrnsteinZernike

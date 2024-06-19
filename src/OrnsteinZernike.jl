
"""
A generic solver package for Ornstein-Zernike equations from liquid state theory

"""
module OrnsteinZernike
using FFTW, StaticArrays, LinearAlgebra, Hankel, SpecialFunctions, Dierckx 
using Bessels: besselj
using FunctionZeros
using Roots: find_zero

export solve
export SimpleLiquid
export OZSolution
export Exact, FourierIteration, NgIteration, DensityRamp, TemperatureRamp
export PercusYevick,  HypernettedChain, MeanSpherical, ModifiedHypernettedChain, Verlet, MartynovSarkisov
export SoftCoreMeanSpherical, RogersYoung, ZerahHansen, DuhHaymet, Lee, ChoudhuryGhosh, BallonePastoreGalliGazzillo
export VompeMartynov, CharpentierJackse, BomontBretonnet, Khanpour, ModifiedVerlet, CarbajalTinoko, ExtendedRogersYoung
export CustomPotential, PowerLaw, HardSpheres, LennardJones
export compute_compressibility, compute_excess_energy, compute_virial_pressure
export WCADivision

include("Systems.jl")
include("Potentials.jl")
include("PotentialSplitting.jl")
include("Closures.jl")
include("FourierTransforms.jl")
include("Solutions.jl")
include("Solvers.jl")
include("Thermodynamics.jl")

end # module OrnsteinZernike

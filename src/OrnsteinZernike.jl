
"""
A generic solver package for Ornstein-Zernike equations from liquid state theory

"""
module OrnsteinZernike
using FFTW, StaticArrays, LinearAlgebra, Hankel

export solve
export SimpleLiquid
export OZSolution
export Exact, FourierIteration, NgIteration, DensityRamp
export PercusYevick,  HypernettedChain, MeanSphericalApproximation
export CustomPotential, PowerLaw, HardSpheres, LennardJones
export compute_compressibility, compute_excess_energy, compute_virial_pressure

include("Systems.jl")
include("Potentials.jl")
include("Closures.jl")
include("FourierTransforms.jl")
include("Solutions.jl")
include("Solvers.jl")
include("Thermodynamics.jl")

end # module OrnsteinZernike

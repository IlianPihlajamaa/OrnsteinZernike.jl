using Test, OrnsteinZernike, StaticArrays
using FFTW
using Hankel
using Random
Random.seed!(523)

for target in ["Fourier", "FourierIteration_HS", "NgIteration_HS", "thermodynamics", "DensityRamp"]
    @testset "$target" begin
        include("test_$target.jl")
    end
end
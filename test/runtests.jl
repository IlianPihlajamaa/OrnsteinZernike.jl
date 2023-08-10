using Test, OrnsteinZernike, StaticArrays
using FFTW

for target in ["Fourier", "FourierIteration_HS", "NgIteration_HS", "thermodynamics", "DensityRamp"]
    @testset "$target" begin
        include("test_$target.jl")
    end
end
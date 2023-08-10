using Test, OrnsteinZernike
using FFTW

for target in ["Fourier", "exact_PY", "FourierIteration_HS", "NGIteration_HS", "thermodynamics", "DensityRamp"]
    @testset "$target" begin
        include("test_$target.jl")
    end
end
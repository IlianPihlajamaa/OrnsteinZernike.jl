using Test, OrnsteinZernike, StaticArrays
using FFTW
using Hankel
using Random
using InteractiveUtils: subtypes
import Roots
Random.seed!(523)

function main_test(target)
    @testset "$target" begin
        include("test_$target.jl")
    end
end


for target in ["Fourier", "FourierIteration_HS", "NgIteration_HS", "thermodynamics", "DensityRamp", "dims", "RY", "traits", "closures", "potentials", "charged"]
    main_test(target)
end


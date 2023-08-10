

M = 2^8

F = rand(M)
dr = 0.01
dk = π/(M+1)/dr
r = collect((1:M)*dr)
k = collect((1:M)*dk)

plan1 = FFTW.plan_r2r!(copy(F), FFTW.RODFT00; flags=FFTW.ESTIMATE);
plan2 = FFTW.plan_r2r!(copy(F), FFTW.RODFT00; flags=FFTW.ESTIMATE)

Fhat = OrnsteinZernike.radial_fourier_transform_3d(F, r, k);
Fnew = OrnsteinZernike.inverse_radial_fourier_transform_3d(Fhat, r, k);

@test maximum(abs.(F-Fnew)) .< 1e-13

Fnew2 = copy(F) 
Fhat2 = copy(Fhat) 

OrnsteinZernike.inverse_fourier!(Fnew2, Fhat, plan1, dk);
@test maximum(abs.(F-Fnew2))  .< 1e-13

OrnsteinZernike.fourier!(Fhat2, F, plan2, dr);

@test maximum(abs.(Fhat-Fhat2))  .< 1e-13


for Nspecies = [1, 2, 5] 
    FF = rand(StaticArrays.SMatrix{Nspecies, Nspecies, Float64, Nspecies*Nspecies}, M)
    forwardplan, backwardplan =  OrnsteinZernike.find_fourier_plans_3d(FF)
    FFhat = copy(FF)
    OrnsteinZernike.fourier!(FFhat, FF, forwardplan, dr);
    F2 = copy(FF) 
    OrnsteinZernike.inverse_fourier!(F2, FFhat, backwardplan, dk);
    @test maximum(abs.(maximum.(FF-F2)))  .< 1e-13
    Fhat11 = OrnsteinZernike.radial_fourier_transform_3d(getindex.(FF, 1, 1), r, k);
    @test maximum(abs.(getindex.(FFhat, 1, 1)-Fhat11))  .< 1e-13
end



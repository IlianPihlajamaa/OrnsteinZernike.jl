using FFTW
M = 10000
F = rand(M)
dr = 0.01
dk = π/(M*dr)
r = collect(((1:M) .- 0.5)*dr)
k = collect(((1:M) .- 0.5)*dk)

plan1 = FFTW.plan_r2r!(copy(F), FFTW.RODFT11; flags=FFTW.ESTIMATE);
myplan1 = OrnsteinZernike.My3DPlan(plan1, r, k, dr,dk,M)

Fhat = OrnsteinZernike.radial_fourier_transform_3d(F, r, k);
Fnew = OrnsteinZernike.inverse_radial_fourier_transform_3d(Fhat, r, k);
maximum(abs.(F-Fnew)) .< 1e-13

Fnew2 = copy(F) 
Fhat2 = copy(Fhat) 

OrnsteinZernike.inverse_fourier!(Fnew2, Fhat, myplan1);
maximum(abs.(F-Fnew2))  .< 1e-13
OrnsteinZernike.fourier!(Fhat2, F, myplan1);

maximum(abs.(Fhat-Fhat2))  .< 1e-13
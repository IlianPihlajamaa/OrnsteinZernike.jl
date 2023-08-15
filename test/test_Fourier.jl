
for M = [2^8]
    F = rand(M)
    dr = 0.01
    dk = π/(M+1)/dr
    r = collect((1:M)*dr)
    k = collect((1:M)*dk)

    plan1 = FFTW.plan_r2r!(copy(F), FFTW.RODFT00; flags=FFTW.ESTIMATE);
    myplan1 = OrnsteinZernike.My3DPlan(plan1, r, k, dr,dk,M)

    Fhat = OrnsteinZernike.radial_fourier_transform_3d(F, r, k);
    Fnew = OrnsteinZernike.inverse_radial_fourier_transform_3d(Fhat, r, k);

    @test maximum(abs.(F-Fnew)) .< 1e-13

    Fnew2 = copy(F) 
    Fhat2 = copy(Fhat) 

    OrnsteinZernike.inverse_fourier!(Fnew2, Fhat, myplan1);
    @test maximum(abs.(F-Fnew2))  .< 1e-13

    OrnsteinZernike.fourier!(Fhat2, F, myplan1);

    @test maximum(abs.(Fhat-Fhat2))  .< 1e-13
end

for Nspecies = [1, 2, 5] 
    M = 2^8
    F = rand(M)
    dr = 0.01
    dk = π/(M+1)/dr
    r = collect((1:M)*dr)
    k = collect((1:M)*dk)
    FF = rand(StaticArrays.SMatrix{Nspecies, Nspecies, Float64, Nspecies*Nspecies}, M)
    plan =  OrnsteinZernike.find_fourier_plan_3d(FF)
    myplan = OrnsteinZernike.My3DPlan(plan,r,k, dr,dk,M)
    FFhat = copy(FF)
    OrnsteinZernike.fourier!(FFhat, FF, myplan);
    F2 = copy(FF) 
    OrnsteinZernike.inverse_fourier!(F2, FFhat, myplan);
    @test maximum(abs.(maximum.(FF-F2)))  .< 1e-13
    Fhat11 = OrnsteinZernike.radial_fourier_transform_3d(getindex.(FF, 1, 1), r, k);
    @test maximum(abs.(getindex.(FFhat, 1, 1)-Fhat11))  .< 1e-13
end


## Test ND Fourier
for dims = [3]
    R = 10.0
    N = 100
    Q = QDHT(dims/2-1, 1, R, N)
    r = Q.r
    k = Q.k
    F = @. exp(-r^2)
    fk = OrnsteinZernike.radial_fourier_transform_3d(F.*r, r, k)./k
    fk2 = (2π)^(dims/2)*(Q*(F.*r.^(dims/2-1)) ./ k.^(dims/2-1))
    fk3 = @. exp(-k^2 / 4) * π^(3/2) # exact
    fk4 = similar(fk3)
    OrnsteinZernike.fourier!(fk4, F.*r, Q)
    fk4 ./= k
    @test fk ≈ fk3
    @test fk2 ≈ fk3
    @test fk4 ≈ fk3

    F2 = similar(F)
    OrnsteinZernike.inverse_fourier!(F2, fk.*k, Q) 
    F2 ./= r
    F3 = (2π)^(-dims/2)*(Q\(fk.*k.^(dims/2-1)) ./ r.^(dims/2-1))
    F4 = OrnsteinZernike.inverse_radial_fourier_transform_3d(fk.*k, r, k)./r
    @test F ≈ F2
    @test F ≈ F3
    @test F ≈ F4
end
for Nspecies = [1, 2, 3]
    for dims = [2,3,4,7] 
        R = 3.0
        N = 100
        Q = QDHT(dims/2-1, 1, R, N)
        r = Q.r
        dr = R/N
        k = Q.k
        α = rand(StaticArrays.SMatrix{Nspecies, Nspecies, Float64, Nspecies*Nspecies})
        F = [exp.(-0α.+-ri^2) for ri in r]

        F11 = getindex.(F, 1, 1)
        F11k = similar(F11)
        OrnsteinZernike.fourier!(F11k, F11.*r, Q)
        F11k ./= k
        Fhat = similar(F)
        OrnsteinZernike.fourier!(Fhat, F.*r, Q)
        Fhat ./= k

        F11k2 = getindex.(Fhat, 1, 1)
        @test F11k2 ≈ F11k

        F2 = similar(F)
        OrnsteinZernike.inverse_fourier!(F2, Fhat.*k, Q)
        F2 ./= r

        F3 = similar(F11)
        OrnsteinZernike.inverse_fourier!(F3, getindex.(Fhat.*k, 1, 1), Q)
        F3 ./= r

        @test F2 ≈ F
        @test F3 ≈ getindex.(F, 1, 1)
    end
end

## test with 1d tranform

for M = [100, 200, 300]
    R = 10.0
    dr = R/M
    dims = 1

    dk =  π/((M+0.5)*dr)
    r = [i*dr for i = 0.5:(M-0.5)]
    k = [j*dk for j = 0.5:(M-0.5)]
    plan = OrnsteinZernike.My1DPlan(r, k, dr, dk, M)
    F = @. exp(-r^2)

    fk3 = @. exp(-k^2 / 4) * π^(dims/2) # exact
    fk4 = similar(fk3)
    OrnsteinZernike.fourier!(fk4, F.*r, plan)
    fk4 ./= k
    @test fk4 ≈ fk3

    F2 = similar(F)
    OrnsteinZernike.inverse_fourier!(F2, fk4.*k, plan) 
    F2 ./= r
    @test F ≈ F2
end
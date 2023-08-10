

function main()
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
    
    println("Testing out-of-place")
    
    @test maximum(abs.(F-Fnew)) .< 1e-15
    
    Fnew2 = copy(F) 
    Fhat2 = copy(Fhat) 
    
    OrnsteinZernike.inverse_fourier!(Fnew2, Fhat, plan1, dk);
    println("Testing in-place inverse")
    @test maximum(abs.(F-Fnew2))  .< 1e-15

    @time OrnsteinZernike.fourier!(Fhat2, F, plan2, dr);
    
    println("Testing in-place forward")
    @test maximum(abs.(Fhat-Fhat2))  .< 1e-15

    println("testing multicomponent transform")

    for Nspecies = [1, 2, 5] 
        F = rand(StaticArrays.SMatrix{Nspecies, Nspecies, Float64, Nspecies*Nspecies}, M)
        forwardplan, backwardplan =  OrnsteinZernike.find_fourier_plans_3d(F)
        Fhat = copy(F)
        OrnsteinZernike.fourier!(Fhat, F, forwardplan, dr);
        F2 = copy(F) 
        OrnsteinZernike.inverse_fourier!(F2, Fhat, backwardplan, dk);
        @test maximum(abs.(maximum.(F-F2)))  .< 1e-15
        Fhat11 = OrnsteinZernike.radial_fourier_transform_3d(getindex.(F, 1, 1), r, k);
        @test maximum(abs.(getindex.(Fhat, 1, 1)-Fhat11))  .< 1e-15
    end
end
main()


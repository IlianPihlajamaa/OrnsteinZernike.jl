using Test
using OrnsteinZernike

const OZ = OrnsteinZernike

@testset "Closure renormalization trait" begin
    @test !OZ.uses_renormalized_gamma(OZ.PercusYevick())
    @test OZ.uses_renormalized_gamma(OZ.SMSA())
    @test OZ.uses_renormalized_gamma(OZ.ZerahHansen())
    @test OZ.uses_renormalized_gamma(OZ.BomontBretonnet())
end

struct RecordingClosure <: OZ.Closure
    γ_ref::Base.RefValue{Float64}
    renorm::Bool
end

RecordingClosure(; renorm=false) = RecordingClosure(Base.RefValue(0.0), renorm)

OZ.uses_renormalized_gamma(c::RecordingClosure) = c.renorm

function OZ.bridge_function(c::RecordingClosure, _, _, γ)
    c.γ_ref[] = γ
    return γ
end

function expected_cmulr_uncharged(γ_SR, mayer_f, γ_bridge)
    return -1 - γ_SR + (mayer_f + 1) * exp(γ_SR) * exp(γ_bridge)
end

function expected_cmulr_charged(γ_SR, mayer_f, γ_bridge, βu_LR_coul, q)
    return -1 - γ_SR - q + (mayer_f + 1) * exp(βu_LR_coul + γ_SR + q) * exp(γ_bridge)
end

@testset "closure_apply! contexts" begin
    r = [1.0]
    Γ = [0.5]
    βu = [0.3]
    βu_disp = [0.1]
    mayer_f = [exp(-βu[1]) - 1]
    dest = similar(βu)

    for renorm in (false, true)
        closure = RecordingClosure(renorm=renorm)
        ctx = OZ.UnchargedClosureEvalContext(r, mayer_f, Γ, βu, βu_disp)
        OZ.closure_apply!(dest, closure, ctx)
        γ_SR = Γ[1] / r[1]
        γ_expected = renorm ? γ_SR - βu_disp[1] : γ_SR
        @test isapprox(closure.γ_ref[], γ_expected; atol=1e-12)
        expected = expected_cmulr_uncharged(γ_SR, mayer_f[1], γ_expected) * r[1]
        @test isapprox(dest[1], expected; atol=1e-12)
    end

    Γc = [0.6]
    βu_c = [0.4]
    βu_disp_c = [0.15]
    βu_LR_coul = [0.05]
    q = [0.07]
    mayer_fc = [exp(-βu_c[1]) - 1]
    destc = similar(βu_c)

    for renorm in (false, true)
        closure = RecordingClosure(renorm=renorm)
        ctx = OZ.ChargedClosureEvalContext(r, mayer_fc, Γc, βu_c, βu_disp_c, βu_LR_coul, q)
        OZ.closure_apply!(destc, closure, ctx)
        γ_SR = Γc[1] / r[1]
        γ_expected = renorm ? γ_SR - βu_disp_c[1] : γ_SR
        @test isapprox(closure.γ_ref[], γ_expected; atol=1e-12)
        expected = expected_cmulr_charged(γ_SR, mayer_fc[1], γ_expected, βu_LR_coul[1], q[1]) * r[1]
        @test isapprox(destc[1], expected; atol=1e-12)
    end
end

@testset "Dispersion tail helpers" begin
    pot = OZ.LennardJones(1.0, 1.0)
    @test_throws ErrorException OZ.evaluate_long_range_potential(pot, 1.0, 1.5)

    rc = 2^(1/6)
    wca = OZ.WCADivision(pot, rc)
    βu_core, βu_tail_core = OZ.evaluate_long_range_potential(wca, 1.0, rc/2)
    @test isapprox(βu_tail_core, OZ.evaluate_potential(pot, rc); atol=1e-12)
    βu_far, βu_tail_far = OZ.evaluate_long_range_potential(wca, 1.0, rc*1.5)
    @test isapprox(βu_tail_far, βu_far; atol=1e-12)
end

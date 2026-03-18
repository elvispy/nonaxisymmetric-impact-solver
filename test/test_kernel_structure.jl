@testset "Corrected-kernel structure" begin
    op = build_toeplitz_dtn(9, 9, 0.1)
    K = op.kernel
    @test K ≈ reverse(K, dims=1) atol=1e-12 rtol=1e-12
    @test K ≈ reverse(K, dims=2) atol=1e-12 rtol=1e-12
    @test K[9, 9] > 0.0
    @test any(K .< 0.0)
    @test abs(sum(K)) < 1e-12
end

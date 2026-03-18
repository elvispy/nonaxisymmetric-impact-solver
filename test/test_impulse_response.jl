@testset "Impulse response is nonperiodic" begin
    op = build_toeplitz_dtn(9, 9, 0.1)
    x = zeros(9, 9)
    x[1, 1] = 1.0
    y = apply_direct(op, x)

    @test y[1, 1] ≈ op.kernel[9, 9] atol=1e-12 rtol=1e-12
    @test y[9, 9] ≈ op.kernel[1, 1] atol=1e-12 rtol=1e-12
    @test abs(y[9, 1]) < abs(op.kernel[9, 9])
end

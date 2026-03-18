@testset "Corrected-kernel FFT apply matches direct apply" begin
    op = build_toeplitz_dtn(9, 9, 0.2)
    x = reshape(collect(1.0:81.0), 9, 9)
    y_direct = apply_direct(op, x)
    y_fft = apply(op, x)
    @test y_fft ≈ y_direct atol=1e-10 rtol=1e-10
end

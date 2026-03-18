@testset "DtN API" begin
    corrected = build_toeplitz_dtn(8, 8, 0.1)
    periodic = build_periodic_dtn(8, 8, 0.1)

    @test corrected isa CorrectedKernelDtN2D
    @test periodic isa PeriodicDtN2D
    @test size(corrected.kernel) == (15, 15)
    @test size(periodic.kernel) == (8, 8)
end

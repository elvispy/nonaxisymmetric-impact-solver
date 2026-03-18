@testset "Periodic constants map to zero" begin
    op = build_periodic_dtn(16, 16, 0.1)
    x = ones(16, 16)
    y = apply(op, x)
    @test maximum(abs, y) < 1e-12
end

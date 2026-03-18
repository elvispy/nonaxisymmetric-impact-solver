@testset "Periodic Fourier mode check" begin
    n = 63
    dx = 2π / n
    op = build_periodic_dtn(n, n, dx)
    xs = reshape((0:n-1) .* dx, n, 1)
    ys = reshape((0:n-1) .* dx, 1, n)
    kx = 2
    ky = 1
    f = @. cos(kx * xs + ky * ys)
    target = sqrt(kx^2 + ky^2) .* f
    got = apply(op, f)
    relerr = norm(got - target) / norm(target)
    @test relerr < 1e-10
end

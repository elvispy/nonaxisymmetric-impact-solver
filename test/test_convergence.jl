@testset "Corrected-kernel center self-convergence" begin
    box_length = 20.0
    a = 1.0
    sizes = [33, 65, 129]
    nref = 257

    dx_ref = box_length / nref
    xs_ref, ys_ref = grid_coordinates(nref, nref, dx_ref)
    r2_ref = @. xs_ref^2 + ys_ref^2
    f_ref = @. a / (2 * π * (r2_ref + a^2)^(3 / 2))
    op_ref = build_toeplitz_dtn(nref, nref, dx_ref)
    y_ref = apply(op_ref, f_ref)
    center_ref = y_ref[div(nref + 1, 2), div(nref + 1, 2)]

    errors = Float64[]
    for n in sizes
        dx = box_length / n
        xs, ys = grid_coordinates(n, n, dx)
        r2 = @. xs^2 + ys^2
        f = @. a / (2 * π * (r2 + a^2)^(3 / 2))
        op = build_toeplitz_dtn(n, n, dx)
        y = apply(op, f)
        center = y[div(n + 1, 2), div(n + 1, 2)]
        push!(errors, abs(center - center_ref))
    end

    orders = [
        log(errors[i] / errors[i + 1]) / log(2.0)
        for i in 1:(length(errors) - 1)
    ]

    @test errors[end] < errors[1]
    @test minimum(orders) > 1.8

    exact_center = 1 / π
    @test abs(center_ref - exact_center) < 1e-2
end

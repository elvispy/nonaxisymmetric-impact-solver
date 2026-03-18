@testset "Solver domain moments" begin
    nx = 5
    ny = 5
    contact = falses(nx, ny)
    contact[2:4, 2:4] .= true
    density = ones(nx, ny)

    cfg = SolverConfig(
        nx=nx,
        ny=ny,
        dx=1.0,
        dt=0.1,
        steps=1,
        Fr=4.0,
        We=2.0,
        Re=8.0,
        contact_mask=contact,
        mass_density=density,
    )
    dom = build_solver_domain(cfg)

    @test isapprox(dom.mass_matrix[1, 2], 0.0; atol=1e-12)
    @test isapprox(dom.mass_matrix[1, 3], 0.0; atol=1e-12)
    @test isapprox(dom.mass_matrix[2, 3], 0.0; atol=1e-12)
    @test rank(dom.plane_basis) == 3
    @test maximum(abs, dom.surface_tension_map * ones(nx * ny)) < 1e-12
end

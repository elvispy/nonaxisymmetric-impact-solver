@testset "Solver assembly is square" begin
    nx = 4
    ny = 4
    contact = falses(nx, ny)
    contact[2:3, 2:3] .= true
    density = ones(nx, ny)

    cfg = SolverConfig(
        nx=nx,
        ny=ny,
        dx=0.5,
        dt=0.1,
        steps=1,
        Fr=3.0,
        We=2.0,
        Re=10.0,
        contact_mask=contact,
        mass_density=density,
    )
    system = build_solver(cfg)
    state = zero_state(system)
    A, b = assemble_step_system(system, state)

    N = nx * ny
    @test size(A, 1) == size(A, 2) == 2 * N + 6
    @test length(b) == 2 * N + 6
end

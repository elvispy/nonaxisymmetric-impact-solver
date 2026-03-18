@testset "Solver API" begin
    nx = 5
    ny = 5
    contact = falses(nx, ny)
    contact[2:4, 2:4] .= true
    density = ones(nx, ny)

    cfg = SolverConfig(
        nx=nx,
        ny=ny,
        dx=0.5,
        dt=0.1,
        steps=2,
        Fr=2.0,
        We=3.0,
        Re=10.0,
        contact_mask=contact,
        mass_density=density,
    )

    system = build_solver(cfg)
    state = zero_state(system)

    @test system isa SolverSystem
    @test state isa SolverState
    @test size(system.domain.laplacian) == (nx * ny, nx * ny)
    @test length(system.domain.contact_indices) == count(contact)
end

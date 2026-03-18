@testset "Matrix-free operator matches assembled system" begin
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
        linear_solver=:direct,
        contact_mask=contact,
        mass_density=density,
    )
    system = build_solver(cfg)
    state = zero_state(system)
    A, _ = assemble_step_system(system, state)
    x = collect(range(-1.0, 1.0; length=size(A, 2)))

    @test apply_step_operator(system, x) ≈ A * x atol=1e-10 rtol=1e-10
end

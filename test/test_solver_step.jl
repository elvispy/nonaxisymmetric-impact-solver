@testset "Solver step behavior" begin
    nx = 5
    ny = 5
    contact = falses(nx, ny)
    contact[2:4, 2:4] .= true
    density = ones(nx, ny)

    cfg = SolverConfig(
        nx=nx,
        ny=ny,
        dx=0.5,
        dt=0.05,
        steps=1,
        Fr=Inf,
        We=2.0,
        Re=20.0,
        contact_mask=contact,
        mass_density=density,
    )
    system = build_solver(cfg)
    state0 = zero_state(system)
    state1 = advance_one_step(state0, system)

    @test maximum(abs, state1.eta) < 1e-12
    @test maximum(abs, state1.phi) < 1e-12
    @test maximum(abs, state1.pressure) < 1e-12
    @test maximum(abs, state1.q) < 1e-12
    @test maximum(abs, state1.v) < 1e-12

    cfg2 = SolverConfig(
        nx=nx,
        ny=ny,
        dx=0.5,
        dt=0.05,
        steps=1,
        Fr=2.0,
        We=2.0,
        Re=20.0,
        body_damping=0.1,
        contact_mask=contact,
        mass_density=density,
        external_load=[1.0, 0.0, 0.0],
    )
    system2 = build_solver(cfg2)
    next_state = advance_one_step(zero_state(system2), system2)

    contact_indices = system2.domain.contact_indices
    plane_values = system2.domain.plane_basis * next_state.q
    @test vec(next_state.eta)[contact_indices] ≈ plane_values atol=1e-10 rtol=1e-10
    @test maximum(abs, next_state.pressure[.!contact]) < 1e-12

    cfg3 = SolverConfig(
        nx=nx,
        ny=ny,
        dx=0.5,
        dt=0.05,
        steps=1,
        Fr=2.0,
        We=2.0,
        Re=20.0,
        body_damping=0.1,
        contact_mask=contact,
        mass_density=density,
        external_load=[1.0, 0.0, 0.0],
        linear_solver=:direct,
    )
    direct_state = advance_one_step(zero_state(build_solver(cfg3)), build_solver(cfg3))
    @test next_state.eta ≈ direct_state.eta atol=1e-8 rtol=1e-8
    @test next_state.phi ≈ direct_state.phi atol=1e-8 rtol=1e-8
    @test next_state.pressure ≈ direct_state.pressure atol=1e-8 rtol=1e-8
    @test next_state.q ≈ direct_state.q atol=1e-8 rtol=1e-8
    @test next_state.v ≈ direct_state.v atol=1e-8 rtol=1e-8
end

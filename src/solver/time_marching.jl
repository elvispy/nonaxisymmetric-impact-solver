"""
    advance_one_step(state, system)

Advance the coupled solver by one implicit Euler step using either the dense reference
solve or the matrix-free GMRES path.
"""
function advance_one_step(state::SolverState, system::SolverSystem)
    sol = if system.config.linear_solver == :direct
        A, b = assemble_step_system(system, state)
        A \ b
    elseif system.config.linear_solver == :gmres
        b = build_step_rhs(system, state)
        gmres_solve(system, b)
    else
        throw(ArgumentError("unsupported linear solver: $(system.config.linear_solver)"))
    end

    cfg = system.config
    dom = system.domain
    layout = step_layout(system)

    eta_free = sol[layout.eta_rng]
    phi_next = reshape(sol[layout.phi_rng], cfg.nx, cfg.ny)
    pressure_contact = sol[layout.p_rng]
    v_next = sol[layout.v_rng]
    q_next = sol[layout.q_rng]

    eta_vec = dom.eta_free_map * eta_free + dom.q_map * q_next
    eta_next = reshape(eta_vec, cfg.nx, cfg.ny)

    pressure_vec = dom.pressure_map * pressure_contact
    pressure_next = reshape(pressure_vec, cfg.nx, cfg.ny)

    return SolverState(
        state.time + cfg.dt,
        eta_next,
        phi_next,
        pressure_next,
        q_next,
        v_next,
    )
end

"""
    solve_motion(config; initial_state=nothing)

Advance the prescribed-contact impact model for `config.steps` timesteps and return
the full state history.
"""
function solve_motion(config::SolverConfig; initial_state::Union{Nothing, SolverState}=nothing)
    system = build_solver(config)
    state = isnothing(initial_state) ? zero_state(system) : initial_state
    history = Vector{SolverState}(undef, config.steps + 1)
    history[1] = state

    for step in 1:config.steps
        state = advance_one_step(state, system)
        history[step + 1] = state
    end

    return history
end

"""
    assemble_step_system(system, state)

Assemble the dense implicit Euler step system. This is the correctness-reference path
used to validate the matrix-free operator on small problems.
"""
function assemble_step_system(system::SolverSystem, state::SolverState)
    cfg = system.config
    dom = system.domain

    N = dom.N
    F = length(dom.free_indices)
    M = length(dom.contact_indices)

    eta_prev = vec(state.eta)
    phi_prev = vec(state.phi)

    total = F + N + M + 3 + 3
    A = zeros(Float64, total, total)
    b = zeros(Float64, total)

    eta_rng = 1:F
    phi_rng = (last(eta_rng) + 1):(last(eta_rng) + N)
    p_rng = (last(phi_rng) + 1):(last(phi_rng) + M)
    v_rng = (last(p_rng) + 1):(last(p_rng) + 3)
    q_rng = (last(v_rng) + 1):(last(v_rng) + 3)

    Tη = dom.eta_free_map
    Tq = dom.q_map
    Sp = dom.pressure_map
    L = dom.laplacian
    D = build_dense_dtn_matrix(
        cfg.nx,
        cfg.ny,
        cfg.dx;
        near_radius=cfg.near_radius,
        quadrature_order=cfg.quadrature_order,
    )

    η_all_coeff = Tη
    q_all_coeff = Tq
    IN = Matrix{Float64}(I, N, N)
    I3 = Matrix{Float64}(I, 3, 3)

    kin_rows = 1:N
    A[kin_rows, eta_rng] .= Tη - cfg.dt * (2.0 / cfg.Re) .* (L * η_all_coeff)
    A[kin_rows, phi_rng] .= -cfg.dt .* D
    A[kin_rows, q_rng] .= Tq - cfg.dt * (2.0 / cfg.Re) .* (L * q_all_coeff)
    b[kin_rows] .= eta_prev

    dyn_rows = (last(kin_rows) + 1):(last(kin_rows) + N)
    A[dyn_rows, eta_rng] .= (cfg.dt / cfg.Fr) .* η_all_coeff - (cfg.dt / cfg.We) .* (L * η_all_coeff)
    A[dyn_rows, phi_rng] .= IN - cfg.dt * (2.0 / cfg.Re) .* L
    A[dyn_rows, p_rng] .= cfg.dt .* Sp
    A[dyn_rows, q_rng] .= (cfg.dt / cfg.Fr) .* q_all_coeff - (cfg.dt / cfg.We) .* (L * q_all_coeff)
    b[dyn_rows] .= phi_prev

    body_rows = (last(dyn_rows) + 1):(last(dyn_rows) + 3)
    damping = cfg.dt * cfg.body_damping .* I3
    if cfg.include_surface_tension_boundary
        A[body_rows, eta_rng] .-= (cfg.dt / cfg.We) .* (dom.surface_tension_map * η_all_coeff)
        A[body_rows, q_rng] .-= (cfg.dt / cfg.We) .* (dom.surface_tension_map * q_all_coeff)
    end
    A[body_rows, p_rng] .= -cfg.dt .* dom.force_map
    A[body_rows, v_rng] .= dom.mass_matrix + damping
    b[body_rows] .= dom.mass_matrix * state.v .+ cfg.dt .* cfg.external_load .+ cfg.dt .* dom.gravity_load

    q_rows = (last(body_rows) + 1):(last(body_rows) + 3)
    A[q_rows, v_rng] .= -cfg.dt .* I3
    A[q_rows, q_rng] .= I3
    b[q_rows] .= state.q

    return A, b
end

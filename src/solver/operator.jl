"""
    step_layout(system)

Return the unknown and row ranges used by the coupled timestep operator.
"""
function step_layout(system::SolverSystem)
    dom = system.domain
    F = length(dom.free_indices)
    N = dom.N
    M = length(dom.contact_indices)

    eta_rng = 1:F
    phi_rng = (last(eta_rng) + 1):(last(eta_rng) + N)
    p_rng = (last(phi_rng) + 1):(last(phi_rng) + M)
    v_rng = (last(p_rng) + 1):(last(p_rng) + 3)
    q_rng = (last(v_rng) + 1):(last(v_rng) + 3)
    total = last(q_rng)
    kin_rows = 1:N
    dyn_rows = (last(kin_rows) + 1):(last(kin_rows) + N)
    body_rows = (last(dyn_rows) + 1):(last(dyn_rows) + 3)
    q_rows = (last(body_rows) + 1):(last(body_rows) + 3)

    return (; eta_rng, phi_rng, p_rng, v_rng, q_rng, kin_rows, dyn_rows, body_rows, q_rows, total)
end

"""
    build_step_rhs(system, state)

Build the right-hand side for one implicit Euler step.
"""
function build_step_rhs(system::SolverSystem, state::SolverState)
    cfg = system.config
    dom = system.domain
    layout = step_layout(system)

    b = zeros(Float64, layout.total)
    b[layout.kin_rows] .= vec(state.eta)
    b[layout.dyn_rows] .= vec(state.phi)
    b[layout.body_rows] .= dom.mass_matrix * state.v .+ cfg.dt .* cfg.external_load .+ cfg.dt .* dom.gravity_load
    b[layout.q_rows] .= state.q
    return b
end

"""
    apply_step_operator(system, x)

Apply the full coupled timestep operator without assembling the dense matrix.
"""
function apply_step_operator(system::SolverSystem, x::AbstractVector{<:Real})
    y = zeros(Float64, length(x))
    apply_step_operator!(y, system, x)
    return y
end

"""
    apply_step_operator!(y, system, x)

In-place matrix-free application of the coupled timestep operator.
"""
function apply_step_operator!(y::AbstractVector{Float64}, system::SolverSystem, x::AbstractVector{<:Real})
    cfg = system.config
    dom = system.domain
    layout = step_layout(system)
    length(x) == layout.total || throw(DimensionMismatch("input vector has wrong length"))
    length(y) == layout.total || throw(DimensionMismatch("output vector has wrong length"))

    fill!(y, 0.0)

    η_free = view(x, layout.eta_rng)
    φ = reshape(view(x, layout.phi_rng), cfg.nx, cfg.ny)
    p = view(x, layout.p_rng)
    v = view(x, layout.v_rng)
    q = view(x, layout.q_rng)

    # Reconstruct the full surface field from free-surface unknowns and rigid-plane contact values.
    η_all = dom.eta_free_map * η_free + dom.q_map * q
    Lη = dom.laplacian * η_all
    Dφ = vec(apply(dom.dtn_operator, φ))

    y[layout.kin_rows] .= η_all .- cfg.dt .* (2.0 / cfg.Re) .* Lη .- cfg.dt .* Dφ
    y[layout.dyn_rows] .=
        (cfg.dt / cfg.Fr) .* η_all .- (cfg.dt / cfg.We) .* Lη .+
        vec(φ) .- cfg.dt .* (2.0 / cfg.Re) .* (dom.laplacian * vec(φ)) .+
        cfg.dt .* (dom.pressure_map * p)
    y[layout.body_rows] .=
        (dom.mass_matrix + cfg.dt * cfg.body_damping .* Matrix{Float64}(I, 3, 3)) * v .-
        cfg.dt .* (dom.force_map * p)
    if cfg.include_surface_tension_boundary
        y[layout.body_rows] .-= (cfg.dt / cfg.We) .* (dom.surface_tension_map * η_all)
    end
    y[layout.q_rows] .= q .- cfg.dt .* v

    return y
end

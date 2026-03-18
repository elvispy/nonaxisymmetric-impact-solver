#!/usr/bin/env julia

ENV["GKSwstype"] = "100"

using Dates
using VerticalDynamics

"""Build a small centered rectangular contact mask."""
function rectangular_contact_mask(nx::Int, ny::Int; halfwidth_x::Int=3, halfwidth_y::Int=3)
    mask = falses(nx, ny)
    ic = fld(nx, 2) + 1
    jc = fld(ny, 2) + 1
    mask[(ic - halfwidth_x):(ic + halfwidth_x), (jc - halfwidth_y):(jc + halfwidth_y)] .= true
    return mask
end

"""Run one representative case and write an MP4 heatmap of the free-surface elevation."""
function main()
    nx = 25
    ny = 25
    dx = 0.08
    dt = 0.02
    steps = 40

    contact_mask = rectangular_contact_mask(nx, ny; halfwidth_x=2, halfwidth_y=2)
    mass_density = ones(nx, ny)

    cfg = SolverConfig(
        nx=nx,
        ny=ny,
        dx=dx,
        dt=dt,
        steps=steps,
        Fr=1.0,
        We=1.0,
        Re=25.0,
        body_damping=0.2,
        contact_mask=contact_mask,
        mass_density=mass_density,
        external_load=zeros(3),
        near_radius=2,
        quadrature_order=8,
        linear_solver=:gmres,
        krylov_tol=1e-9,
        krylov_maxiter=200,
        include_surface_tension_boundary=true,
    )

    system = build_solver(cfg)
    state0 = zero_state(system)
    state0 = SolverState(state0.time, state0.eta, state0.phi, state0.pressure, [0.06, 0.0, 0.0], zeros(3))

    history = solve_motion(cfg; initial_state=state0)
    x, y = grid_coordinates(nx, ny, dx)
    x_full = repeat(x, 1, ny)
    y_full = repeat(y, nx, 1)

    tag = Dates.format(now(), "yyyymmdd_HHMMSS")
    output_dir = joinpath(pwd(), "outputs")
    mkpath(output_dir)
    video_path = joinpath(output_dir, "case_eta_$(tag).mp4")

    render_eta_video(history, x_full, y_full, video_path; fps=20, title_prefix="eta")
    println(video_path)
    return nothing
end

main()

# vertical-dynamics

`vertical-dynamics` is a Julia testbed for non-wetting impact models on a 2D surface grid. The current codebase focuses on a non-axisymmetric bath with a prescribed contact mask, a rigid-plane contact model, and structured Dirichlet-to-Neumann operators for the fluid coupling.

## Current Scope

The repository currently includes:

- a corrected nonperiodic DtN operator on a uniform Cartesian grid
- a prescribed-contact rigid-plane solver for non-axisymmetric bath dynamics
- dense small-grid reference assembly and a matrix-free GMRES solve path
- a demo script that runs one case and renders a top-view heatmap video of `η(x,y,t)`

The current solver is intended as a research testbed, not a fully validated production code.

## Repository Map

Start with:

- [docs/project-map.md](/Users/eaguerov/Documents/Github/vertical-dynamics/docs/project-map.md)
- [docs/source-tree.md](/Users/eaguerov/Documents/Github/vertical-dynamics/docs/source-tree.md)
- [docs/test-guide.md](/Users/eaguerov/Documents/Github/vertical-dynamics/docs/test-guide.md)

## Requirements

- Julia `1.12`
- `ffmpeg` in `PATH` if you want to render MP4 videos

Julia dependencies are managed through the local project:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## Running Tests

Run the full test suite with:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Running the Demo Case

Run one demo case and produce an MP4 heatmap video with:

```bash
julia --project=. scripts/run_case_and_make_video.jl
```

The script writes output under `outputs/` and prints the final MP4 path.

## Customizing a Case

The easiest way to run your own case is to copy or edit:

- [scripts/run_case_and_make_video.jl](/Users/eaguerov/Documents/Github/vertical-dynamics/scripts/run_case_and_make_video.jl)

That script is the current reference example for how to:

- build a `contact_mask`
- set grid and timestep parameters such as `nx`, `ny`, `dx`, `dt`, and `steps`
- choose fluid/body parameters such as `Fr`, `We`, `Re`, and `body_damping`
- set the mass distribution through `mass_density`
- choose the initial rigid-plane state through `SolverState(..., q, v)`
- call `solve_motion(cfg; initial_state=state0)`
- render the result with `render_eta_video(...)`

In practice, the main block to edit is:

```julia
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
```

If you only want to change parameters, start from that script rather than calling low-level functions directly.

## Current Model Assumptions

The current coupled model assumes:

- the contact mask is prescribed
- the contact region is constrained to remain coplanar
- the contact pressure is solved as an unknown field on the mask
- the body motion is reduced to rigid-plane generalized coordinates
- the bath is discretized on a uniform Cartesian grid

The code does not yet solve for an unknown contact region.

## Example Output

After running the demo script, a typical output file looks like:

- [outputs/case_eta_20260318_095813.mp4](/Users/eaguerov/Documents/Github/vertical-dynamics/outputs/case_eta_20260318_095813.mp4)

## Notes

- The periodic DtN oracle was intentionally removed from the active API to keep the codebase focused on the nonperiodic path.
- The visualization path uses `Plots.jl` for frame rendering and `ffmpeg` for MP4 assembly.

# Master Roadmap

This is the internal planning document for the repository. It is meant to capture the current project vision, what has already been built, and what should happen next. It does not override the source code or tests.

## Project Vision

The long-term goal is a research-grade Julia codebase for fluid-bath impact problems in which the bath dynamics are fully non-axisymmetric, the fluid coupling is handled through a structured Dirichlet-to-Neumann operator, and the solver architecture remains scalable enough for realistic 2D surface grids.

The current direction is deliberately staged:

1. establish a reliable nonperiodic DtN operator
2. build a first prescribed-contact non-axisymmetric solver
3. make the linear algebra scalable through matrix-free Krylov methods
4. add enough diagnostics and visualization to run and inspect actual cases
5. later move toward unknown contact sets, stronger validation, and richer physics

## Current Implemented State

### DtN Layer

The active DtN path is a corrected nonperiodic kernel operator on a uniform Cartesian grid. It uses:

- real-space singular-kernel construction
- local correction near the origin
- FFT-based linear convolution for application

The older periodic oracle has been removed from the active API to keep the codebase focused on the nonperiodic route.

### Solver Layer

The implemented solver is a first prescribed-contact rigid-plane model. It advances:

- free-surface elevation `η(x,y,t)`
- surface potential `φ(x,y,t)`
- contact pressure on a prescribed mask
- rigid-plane generalized coordinates `q = (a, b, c)` and velocities

The step solve supports:

- dense assembled reference solves on small grids
- matrix-free GMRES for the scalable path

The rigid-plane body model also includes:

- generalized mass matrix from a prescribed mass-density field
- linear surface-tension forcing induced by the contact-mask boundary

### Visualization Layer

The repository now includes a minimal visualization workflow:

- `src/visualization.jl` computes global color limits and renders `η` heatmap videos
- `scripts/run_case_and_make_video.jl` runs one representative case and writes an MP4 through `ffmpeg`

This is intended as a practical inspection tool, not a full post-processing framework.

## Completed Milestones

- [x] Study the Galeano-Rios 2017 formulation and extract the operative boundary/contact system.
- [x] Separate the three main design axes: vorticity elimination, symmetry assumptions, and DtN implementation.
- [x] Choose the real-space corrected-kernel DtN route as the active implementation direction.
- [x] Implement the corrected nonperiodic DtN operator in Julia.
- [x] Add direct and FFT-based DtN application paths and verify their consistency.
- [x] Build a first prescribed-contact non-axisymmetric finite-difference solver.
- [x] Keep contact pressure as an unknown reaction field rather than prescribing its shape.
- [x] Add generalized rigid-plane body coordinates and nonuniform mass-density support.
- [x] Add linear contact-boundary surface-tension forcing.
- [x] Replace dense full-system solves with a matrix-free GMRES path while retaining dense assembly as a correctness reference.
- [x] Add project-level documentation and source-level docstrings.
- [x] Add a demo script and MP4 visualization path for top-view `η` output.
- [x] Remove the periodic DtN oracle from the active API and test suite.

## Immediate Next Milestones

- [ ] Create a clean preconditioner for the matrix-free GMRES solve so larger grids become practical.
- [ ] Define one or two benchmark cases with expected qualitative behavior and use them to sanity-check pressure, force, and rigid-plane dynamics.
- [ ] Add a cleaner output/data product path beyond raw video, for example saved histories or lightweight diagnostics files.
- [ ] Decide whether the next modeling step should be unknown contact-set dynamics or stronger viscous correction.
- [ ] Refactor the demo script into a more reusable case-definition pattern if multiple scenarios are about to be added.

## Medium-Term Research And Engineering Goals

- [ ] Move from prescribed contact masks to an unknown contact-set formulation.
- [ ] Design the contact update as an active-set or complementarity problem rather than a fixed mask.
- [ ] Add a reusable sparse preconditioner that can be cached for fixed `dt` and fixed contact topology.
- [ ] Run mesh and timestep studies on the full coupled solver, not only on the isolated DtN operator.
- [ ] Compare the non-axisymmetric solver behavior against the neighboring `km-*` repositories on matched simplified cases.
- [ ] Decide whether the viscous upgrade should be a boundary-layer correction or a more radical bulk reformulation.

## Open Questions

### Contact Modeling

The largest missing piece is the contact-set update. The current solver is square and operational only because the mask is prescribed. The moment the contact area becomes unknown, the formulation will need an explicit active-set strategy or a related constrained solve.

### Viscosity

The current weakly viscous reduction is good enough to support the first solver architecture, but it is not the final word. The main unresolved question is whether the next viscosity improvement should stay compatible with the boundary-only DtN framework or whether the underlying model class needs to change.

### Validation

The code is now structurally predictive for the prescribed-contact model, but it is not yet scientifically validated. That validation step needs benchmark cases, diagnostic checks, and comparisons that go beyond “the code runs and the tests pass.”

## Working Rule

When the architecture changes materially, update this file first. More detailed implementation notes can stay in the narrower plan and memo files, but this document should remain the best single internal snapshot of where the project stands.

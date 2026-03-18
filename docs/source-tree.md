# Source Tree Guide

## Package Entrypoint

### [src/VerticalDynamics.jl](/Users/eaguerov/Documents/Github/vertical-dynamics/src/VerticalDynamics.jl)

Purpose: loads the DtN and solver subsystems and exports the public API.

Key exported symbols:

- DtN: `CorrectedKernelDtN2D`, `build_corrected_kernel_dtn`, `build_toeplitz_dtn`, `apply`, `apply_direct`, `grid_coordinates`
- Solver: `SolverConfig`, `SolverDomain`, `SolverSystem`, `SolverState`, `build_solver`, `build_solver_domain`, `build_dense_dtn_matrix`, `assemble_step_system`, `step_layout`, `build_step_rhs`, `apply_step_operator`, `apply_step_operator!`, `gmres_solve`, `advance_one_step`, `solve_motion`, `zero_state`

## DtN Subsystem

The `src/dtn/` files implement structured Dirichlet-to-Neumann operators for a uniform 2D grid.

### [src/dtn/toeplitz_operator.jl](/Users/eaguerov/Documents/Github/vertical-dynamics/src/dtn/toeplitz_operator.jl)

Purpose: defines the DtN operator types and constructs their kernels.

Key functions and types:

- `CorrectedKernelDtN2D`
  Nonperiodic corrected kernel with FFT-based linear convolution support.
- `build_corrected_kernel_dtn(nx, ny, dx; near_radius, quadrature_order)`
- `build_toeplitz_dtn(...)`
  Current default alias for the corrected nonperiodic operator.
- `grid_coordinates(nx, ny, dx)`

Internal focus:

- origin-cell correction from quadratic Taylor exactness
- corrected axis-neighbor weights
- kernel FFT embedding for linear convolution

### [src/dtn/apply.jl](/Users/eaguerov/Documents/Github/vertical-dynamics/src/dtn/apply.jl)

Purpose: applies DtN operators either directly or by FFT.

Key functions:

- `apply_direct(op::CorrectedKernelDtN2D, x)`
- `apply(op::CorrectedKernelDtN2D, x)`

Direct dependencies:

- operator types and stored kernels from [src/dtn/toeplitz_operator.jl](/Users/eaguerov/Documents/Github/vertical-dynamics/src/dtn/toeplitz_operator.jl)

## Solver Subsystem

The `src/solver/` files implement a prescribed-contact rigid-plane solver for a 2D bath.

### [src/solver/types.jl](/Users/eaguerov/Documents/Github/vertical-dynamics/src/solver/types.jl)

Purpose: defines the solver configuration and state containers.

Key types:

- `SolverConfig`
- `SolverDomain`
- `SolverSystem`
- `SolverState`

Key functions:

- `SolverConfig(; ...)`
- `zero_state(system)`

### [src/solver/domain.jl](/Users/eaguerov/Documents/Github/vertical-dynamics/src/solver/domain.jl)

Purpose: builds spatial operators and geometry-dependent maps.

Key functions:

- `build_solver(config)`
- `build_solver_domain(config)`
- `build_dense_dtn_matrix(nx, ny, dx; ...)`

Internal responsibilities:

- Cartesian 5-point Laplacian
- contact/free index maps
- rigid-plane basis over the contact mask
- generalized body mass matrix from `mass_density`
- linear contact-boundary surface-tension map

### [src/solver/assembly.jl](/Users/eaguerov/Documents/Github/vertical-dynamics/src/solver/assembly.jl)

Purpose: assembles the dense implicit Euler step matrix for correctness checks and small-grid reference solves.

Key function:

- `assemble_step_system(system, state)`

### [src/solver/operator.jl](/Users/eaguerov/Documents/Github/vertical-dynamics/src/solver/operator.jl)

Purpose: defines the matrix-free timestep operator.

Key functions:

- `step_layout(system)`
- `build_step_rhs(system, state)`
- `apply_step_operator(system, x)`
- `apply_step_operator!(y, system, x)`

This file is the core scalable solver layer: it combines sparse finite-difference terms, FFT-based DtN application, contact-pressure injection, and rigid-plane coupling without assembling the full matrix.

### [src/solver/krylov.jl](/Users/eaguerov/Documents/Github/vertical-dynamics/src/solver/krylov.jl)

Purpose: provides the current Krylov linear solver.

Key function:

- `gmres_solve(system, b; x0=nothing)`

Current limitation:

- no preconditioner yet

### [src/solver/time_marching.jl](/Users/eaguerov/Documents/Github/vertical-dynamics/src/solver/time_marching.jl)

Purpose: advances the coupled model in time.

Key functions:

- `advance_one_step(state, system)`
- `solve_motion(config; initial_state=nothing)`

Behavior:

- uses `:direct` or `:gmres` according to `SolverConfig.linear_solver`
- reconstructs `η` and pressure fields from reduced unknown vectors after each solve

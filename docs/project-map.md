# Project Map

This repository is a Julia package for building and testing free-surface impact models. The current implementation has two main layers:

- `src/dtn/`: structured Dirichlet-to-Neumann operators for the half-space free surface
- `src/solver/`: a prescribed-contact finite-difference solver with rigid-plane body kinematics

The project currently supports a fixed contact mask. It does not yet solve for the contact set as an unknown.

## Main Directories

- [src/VerticalDynamics.jl](/Users/eaguerov/Documents/Github/vertical-dynamics/src/VerticalDynamics.jl)
  Package entrypoint and exports.

- [src/dtn/toeplitz_operator.jl](/Users/eaguerov/Documents/Github/vertical-dynamics/src/dtn/toeplitz_operator.jl)
  DtN operator types and constructors.

- [src/dtn/apply.jl](/Users/eaguerov/Documents/Github/vertical-dynamics/src/dtn/apply.jl)
  Direct and FFT-based DtN application routines.

- [src/solver/types.jl](/Users/eaguerov/Documents/Github/vertical-dynamics/src/solver/types.jl)
  Solver configuration, domain, system, and state types.

- [src/solver/domain.jl](/Users/eaguerov/Documents/Github/vertical-dynamics/src/solver/domain.jl)
  Cartesian-grid geometry, Laplacian, mask maps, body mass matrix, and boundary-force maps.

- [src/solver/assembly.jl](/Users/eaguerov/Documents/Github/vertical-dynamics/src/solver/assembly.jl)
  Dense assembled step matrix used as a correctness reference.

- [src/solver/operator.jl](/Users/eaguerov/Documents/Github/vertical-dynamics/src/solver/operator.jl)
  Matrix-free block operator for the timestep solve.

- [src/solver/krylov.jl](/Users/eaguerov/Documents/Github/vertical-dynamics/src/solver/krylov.jl)
  GMRES implementation for matrix-free step solves.

- [src/solver/time_marching.jl](/Users/eaguerov/Documents/Github/vertical-dynamics/src/solver/time_marching.jl)
  One-step advance and `solve_motion` loop.

- [test](/Users/eaguerov/Documents/Github/vertical-dynamics/test)
  DtN and solver verification suite.

- [docs/plans](/Users/eaguerov/Documents/Github/vertical-dynamics/docs/plans)
  Research, design, and implementation notes.

- [.requirements](/Users/eaguerov/Documents/Github/vertical-dynamics/.requirements)
  Per-task requirement snapshots from development turns.

## Current Solver Model

The implemented solver advances:

- free-surface elevation `η(x,y,t)`
- surface potential `φ(x,y,t)`
- contact pressure over a prescribed mask
- rigid-plane body variables `q = (a,b,c)` and velocities

with:

- 2D Cartesian finite differences
- implicit Euler time stepping
- FFT-based DtN application
- dense direct solve or matrix-free GMRES

The surface-tension contribution from the contact-mask boundary is also included as a linear generalized load on the rigid-plane degrees of freedom.

## Reference vs Scalable Paths

Two solve paths coexist on purpose:

- dense assembled reference:
  [src/solver/assembly.jl](/Users/eaguerov/Documents/Github/vertical-dynamics/src/solver/assembly.jl)

- scalable matrix-free path:
  [src/solver/operator.jl](/Users/eaguerov/Documents/Github/vertical-dynamics/src/solver/operator.jl)
  and
  [src/solver/krylov.jl](/Users/eaguerov/Documents/Github/vertical-dynamics/src/solver/krylov.jl)

The dense path is used to verify correctness on small grids. The matrix-free path is the intended direction for larger problems.

## Related Docs

- [docs/source-tree.md](/Users/eaguerov/Documents/Github/vertical-dynamics/docs/source-tree.md)
- [docs/test-guide.md](/Users/eaguerov/Documents/Github/vertical-dynamics/docs/test-guide.md)
- [docs/plans/README.md](/Users/eaguerov/Documents/Github/vertical-dynamics/docs/plans/README.md)

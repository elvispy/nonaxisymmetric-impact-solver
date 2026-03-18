# As Is

The repository currently contains a Julia package with Dirichlet-to-Neumann operators and tests for periodic and nonperiodic application. There is no time-dependent impact solver, no finite-difference bath discretization beyond the DtN operator, no rigid-body/contact model, and no analogue of `solve_motion.m` from `km-disk-vibrating`.

# To Be

The repository should contain a first non-axisymmetric finite-difference solver in Julia that mirrors the structure of `km-disk-vibrating/matlab/1_code/simulation_code/solve_motion.m`. The bath will be discretized on a 2D Cartesian grid, the contact mask will be prescribed, the contact surface will be constrained to lie on a rigid plane, the contact pressure will remain an unknown field over the mask, and the time step will be implicit Euler. The body will support prescribed non-uniform mass distributions from the start through a generalized rigid-plane inertia matrix.

# Requirements

1. The package shall expose a solver entry point analogous to `solve_motion`.
2. The solver shall discretize the bath on a uniform Cartesian grid with finite differences for the horizontal Laplacian.
3. The solver shall accept a prescribed contact mask and enforce a rigid-plane constraint on that mask.
4. The solver shall treat the contact pressure on the prescribed mask as an unknown field solved implicitly at each step.
5. The solver shall include rigid-plane body dynamics with generalized coordinates `q = (a, b, c)` and generalized velocities.
6. The solver shall support prescribed non-uniform mass distributions through a generalized mass matrix derived from the grid density field.
7. The solver shall assemble one square implicit Euler linear system per time step.
8. The solver shall provide a one-step advance function and a multi-step `solve_motion` loop.
9. The package shall support a matrix-free step operator using FFT-based DtN application and a Krylov linear solve.
10. The package shall retain the dense assembled solve as a correctness reference for small grids.
11. The package shall include tests for system dimension consistency, rigid-plane enforcement, body mass-matrix construction, operator consistency, and zero-solution invariance.

# Acceptance Criteria

1.1 A function `solve_motion(config)` is exported and returns a history of states.
1.2 The solver code is organized in dedicated solver files rather than folded into the DtN modules.

2.1 A domain builder constructs Cartesian coordinates and a 5-point Laplacian on the full 2D grid.
2.2 The Laplacian acts on flattened grid vectors with consistent indexing.

3.1 The contact mask is represented explicitly on the grid.
3.2 At every advanced state, the contact-point values of `η` equal an affine plane `a + bx + cy`.

4.1 The step solve includes one pressure unknown per contact grid point.
4.2 Pressure is zero outside the prescribed contact mask by construction.

5.1 The body unknowns include `q_next` and `v_next`, each of length 3.
5.2 The rigid-plane kinematic update `q_next = q_prev + dt * v_next` is enforced in the linear system.

6.1 A non-uniform density field over the body grid produces the generalized mass matrix
`M = ∫ ρ [1, x, y]^T [1, x, y] dA`.
6.2 For a centered uniform density field, the mixed first moments vanish numerically.

7.1 The assembled step matrix is square.
7.2 The matrix row/column count matches `(N - M) + N + M + 3 + 3 = 2N + 6`, where `N` is the total number of bath nodes and `M` the number of contact nodes.

8.1 `advance_one_step(state, system)` returns a new state with updated `η`, `φ`, `p`, `q`, `v`, and time.
8.2 `solve_motion(config)` repeatedly applies `advance_one_step`.

9.1 A helper applies the DtN block through the existing FFT-based operator without assembling it densely.
9.2 A matrix-free `A*v` implementation matches the assembled reference operator on small grids.

10.1 The dense assembled solve remains available for small-grid verification.
10.2 The matrix-free Krylov solve is the default solve path.

11.1 Tests pass for the solver API and system sizing.
11.2 A zero initial state with zero gravity and zero forcing remains zero after one step.
11.3 A one-step solve enforces the plane constraint on the prescribed contact mask.
11.4 Matrix-free operator application matches the assembled operator on a test vector to numerical tolerance.

# Testing Plan

1. Add solver API tests covering domain/system construction and exported entry points.
2. Add tests for generalized rigid-body mass matrix moments using centered density fields.
3. Add a test that the assembled implicit system is square with the expected dimension count.
4. Add a one-step test verifying that contact values of `η` lie exactly on the solved rigid plane.
5. Add a zero-invariance test with no forcing, no gravity, and zero initial state.
6. Run the full Julia test suite after integrating the solver modules.

# Implementation Plan

1. Create solver module files for types, grid/domain construction, dense DtN matrix assembly, step assembly, and time marching.
   Test: add an API smoke test that imports and constructs a solver configuration.
2. Implement Cartesian grid and Laplacian construction plus contact/body geometry helpers.
   Test: verify dimensions and centered-moment properties.
3. Implement the implicit Euler system assembly for the prescribed-mask rigid-plane formulation.
   Test: verify the matrix is square with size `2N + 6`.
4. Implement `advance_one_step` and `solve_motion`.
   Test: verify zero-solution invariance and rigid-plane enforcement.
5. Export the new solver API and update `runtests.jl`.
   Test: run `julia --project=. -e 'using Pkg; Pkg.test()'`.

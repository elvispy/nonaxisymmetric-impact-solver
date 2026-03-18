# As Is

The repository contains a Julia DtN operator implementation under `src/dtn/` that builds the exact periodic Fourier-symbol operator on a rectangular grid. The exported API `build_toeplitz_dtn(nx, ny, dx)` currently returns this periodic operator, and `apply(op, x)` performs periodic FFT application. The existing tests verify periodic properties such as exact Fourier-mode action and machine-precision agreement between direct circular convolution and FFT application.

# To Be

The main API should construct a nonperiodic, zero-extended 2D half-space DtN operator based on the real-space singular-integral kernel. The singular origin-cell contribution should be corrected from quadratic Taylor exactness, yielding a discrete operator that is second-order accurate for smooth fields that decay to near zero at the computational boundaries. The current periodic operator should remain available as a reference/oracle, but it should no longer be the default `build_toeplitz_dtn` path.

# Requirements

1. The package shall expose a nonperiodic corrected-kernel operator as the default `build_toeplitz_dtn(nx, ny, dx)` result.
2. The package shall retain a periodic reference operator accessible through a distinct API.
3. The nonperiodic operator shall use zero-extended Toeplitz application rather than periodic wrap-around.
4. The nonperiodic operator shall include an explicit singular origin-cell correction derived from quadratic exactness.
5. The nonperiodic operator shall preserve the constant nullspace exactly.
6. The package shall provide both direct and FFT-based application paths for the nonperiodic operator.
7. The test suite shall distinguish the periodic reference behavior from the nonperiodic corrected-kernel behavior.
8. The nonperiodic operator shall demonstrate approximately second-order convergence on smooth fields that decay rapidly near the computational boundaries, using the periodic oracle only as a validation reference.

# Acceptance Criteria

1.1 `build_toeplitz_dtn(nx, ny, dx)` returns the nonperiodic corrected-kernel type.
1.2 Applying the default operator to a compactly supported or rapidly decaying field does not wrap values across opposite boundaries.

2.1 A separate constructor such as `build_periodic_dtn(nx, ny, dx)` is exported.
2.2 The periodic constructor preserves the current exact-symbol behavior for Fourier-mode tests.

3.1 `apply(op, x)` for the corrected-kernel operator agrees with `apply_direct(op, x)` up to numerical tolerance on small grids.
3.2 An impulse located near one boundary does not create mirrored contributions on the opposite boundary under the corrected-kernel operator.

4.1 The origin-cell correction coefficient is implemented explicitly from a quadratic Taylor argument.
4.2 The corrected-kernel operator includes this correction through a local stencil contribution.

5.1 Applying the corrected-kernel operator to a constant field yields zero to machine precision.

6.1 The corrected-kernel operator supports a direct reference path and an FFT-accelerated path.
6.2 The FFT path uses linear convolution with zero extension rather than circular convolution.

7.1 Tests for exact Fourier-mode action target the periodic operator only.
7.2 Tests for nonperiodic behavior target the corrected-kernel operator only.

8.1 A refinement study on a smooth decaying field shows observed convergence near order 2 against the periodic-oracle reference on a sufficiently large box.

# Testing Plan

1. Add API tests covering the default corrected-kernel constructor and the separate periodic constructor.
2. Add nullspace and no-wrap tests for the corrected-kernel operator.
3. Add FFT-vs-direct consistency tests for the corrected-kernel operator on small grids.
4. Keep a periodic Fourier-mode exactness test for the periodic operator.
5. Add a corrected-kernel convergence test on a smooth Gaussian field over a large fixed box, using the periodic oracle as reference.

# Implementation Plan

1. Split the current operator into periodic and nonperiodic types, export both constructors, and update API tests.
   Test: run API and periodic Fourier-mode tests.
2. Implement nonperiodic kernel construction with midpoint far field, high-order near-cell quadrature, and symbolic origin-cell correction.
   Test: run nullspace and kernel-structure tests.
3. Implement direct zero-extended application for the nonperiodic operator.
   Test: add and run no-wrap and impulse tests.
4. Implement FFT-based linear convolution for the nonperiodic operator and verify against direct application.
   Test: run FFT-vs-direct consistency tests.
5. Add convergence tests on smooth decaying fields against the periodic oracle and tune near-field defaults if needed.
   Test: run the full Julia test suite.

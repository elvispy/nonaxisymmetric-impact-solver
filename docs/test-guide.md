# Test Guide

Run the full suite with:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

A green run means:

- DtN operator construction and application remain internally consistent
- the prescribed-contact solver still matches its dense and matrix-free formulations on the covered small-grid cases

## DtN Tests

### [test/test_api.jl](/Users/eaguerov/Documents/Github/vertical-dynamics/test/test_api.jl)

- Covers: [src/dtn/toeplitz_operator.jl](/Users/eaguerov/Documents/Github/vertical-dynamics/src/dtn/toeplitz_operator.jl)
- Verifies: constructor availability, returned types, kernel dimensions

### [test/test_kernel_structure.jl](/Users/eaguerov/Documents/Github/vertical-dynamics/test/test_kernel_structure.jl)

- Covers: [src/dtn/toeplitz_operator.jl](/Users/eaguerov/Documents/Github/vertical-dynamics/src/dtn/toeplitz_operator.jl)
- Verifies: symmetry, sign structure, zero-row-sum kernel property

### [test/test_impulse_response.jl](/Users/eaguerov/Documents/Github/vertical-dynamics/test/test_impulse_response.jl)

- Covers: [src/dtn/apply.jl](/Users/eaguerov/Documents/Github/vertical-dynamics/src/dtn/apply.jl)
- Verifies: nonperiodic impulse response behavior

### [test/test_fft_apply.jl](/Users/eaguerov/Documents/Github/vertical-dynamics/test/test_fft_apply.jl)

- Covers: [src/dtn/apply.jl](/Users/eaguerov/Documents/Github/vertical-dynamics/src/dtn/apply.jl)
- Verifies: FFT and direct application agree

### [test/test_convergence.jl](/Users/eaguerov/Documents/Github/vertical-dynamics/test/test_convergence.jl)

- Covers: [src/dtn/toeplitz_operator.jl](/Users/eaguerov/Documents/Github/vertical-dynamics/src/dtn/toeplitz_operator.jl), [src/dtn/apply.jl](/Users/eaguerov/Documents/Github/vertical-dynamics/src/dtn/apply.jl)
- Verifies: corrected-kernel self-convergence at the center on a fixed domain

## Solver Tests

### [test/test_solver_api.jl](/Users/eaguerov/Documents/Github/vertical-dynamics/test/test_solver_api.jl)

- Covers: [src/solver/types.jl](/Users/eaguerov/Documents/Github/vertical-dynamics/src/solver/types.jl), [src/solver/domain.jl](/Users/eaguerov/Documents/Github/vertical-dynamics/src/solver/domain.jl)
- Verifies: solver construction and basic domain bookkeeping

### [test/test_solver_domain.jl](/Users/eaguerov/Documents/Github/vertical-dynamics/test/test_solver_domain.jl)

- Covers: [src/solver/domain.jl](/Users/eaguerov/Documents/Github/vertical-dynamics/src/solver/domain.jl)
- Verifies: generalized mass moments, plane basis rank, and zero boundary-force response on constant `η`

### [test/test_solver_assembly.jl](/Users/eaguerov/Documents/Github/vertical-dynamics/test/test_solver_assembly.jl)

- Covers: [src/solver/assembly.jl](/Users/eaguerov/Documents/Github/vertical-dynamics/src/solver/assembly.jl)
- Verifies: dense step matrix is square with the expected size

### [test/test_solver_operator.jl](/Users/eaguerov/Documents/Github/vertical-dynamics/test/test_solver_operator.jl)

- Covers: [src/solver/operator.jl](/Users/eaguerov/Documents/Github/vertical-dynamics/src/solver/operator.jl), [src/solver/assembly.jl](/Users/eaguerov/Documents/Github/vertical-dynamics/src/solver/assembly.jl)
- Verifies: matrix-free operator action matches the assembled operator on a test vector

### [test/test_solver_step.jl](/Users/eaguerov/Documents/Github/vertical-dynamics/test/test_solver_step.jl)

- Covers: [src/solver/time_marching.jl](/Users/eaguerov/Documents/Github/vertical-dynamics/src/solver/time_marching.jl), [src/solver/krylov.jl](/Users/eaguerov/Documents/Github/vertical-dynamics/src/solver/krylov.jl)
- Verifies: zero-state invariance, rigid-plane enforcement, pressure locality, and GMRES-vs-direct one-step agreement

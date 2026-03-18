# Julia Toeplitz DtN Operator Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a Julia package that constructs and applies a 2D half-space Dirichlet-to-Neumann operator using a real-space block-Toeplitz kernel with FFT-accelerated application and tests for second-order behavior on square grids.

**Architecture:** The package will separate three concerns: kernel construction, FFT/circulant embedding for fast operator application, and verification utilities. The first working slice will prioritize a stable API and correctness properties (`N*1=0`, symmetry, shape), then extend to mode-based convergence tests and near-field correction logic.

**Tech Stack:** Julia 1.12, `FFTW`, Julia `Test` stdlib, lightweight package skeleton under `src/` and `test/`.

---

## Chunk 1: Package Skeleton And Operator API

### Task 1: Create the Julia package layout

**Files:**
- Create: `Project.toml`
- Create: `src/VerticalDynamics.jl`
- Create: `src/dtn/toeplitz_operator.jl`
- Create: `test/runtests.jl`
- Create: `test/test_api.jl`

- [ ] **Step 1: Write the failing test**

```julia
using Test
using VerticalDynamics

@testset "Toeplitz DtN API" begin
    op = build_toeplitz_dtn(8, 8, 0.1)
    @test size(op.kernel) == (8, 8)
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: FAIL because package files and `build_toeplitz_dtn` do not exist yet.

- [ ] **Step 3: Write minimal implementation**

Create a minimal module exporting `build_toeplitz_dtn` and a placeholder operator type with a stored kernel matrix.

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: PASS for API construction test.

- [ ] **Step 5: Commit**

```bash
git add Project.toml src test
git commit -m "feat: scaffold Julia Toeplitz DtN package"
```

### Task 2: Enforce constant-nullspace behavior

**Files:**
- Modify: `src/dtn/toeplitz_operator.jl`
- Create: `src/dtn/apply.jl`
- Modify: `src/VerticalDynamics.jl`
- Create: `test/test_nullspace.jl`

- [ ] **Step 1: Write the failing test**

```julia
using Test
using VerticalDynamics

@testset "Constants map to zero" begin
    op = build_toeplitz_dtn(16, 16, 0.1)
    x = ones(16, 16)
    y = apply(op, x)
    @test maximum(abs, y) < 1e-12
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: FAIL because `apply` is missing or constants do not map to zero.

- [ ] **Step 3: Write minimal implementation**

Implement kernel storage and `apply(op, x)` using a direct convolution placeholder or padded FFT placeholder. Define the self-term by zero row sum so constants map exactly to zero.

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project=. -e 'using Pkg; Pkg.test(test_args=["nullspace"])'`
Expected: PASS for nullspace behavior.

- [ ] **Step 5: Commit**

```bash
git add src test
git commit -m "feat: add Toeplitz DtN apply and nullspace condition"
```

## Chunk 2: Kernel Construction And Symmetry

### Task 3: Build a symmetric far-field kernel

**Files:**
- Modify: `src/dtn/toeplitz_operator.jl`
- Create: `src/dtn/kernel_construction.jl`
- Modify: `src/VerticalDynamics.jl`
- Create: `test/test_symmetry.jl`

- [ ] **Step 1: Write the failing test**

```julia
using Test
using VerticalDynamics

@testset "Kernel symmetry" begin
    op = build_toeplitz_dtn(17, 17, 0.1)
    K = op.kernel
    @test K == reverse(K, dims=1)
    @test K == reverse(K, dims=2)
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: FAIL because kernel construction is incomplete or asymmetric.

- [ ] **Step 3: Write minimal implementation**

Construct the reference kernel from offset geometry on a square grid, using symmetric far-field weights and an exact zero-row-sum self-term.

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project=. -e 'using Pkg; Pkg.test(test_args=["symmetry"])'`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src test
git commit -m "feat: add symmetric far-field DtN kernel"
```

### Task 4: Add circulant embedding and cached FFT application

**Files:**
- Create: `src/dtn/fft_apply.jl`
- Modify: `src/dtn/apply.jl`
- Modify: `src/dtn/toeplitz_operator.jl`
- Modify: `src/VerticalDynamics.jl`
- Create: `test/test_fft_apply.jl`

- [ ] **Step 1: Write the failing test**

```julia
using Test
using VerticalDynamics

@testset "FFT apply matches direct apply" begin
    op = build_toeplitz_dtn(8, 8, 0.2)
    x = reshape(collect(1.0:64.0), 8, 8)
    y_direct = apply_direct(op, x)
    y_fft = apply(op, x)
    @test y_fft ≈ y_direct atol=1e-10 rtol=1e-10
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: FAIL because FFT path is missing or disagrees.

- [ ] **Step 3: Write minimal implementation**

Implement circulant embedding, cache the kernel FFT, and route `apply` through FFT multiplication while retaining `apply_direct` as a small-grid reference.

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project=. -e 'using Pkg; Pkg.test(test_args=["fft_apply"])'`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src test
git commit -m "feat: add FFT-accelerated Toeplitz DtN apply"
```

## Chunk 3: Near-Field Correction And Accuracy

### Task 5: Introduce near-field correction hooks

**Files:**
- Create: `src/dtn/nearfield_correction.jl`
- Modify: `src/dtn/kernel_construction.jl`
- Modify: `src/VerticalDynamics.jl`
- Create: `test/test_nearfield.jl`

- [ ] **Step 1: Write the failing test**

```julia
using Test
using VerticalDynamics

@testset "Near-field correction metadata" begin
    op = build_toeplitz_dtn(32, 32, 0.1; near_radius=2)
    @test op.near_radius == 2
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: FAIL because near-field correction support is missing.

- [ ] **Step 3: Write minimal implementation**

Add data structures and hooks for corrected local weights, leaving the first correction as a no-op placeholder if necessary.

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project=. -e 'using Pkg; Pkg.test(test_args=["nearfield"])'`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src test
git commit -m "feat: add near-field correction framework"
```

### Task 6: Validate Fourier-mode action on smooth fields

**Files:**
- Create: `test/test_fourier_modes.jl`
- Create: `src/dtn/test_functions.jl`
- Modify: `src/VerticalDynamics.jl`

- [ ] **Step 1: Write the failing test**

```julia
using Test
using VerticalDynamics

@testset "Fourier mode check" begin
    nx = 32
    dx = 2π / nx
    op = build_toeplitz_dtn(nx, nx, dx)
    kx, ky = 2, 1
    x, y = grid_coordinates(nx, nx, dx)
    f = @. cos(kx * x + ky * y)
    target = sqrt(kx^2 + ky^2) .* f
    got = apply(op, f)
    @test norm(got - target) / norm(target) < 5e-2
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: FAIL because the operator is not yet accurate enough.

- [ ] **Step 3: Write minimal implementation**

Add mode-generation utilities and improve near-field weights until the smooth-mode test passes at coarse resolution.

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project=. -e 'using Pkg; Pkg.test(test_args=["fourier_modes"])'`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src test
git commit -m "feat: validate DtN action on Fourier modes"
```

### Task 7: Prove second-order convergence on square grids

**Files:**
- Create: `test/test_convergence.jl`
- Modify: `src/dtn/test_functions.jl`

- [ ] **Step 1: Write the failing test**

```julia
using Test
using VerticalDynamics

@testset "Second-order convergence" begin
    sizes = [16, 32, 64]
    errors = Float64[]
    for n in sizes
        dx = 2π / n
        op = build_toeplitz_dtn(n, n, dx)
        x, y = grid_coordinates(n, n, dx)
        f = @. cos(2x + y)
        target = sqrt(5.0) .* f
        got = apply(op, f)
        push!(errors, norm(got - target) / norm(target))
    end
    rate = log(errors[1] / errors[end]) / log(4)
    @test rate > 1.8
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: FAIL because the current near-field treatment does not yet achieve quadratic convergence.

- [ ] **Step 3: Write minimal implementation**

Tune or redesign the near-field correction so the observed convergence rate exceeds `1.8` on square grids.

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project=. -e 'using Pkg; Pkg.test(test_args=["convergence"])'`
Expected: PASS with observed rate near 2.

- [ ] **Step 5: Commit**

```bash
git add src test
git commit -m "feat: achieve second-order Toeplitz DtN convergence"
```

## Chunk 4: Documentation And Usage

### Task 8: Document operator assumptions and fallback paths

**Files:**
- Create: `README.md`
- Modify: `docs/plans/julia-toeplitz-dtn-design-brief.md`

- [ ] **Step 1: Write the failing test**

There is no executable test for this task. Instead, verify documentation coverage manually against the design brief.

- [ ] **Step 2: Run manual verification**

Check that the README includes:
- purpose of the operator
- square-grid assumption
- second-order target
- near-field correction strategy
- rollback alternatives

- [ ] **Step 3: Write minimal implementation**

Add a short package README and update the brief if implementation details diverge.

- [ ] **Step 4: Run manual verification**

Re-read the README and confirm it matches the implemented API.

- [ ] **Step 5: Commit**

```bash
git add README.md docs/plans/julia-toeplitz-dtn-design-brief.md
git commit -m "docs: describe Toeplitz DtN package and assumptions"
```

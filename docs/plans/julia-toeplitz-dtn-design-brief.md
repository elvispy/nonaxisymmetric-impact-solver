# Julia design brief: block-Toeplitz DtN operator

## Objective

Implement the half-space Dirichlet-to-Neumann operator in Julia for a scalar field on a uniform 2D boundary grid, using the real-space singular-kernel formulation and FFT-accelerated application. The first target is not the full impact solver, only a reliable operator `y = N*x`.

The accuracy target is explicit: for smooth test fields on a square grid with `hx = hy = dx`, the discrete operator should exhibit quadratic convergence with respect to `dx`.

## Governing operator

For a boundary field `f(r)` on the plane,

```text
Nf(r) = (1 / 2π) PV ∫ (f(r) - f(s)) / |r - s|^3 dA(s).
```

On a uniform Cartesian grid, the discrete operator should inherit translation invariance away from boundary truncation. On a periodic embedded grid, this becomes block-circulant after embedding, so matrix-vector products can be evaluated by FFT.

## Design choice

Use the real-space kernel as the primary discretization, not the Fourier-symbol representation. The reason is strategic:

- it stays closer to the paper's DtN construction
- it makes the singular treatment explicit rather than hidden in spectral machinery
- it avoids making periodic-box and Gibbs behavior the first debugging problem

The FFT is used as an accelerator for convolution-like application, not as the definition of the operator.

## Discrete representation

Let the boundary field live on a uniform grid `(x_i, y_j)` with spacings `hx`, `hy`. The operator will be represented by a translation-invariant stencil/kernel:

```text
(N_h f)_{ij} = Σ_{m,n} W_{m,n} (f_{ij} - f_{i-m,j-n}),
```

where `W_{m,n}` approximates

```text
(1 / 2π) ∫_{cell(m,n)} 1 / |r|^3 dA.
```

The kernel is block-Toeplitz with Toeplitz blocks. Only one reference kernel tile needs to be stored.

## Main numerical issue: singular and near-singular treatment

This is the core design problem.

The continuous operator is well-defined because the difference `f(r) - f(s)` cancels the singularity in principal-value sense. A naive discrete kernel will generally destroy that cancellation and produce poor accuracy or instability.

The implementation therefore needs an explicit local treatment.

### Required properties

The discrete operator should satisfy:

- constants map to zero
- symmetry is preserved
- the discrete kernel is even in both coordinates
- the self-term is chosen so the row sum is zero
- near-neighbor weights are computed with higher care than far-field weights

### Proposed strategy

Split the kernel construction into near field and far field.

#### Far field

For offsets with `|(m,n)|` larger than a small cutoff, use cell-center or cell-averaged quadrature:

```text
W_{m,n} ≈ (hx hy) / (2π |r_{m,n}|^3).
```

This is cheap and adequate away from the singularity.

To preserve second-order convergence, the far-field rule should be at least second-order accurate. Plain cell-center sampling is acceptable only if convergence tests confirm that the near-field error dominates and the global operator still converges quadratically. Otherwise, the far field should also use a higher-order cell average.

#### Near field

For a small local neighborhood around the origin, compute weights by a better quadrature over each cell, or by singularity-subtracted integration. The local neighborhood should be symmetric and include at least first and second neighbors.

This near-field construction is the decisive part of achieving `O(dx^2)` convergence. A first-order local treatment is not acceptable even if the FFT application is exact.

The self-weight is not approximated directly from the singular kernel. Instead it is defined by the zero-row-sum condition:

```text
W_{0,0} = - Σ_{(m,n) ≠ (0,0)} W_{m,n}.
```

This enforces `N_h 1 = 0` exactly for the translation-invariant infinite-grid kernel. On a finite zero-extended domain, only the periodic oracle preserves the global constant nullspace exactly; the corrected nonperiodic operator is intended for fields that decay near the box boundary.

However, zero-row-sum alone does not guarantee second-order accuracy. The near-field weights must also reproduce the correct action of the operator on smooth low-order test fields up to the desired truncation order.

### Chosen near-field strategy

The preferred first implementation is a corrected near-field stencil obtained from local polynomial exactness.

The plan is:

- use the real-space singular kernel to define the far-field weights
- choose a small symmetric near-field patch around the origin
- compute the origin-cell contribution analytically from the quadratic Taylor term, which yields a local Laplacian correction
- calibrate the nearest axis-neighbor weights so that the discrete operator is exact on the radial quadratic polynomial over the full `3 x 3` local patch
- keep diagonal and farther weights kernel-based

This is preferred over pure high-order cell integration because it gives direct control over the singular cancellation while preserving the real-space nonperiodic formulation.

### Rollback alternatives

If the polynomial-exact near field is unstable, overfit, or difficult to tune, the alternatives are:

#### Alternative A: high-order near-cell quadrature

Compute the near-field cell weights directly with higher-order quadrature or singularity-subtracted cell integration, while keeping the far field and FFT application unchanged.

Use this if:

- the fitted local weights are too sensitive to calibration choices
- the linear system for correction weights is poorly conditioned
- the quadrature route turns out to be simpler than expected

#### Alternative B: larger near-field patch with fewer matching constraints

Keep the correction idea but increase the local patch and reduce the number of exact matching conditions.

Use this if:

- second-order convergence is nearly achieved but not robust
- local anisotropy remains visible
- the smallest patch is too restrictive

#### Alternative C: spectral DtN as validation oracle or fallback

Keep the Toeplitz/kernel implementation as the main route, but add the Fourier-symbol DtN

```text
\widehat{Nf}(k_x,k_y) = sqrt(k_x^2 + k_y^2) \hat{f}(k_x,k_y)
```

as either:

- a validation oracle for smooth fields
- or the primary operator if the kernel discretization fails to achieve stable second-order accuracy

## FFT-accelerated application

Once the reference kernel is built, applying the operator is a discrete convolution-like operation. The implementation should use circulant embedding and FFTs:

1. embed the kernel in a padded circulant array
2. FFT the kernel once and cache it
3. FFT the input field
4. multiply in Fourier space
5. inverse FFT and crop to the physical grid

This yields `O(N log N)` application cost instead of `O(N^2)`.

## Proposed Julia structure

### Core type

```text
struct ToeplitzDtN2D
    nx::Int
    ny::Int
    hx::Float64
    hy::Float64
    kernel::Matrix{Float64}
    kernel_fft::Matrix{ComplexF64}
    padding_shape::Tuple{Int,Int}
end
```

### Primary API

```text
build_toeplitz_dtn(nx, ny, hx, hy; near_radius=...)
apply!(y, op, x)
apply(op, x)
```

### Internal helpers

```text
build_kernel(...)
compute_far_weight(...)
compute_near_weight(...)
embed_kernel_circulant(...)
```

## Verification plan

The operator must be validated before any coupling to the impact problem.

### Test 1: constants

Input: constant field.

Expected:

```text
N_h 1 = 0
```

up to machine precision.

### Test 2: symmetry

For symmetric inputs, the output should preserve expected even/odd symmetry.

### Test 3: translation invariance

Applying the operator to a shifted field should agree with shifting the output, modulo boundary/padding conventions.

### Test 4: Fourier-mode check

Even though the implementation is real-space, Fourier modes are still the cleanest reference test for the periodic oracle. For

```text
f(x,y) = exp(i(k_x x + k_y y)),
```

the exact DtN value is

```text
Nf = |k| f.
```

This gives a direct error metric for the periodic reference implementation.

### Test 5: grid refinement

Measure error decay under refinement for smooth test fields. This is the main evidence that the singular treatment is correct.

### Test 6: explicit second-order target

Run a refinement study on square grids with `dx = dy` and verify

```text
||N_h f - Nf|| ≤ C dx^2
```

for smooth test functions over the pre-asymptotic range before finite-box truncation dominates. For the nonperiodic operator, fixed-domain self-convergence is often a cleaner test than direct comparison to the infinite-domain exact DtN because the latter includes truncation error from zero extension outside the computational box.

## Out-of-scope for first version

The first version should not include:

- contact constraints
- pressure solve
- body dynamics
- nonuniform grids
- boundary-layer viscous correction

Those depend on first having a trustworthy `N*x`.

## Fallback if the kernel route underperforms

If the local singular treatment proves too fragile or too inaccurate, the fallback is not to abandon the project, but to add a reference Fourier-symbol implementation:

```text
\widehat{Nf}(k_x,k_y) = sqrt(k_x^2 + k_y^2) \hat{f}(k_x,k_y),
```

and use it as either:

- a validation oracle
- or the primary operator, with the kernel route retained for comparison

## Recommendation

Proceed with the Toeplitz/kernel DtN as the first Julia target. The singular/local quadrature is the main technical milestone, and it must be designed around a non-negotiable second-order accuracy target for `hx = hy = dx`. The current implementation path is:

- keep a periodic Fourier-symbol DtN as an oracle
- use a zero-extended real-space Toeplitz kernel as the main operator
- correct the singular origin cell analytically
- calibrate the nearest axis-neighbor weights by quadratic exactness on the local `3 x 3` patch

Once that is correct, the rest of the operator application is straightforward and efficient.

# Memo: Fourier-space DtN versus real-space kernel DtN

## Purpose

This note explains the proposed Dirichlet-to-Neumann (DtN) implementation strategy for the fluid half-space problem. The goal is not to defend Fourier methods on taste, but to explain exactly what is gained, what is lost, and how this compares to the real-space singular-kernel viewpoint that is closer to Galeano-Rios et al. (2017).

## 1. What the DtN operator is

Suppose `φ(x, y, z)` is harmonic in the lower half-space,

```text
Δφ = 0,   z < 0,
```

with boundary trace

```text
φ(x, y, 0) = f(x, y).
```

The Dirichlet-to-Neumann map is the operator

```text
Nf = ∂_z φ(x, y, 0),
```

meaning:

- input: boundary values `f`
- output: the normal derivative of the harmonic extension of `f`

So `N` is an operator on boundary functions. In the continuous problem it is not "a matrix" first. It becomes a matrix only after discretization.

## 2. What it means to Fourier transform an operator

This is the key point.

For translation-invariant operators, complex exponentials are eigenfunctions. If

```text
f(x, y) = exp(i k · x),   k = (k_x, k_y),
```

then the harmonic extension into `z < 0` is

```text
φ(x, y, z) = exp(i k · x) exp(|k| z),
```

because `z < 0` and decay at depth requires `exp(|k| z)`.

Therefore

```text
∂_z φ(x, y, 0) = |k| exp(i k · x),
```

so the DtN operator acts on a Fourier mode by multiplication:

```text
N exp(i k · x) = |k| exp(i k · x).
```

That is the entire meaning of saying that the operator has Fourier symbol `|k|`. In Fourier space, the operator is diagonal:

```text
\widehat{Nf}(k) = |k| \hat{f}(k).
```

This is not a different physical model. It is the same half-space DtN operator written in the basis where translation-invariant operators are simplest.

## 3. Relation to the paper's real-space formula

The paper writes the DtN map in real space as the singular integral

```text
Nf(r) = (1 / 2π) PV ∫ (f(r) - f(s)) / |r - s|^3 dA(s).
```

This is the same operator.

So there are two mathematically equivalent representations:

- real space: principal-value singular integral
- Fourier space: multiplier `|k|`

The question is not which operator is correct. Both are. The question is which representation is easier to discretize well for the geometry you care about.

## 4. Why translation invariance matters

On a flat infinite plane, the kernel depends only on `r - s`, not on `r` and `s` separately. That is why:

- in continuous form, the DtN map is a convolution-type operator
- in discrete form on a uniform grid, the matrix depends only on index differences

In 1D that gives Toeplitz structure. In 2D it gives block-Toeplitz with Toeplitz blocks. On a periodic grid it becomes block-circulant, which is exactly diagonalized by the discrete Fourier transform.

This is why Fourier keeps reappearing. It is the natural basis for translation-invariant operators.

## 5. Clarifying the Toeplitz point

You were right to push back on the "explicit matrix" language. If you store only one row because the operator is Toeplitz or block-Toeplitz, that is already a structured operator representation, not a dense explicit matrix build.

The distinction I meant was this:

### Toeplitz-first view

You start from the real-space kernel, discretize it into one reference row, and use FFT-based convolution to apply the operator efficiently.

### Fourier-symbol-first view

You start from the exact spectral symbol `|k|` and apply the operator by:

1. FFT the grid values
2. multiply by `|k|`
3. inverse FFT

Both are matrix-free in practice. Neither requires forming a dense matrix. Neither requires inverting the DtN operator.

The real difference is where the discretization difficulty sits:

- Toeplitz/kernel-first: difficulty is near the singular kernel and the self-term
- Fourier-symbol-first: difficulty shifts to domain truncation, periodicity, and contact-edge resolution

## 6. Why I suggested the Fourier-symbol-first route

For the flat half-space problem, the symbol `|k|` is exact. That gives several real advantages.

### Advantage 1: the singularity disappears analytically

In real space, the kernel is singular at `r = s`, so the discrete self-interaction needs special treatment. In Fourier space, that singularity is already encoded in the multiplier `|k|`; there is no principal-value quadrature to design.

### Advantage 2: verification is cleaner

If the operator is implemented as

```text
\widehat{Nf}(k) = |k| \hat{f}(k),
```

then every Fourier mode is an exact test case. This is extremely useful for debugging.

### Advantage 3: non-axisymmetric generalization is immediate

The paper's radial construction is natural in axisymmetry, but awkward in full 2D boundary coordinates. The spectral operator does not care whether the state is axisymmetric or not.

### Advantage 4: no kernel tabulation is needed

The operator is defined directly by its symbol. There is no need to derive, store, or regularize a reference kernel row unless you specifically want a real-space method.

## 7. What the Fourier-symbol-first route does not magically solve

This route only simplifies the fluid operator. It does not solve:

- the moving contact set
- the pressure solve in the constrained region
- sharp contact-edge regularity
- non-periodic or awkward outer-domain geometry

Those are separate issues.

## 8. The three implementation choices I referred to

The earlier list was too compressed. Here is what it actually meant.

### Choice A: exact periodic Fourier DtN

Discretize the boundary on a rectangular periodic grid and define

```text
\widehat{Nf}(k_x, k_y) = sqrt(k_x^2 + k_y^2) \hat{f}(k_x, k_y).
```

Implementation:

1. store `f` on a 2D grid
2. FFT to get `\hat{f}`
3. multiply by the spectral weights
4. inverse FFT

This is the cleanest implementation of the exact half-space DtN on a periodic box.

What can go wrong:

- periodic copies of the solution may interact
- the contact edge can create Gibbs oscillations
- the physical problem is not naturally periodic, so padding or large domains may be needed

### Choice B: real-space singular-kernel discretization

Discretize the principal-value integral directly:

```text
Nf(r_i) ≈ Σ_j K_{i-j} (f_i - f_j),
```

with a specially designed treatment of the near-singular region and self-term.

This is closest in spirit to the paper.

What can go wrong:

- the self-term is delicate
- near-field quadrature accuracy matters a lot
- in 2D the discrete operator is more involved than in the radial case

The upside is that it can be made less dependent on periodic assumptions.

### Choice C: hybrid near-field / far-field method

Split the operator into:

- near field: handled in real space with singularity subtraction or exact local quadrature
- far field: handled spectrally or by FFT convolution

This is more complicated, but it can combine the strengths of both approaches if the pure spectral method struggles at the contact edge.

## 9. Why the singularity issue is different in the two views

This was the "main technical issue" I referred to.

In real space, the continuous formula has a singular kernel:

```text
(f(r) - f(s)) / |r - s|^3.
```

The numerator cancels the singularity in the principal-value sense, but the discrete scheme must preserve that cancellation. That is numerically delicate.

In Fourier space, the same singular structure is already compressed into the multiplier `|k|`. So you no longer have to discretize a singular integral directly. The price is that you accept a periodic-box representation of the boundary operator.

## 10. When deviating from the paper is truly justified

Deviation is justified if the target problem has changed.

If the target is still:

- axisymmetric
- radial mesh
- contact radius search

then the paper's style of DtN is natural.

If the target is:

- non-axisymmetric contact
- 2D boundary field
- efficient repeated operator application

then the Fourier-symbol-first route has a strong and concrete advantage. It is not just "more convoluted." It is the representation naturally matched to a flat 2D boundary and a translation-invariant operator.

## 11. Recommendation

For this project, the most sensible path is:

1. implement the exact periodic Fourier DtN first as the reference operator
2. verify it on individual Fourier modes
3. only build a real-space Toeplitz/kernel version if one of these becomes a real blocker:
   - periodic-image pollution
   - severe edge artifacts at the contact set
   - need for a non-periodic outer-domain treatment

This order is not ideological. It is practical. The spectral route gives the fastest path to a correct non-axisymmetric DtN operator with the least ambiguity in the operator itself.

## 12. Backup ideas if the first route struggles

### A. Padded Fourier DtN

Keep the spectral operator but enlarge the computational box to reduce wrap-around effects.

### B. Kernel-based Toeplitz DtN

Use the real-space singular-integral discretization, store one reference stencil/row, and apply it by FFT-based convolution.

### C. Hybrid near/far split

Treat local singular interactions in real space and use FFT structure for the far field.

### D. Boundary-layer correction

If the true issue turns out to be viscous stress near the contact edge rather than the DtN operator itself, keep the Fourier DtN for the harmonic part and add a boundary-layer viscous correction separately.

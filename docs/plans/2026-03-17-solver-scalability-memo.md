# Solver Scalability Memo

## Context

The current solver scaffold is correct as a testbed, but it is not scalable. The bottleneck is the Dirichlet-to-Neumann block: if it is assembled densely, then the full implicit Euler system

```text
A(dt) x = b
```

will become too large to store or factor for realistic two-dimensional grids.

The right question is therefore not simply whether to use a direct or iterative solver. The real question is which parts of `A(dt)` should ever be assembled explicitly, and which parts should remain operator applications.

## Structure of the linear system

For the current prescribed-contact formulation, the one-step system has three qualitatively different pieces:

1. finite-difference bath blocks
2. contact-pressure and rigid-plane coupling blocks
3. the DtN block

The finite-difference blocks are sparse. The contact and rigid-body blocks are sparse or very low-dimensional. The DtN block is the only nonlocal piece.

This suggests the decomposition

```text
A(dt) = S(dt) + C(dt),
```

where

- `S(dt)` is sparse and local
- `C(dt)` is the DtN contribution, which should be applied matrix-free

## Candidate approaches

### 1. Full matrix-free Krylov solve

The first option is to never assemble the full matrix at all. Instead, define the action

```text
v ↦ A(dt)v
```

by blocks:

- sparse matrix-vector products for Laplacian and time-stepping blocks
- FFT-based DtN application for the nonlocal block
- tiny dense operations for rigid-plane and body coupling

Then solve with GMRES, or MINRES if the final algebraic form is made symmetric enough.

#### Advantages

- no dense DtN storage
- no dense factorization
- naturally compatible with FFT-based Toeplitz DtN application
- easiest path to large grids

#### Disadvantages

- requires a good preconditioner
- convergence may degrade as `dt`, stiffness, or contact complexity grow
- implementation is more algebraic than the current assembled testbed

#### Verdict

This is the most flexible long-term path.

### 2. Schur complement reduction

The second option is to eliminate the small rigid-body/contact unknowns, or eliminate selected bath variables, to form a reduced system. Since the body block is very small, one can expect some benefit from eliminating it analytically.

A typical outcome would be a reduced bath system of the form

```text
K(dt) u = f,
```

where `K(dt)` consists of sparse finite-difference terms, DtN action, and a low-rank correction induced by contact/body coupling.

#### Advantages

- exposes the true large-scale bath operator more clearly
- keeps body coupling exact
- may improve preconditioning design

#### Disadvantages

- algebra becomes more involved
- easy to accidentally densify the reduced operator if done naively
- not obviously simpler than matrix-free GMRES on the full system

#### Verdict

Useful as an analysis tool and possibly as an implementation refinement later, but not the best first scalable solver path.

### 3. Sparse-direct preconditioner plus iterative outer solve

The third option is a hybrid:

1. build a sparse approximation `P(dt)` to `A(dt)` that drops or simplifies the DtN part
2. factor `P(dt)` once
3. use that factorization as a preconditioner inside GMRES for the true operator

The approximation could be:

- the full sparse finite-difference/time-stepping/body block, with DtN removed
- or the same sparse block plus a local approximation to DtN

#### Advantages

- very practical
- strong reuse when `dt` and contact mask are fixed
- keeps the true DtN action matrix-free
- usually much easier to stabilize than an unpreconditioned Krylov solve

#### Disadvantages

- depends on the quality of the sparse approximation
- still requires matrix-free operator infrastructure
- preconditioner quality may worsen if the DtN term dominates

#### Verdict

This is probably the best near-term production route.

## Exploiting matrix structure

The nonlocal DtN operator should not be viewed as a generic dense matrix. On a uniform grid it is structured:

- approximately Toeplitz in the nonperiodic setting
- exactly diagonalizable in Fourier space in the periodic oracle
- naturally applied by FFT after circulant embedding

So the right optimization is not “compress the dense matrix.” The right optimization is:

```text
apply DtN(v) by FFT
```

inside a larger matrix-free block solver.

The remaining pieces of the step operator have exploitable structure as well:

- Laplacian blocks are sparse
- contact mask insert/extract maps are sparse
- rigid-body coupling is low-rank
- the body mass matrix is `3 x 3`

This strongly favors a block-operator implementation over dense assembly.

## Repeated solves for fixed `dt`

If `dt` is fixed and the contact mask is fixed, then there is substantial reuse available even without forming `A^{-1}`.

### Reuse 1: preconditioner

Construct `P(dt)` once and reuse it at every step.

This is the biggest gain.

### Reuse 2: sparse factorization of the preconditioner

If `P(dt)` is sparse, compute its LU or LDLT factorization once and reuse triangular solves inside every Krylov iteration at every time step.

### Reuse 3: FFT plans and DtN embedding

The FFT plan, padded grid shape, and DtN kernel embedding can all be built once and reused.

### Reuse 4: fixed block maps

The free/contact index maps, pressure injection maps, and rigid-plane basis are all fixed if the contact mask is fixed.

Thus, for fixed `dt` and fixed contact mask, the only thing changing step to step is the right-hand side and later, if nonlinearities are added, some explicit source terms.

## Recommendation

The dense assembled solver should remain as the correctness reference only.

The production path should be:

1. refactor the solver into block operator applications
2. keep DtN matrix-free and FFT-based
3. use GMRES on the full step system
4. precondition with a sparse approximation built from the finite-difference and body blocks
5. reuse the preconditioner factorization for all solves with fixed `dt` and fixed contact mask

In short:

- do not assemble dense DtN
- do not store dense `A`
- do not attempt to store `A^{-1}`
- do exploit FFT structure for DtN
- do exploit sparse-direct reuse for the preconditioner

## Practical next step

The next implementation step should be to replace the assembled `assemble_step_system` path by two parallel interfaces:

1. `assemble_step_system` for correctness and testing
2. `apply_step_operator!` plus a preconditioner object for scalable solves

That preserves the current testbed while opening the path to large-scale runs.

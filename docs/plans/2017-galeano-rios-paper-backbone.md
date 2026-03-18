# Galeano-Rios et al. (2017): implementation backbone

Reference: C. A. Galeano-Rios, P. A. Milewski, and J.-M. Vanden-Broeck, "Non-wetting impact of a sphere onto a bath and its application to bouncing droplets," J. Fluid Mech. 826, 97-127 (2017), doi:10.1017/jfm.2017.424.

## What the paper really solves

After the weakly viscous reduction, the operative system is not the full Navier-Stokes problem but the boundary/contact system

```text
η_t = Nφ + (2/Re) Δ_H η                                             (2.29a)

φ_t = -(1/Fr) η + (1/We) Δ_H η + (2/Re) Δ_H φ - p_s
      + (1/We) (κ - Δ_H) η                                          (2.29b)

h_tt = -(1/Fr) - D h_t + (1/M) ∫_{r≤r_c} p_s dA                     (2.29c)
```

with mixed contact conditions

```text
η < h + z_s      for r_c < r < R_o                                  (2.30a)
p_s = 0          for r > r_c                                        (2.30b)
η = h + z_s      for r ≤ r_c                                        (2.30c)
∂_r η(r_c) = ∂_r z_s(r_c)                                           (2.30d)
```

All spatial functions are radial, `r_c` is unknown, and the fluid is reduced to the boundary through the Dirichlet-to-Neumann map `N`.

## Core modeling assumptions

The fluid is linearized around a flat interface. Gravity, capillarity, and weak viscous damping are retained. The solid is smooth, convex, axisymmetric, and perfectly hydrophobic. The air layer is not solved; it is only assumed to transmit pressure. The impact stays in the non-penetrating regime, so the whole method is aimed at early contact and bounce, not cavity formation or violent entry.

## Three design choices worth separating

### 1. Elimination of vorticity

This is the move that makes the boundary-only formulation possible. Once the effective bulk field is potential, `φ` is harmonic in the half-space and the DtN map `Nφ = ∂_z φ|_{z=0}` is well-defined as a boundary operator. If one keeps a genuine vortical component in the bulk, the state is no longer represented by a single harmonic scalar potential, so the reduction to a scalar DtN map is lost or at least no longer clean.

Conclusion: the paper's boundary formulation and the removal of vorticity are structurally linked.

### 2. Axisymmetry

This is the least fundamental assumption and the one most worth removing. It simplifies:

- the geometry of the contact set
- the Laplacian
- the DtN implementation
- the search for `r_c`

But none of those simplifications are essential to the physical idea. If the real goal is non-axisymmetric impacts or later horizontal dynamics, this is the assumption to relax first.

### 3. DtN implementation

The paper implements DtN in radial form through a singular integral and a precomputed matrix. That is reasonable in the axisymmetric setting, but it does not scale naturally to a general 2D boundary.

For a flat interface over infinite depth, the non-axisymmetric DtN operator is especially simple in Fourier space:

```text
φ̂_z(k) = |k| φ̂(k),    k = (k_x, k_y)
```

so

```text
N = |D|
```

on the boundary plane.

This is the main reason a Fourier DtN is attractive.

## Why Fourier DtN is probably better for the next implementation

The paper's non-Fourier DtN has one major advantage: it is geometrically explicit and matches the radial contact geometry directly. But for the general 2D boundary problem, Fourier DtN has stronger engineering advantages:

- It is the exact DtN operator for the linear half-space problem on a flat boundary.
- It avoids assembling and storing a dense singular-integral matrix over a 2D surface mesh.
- It replaces principal-value quadrature with FFTs and a diagonal multiplier `|k|`.
- It generalizes to non-axisymmetric states immediately.
- It is easier to verify: the operator symbol is known analytically.
- It decouples the fluid operator from the contact geometry, which is useful because the hard part is already the moving constrained region.

The tradeoff is equally clear:

- It assumes a flat reference boundary and is most natural on periodic or padded rectangular domains.
- It does not by itself solve the contact problem; it only gives a cleaner fluid operator.
- If the domain geometry becomes non-flat or finite-depth/bounded in a complicated way, the spectral shortcut weakens.

So the comparison is not "paper DtN bad, Fourier DtN good." It is:

- radial singular-integral DtN is natural for the paper's axisymmetric solver
- Fourier DtN is natural for a non-axisymmetric solver on a flat plane

For the next step, the latter matches the real target better.

## Important caveat on the vorticity objection

It is plausible that dropping the vortical contribution weakens the viscous stress description near the contact edge. But the curvature/stress irregularity at the interface between contact and free surface is not caused only by dropping vorticity. The mixed formulation itself is sharp:

- inside contact: `η = h + z_s`
- outside contact: `p_s = 0` and the free surface evolves dynamically
- at the edge: only height and slope are matched explicitly

So even with a better viscous model, a sharp transition can still generate a curvature mismatch or strong localized gradients. In other words, keeping vorticity may help, but it is not guaranteed to remove the edge pathology by itself.

## Recommended implementation direction

The most defensible path is:

1. Keep the boundary-only weakly viscous formulation first.
2. Remove axisymmetry.
3. Implement the fluid operator with Fourier DtN on a 2D boundary grid.
4. Treat the contact region as a moving constrained subset of the boundary.
5. Revisit viscous corrections only after the non-axisymmetric contact solver is working.

This keeps the main innovation focused on the part that matters most: non-axisymmetric contact and wave generation.

## Alternative ideas if the main route fails

### A. Fourier DtN with non-axisymmetric contact

This is the preferred route. Use a rectangular boundary grid, FFT-based `N = |D|`, and represent contact as an active set where `η = h + z_s` and `p_s` acts as a Lagrange multiplier or constrained pressure field.

Failure modes:

- domain truncation or periodic-image contamination
- difficulty enforcing a sharp moving contact set on a spectral grid
- Gibbs-type artifacts near the contact edge

### B. Boundary-layer correction on top of the potential/DtN model

Instead of restoring full bulk vorticity, add a viscous correction that lives on the boundary or in a thin-layer approximation. This is more compatible with a boundary solver than a full vortical bulk model.

Use this if:

- the pure weak-dissipation closure produces poor forces near lift-off/contact-edge events
- curvature/stress behavior near the contact line looks systematically wrong
- you want a viscosity upgrade without abandoning DtN

### C. Full bulk solve with vorticity retained

This is the most expensive and least aligned with the original reduced model. It may eventually be necessary, but it should be treated as a different model class, not as a small extension of Galeano-Rios.

Use this only if:

- boundary-only corrections fail qualitatively
- edge stresses are dominated by effects the reduced model cannot represent
- the main scientific question becomes viscous near-contact structure rather than efficient non-axisymmetric impact prediction

## What to inspect in the neighboring `km-*` repos

When comparing with `/Users/eaguerov/Documents/Github/km-disk-vibrating`, `/Users/eaguerov/Documents/Github/km-dropplet-onto-bath`, and `/Users/eaguerov/Documents/Github/km-dropplet-solidsubstrate-v3`, the right checklist is:

- Do they keep a boundary-only fluid representation or a bulk one?
- Do they assume axisymmetry or encode general 2D boundary geometry?
- How is contact enforced: explicit radius search, active set, penalty, or pressure solve?
- Is the fluid operator radial, matrix-based DtN, Fourier/spectral DtN, or something else?
- Where is viscosity handled: pure damping term, boundary-layer correction, or full bulk diffusion?

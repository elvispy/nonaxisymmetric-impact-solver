# Project File Documentation Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Create durable project documentation that explains the purpose and relationships of the current files in `vertical-dynamics`.

**Architecture:** The documentation should be organized around the code that actually exists: package entrypoint, DtN subsystem, solver subsystem, tests, and research/planning artifacts. The first output should be a top-level project map, followed by focused subsystem docs that stay aligned with the source tree and test suite.

**Tech Stack:** Markdown documentation in `docs/`, Julia package layout under `src/` and `test/`, existing planning notes in `docs/plans/`.

---

## File Structure

The documentation work should treat the repository as five distinct units.

- `Project.toml`, `Manifest.toml`
  Responsibility: package metadata and dependency lockfile.

- `src/VerticalDynamics.jl`
  Responsibility: package entrypoint and public exports.

- `src/dtn/`
  Responsibility: Dirichlet-to-Neumann operators, including periodic reference and corrected nonperiodic operator application.

- `src/solver/`
  Responsibility: prescribed-contact rigid-plane solver, domain construction, step assembly, matrix-free operator path, and Krylov solve.

- `test/`
  Responsibility: verification of DtN operators and solver behavior.

- `docs/plans/`
  Responsibility: design history, research notes, and implementation memos.

- `.requirements/`
  Responsibility: per-task requirements snapshots used during development.

## Documentation Targets

The documentation should be created in layers.

1. A top-level project map
   Purpose: explain what the repository is, what is implemented, and where to look first.

2. A source tree guide
   Purpose: document the responsibility of each current file under `src/`.

3. A test guide
   Purpose: document what each test file verifies and which implementation files it protects.

4. A plans/notes guide
   Purpose: distinguish active implementation docs from historical reasoning notes.

## Chunk 1: Top-Level Project Map

### Task 1: Create the main repository documentation page

**Files:**
- Create: `docs/project-map.md`
- Modify: `README.md` or create `README.md` if absent

- [ ] **Step 1: Write the failing documentation test for presence and coverage**

Create a checklist in the implementation notes requiring:

```text
docs/project-map.md must describe:
- package purpose
- main source directories
- test directories
- planning directories
```

- [ ] **Step 2: Verify the current repo lacks the documentation entrypoint**

Run: `test -f README.md && echo exists || echo missing`
Expected: likely `missing` or incomplete top-level documentation.

- [ ] **Step 3: Write the project map**

Document:

- what `vertical-dynamics` currently implements
- the distinction between DtN infrastructure and solver infrastructure
- where correctness-reference code differs from scalable code paths

- [ ] **Step 4: Link the project map from the repository entrypoint**

If `README.md` is absent, create a short one. If it exists, add a pointer to `docs/project-map.md`.

- [ ] **Step 5: Review for accuracy against the source tree**

Run: `find src test docs/plans .requirements -maxdepth 3 -type f | sort`
Expected: all referenced files exist and are grouped correctly in the new docs.

## Chunk 2: Source Tree Guide

### Task 2: Document every file under `src/`

**Files:**
- Create: `docs/source-tree.md`

- [ ] **Step 1: Write a file inventory from the current source tree**

Use this current structure:

```text
src/VerticalDynamics.jl
src/dtn/apply.jl
src/dtn/toeplitz_operator.jl
src/solver/assembly.jl
src/solver/domain.jl
src/solver/krylov.jl
src/solver/operator.jl
src/solver/time_marching.jl
src/solver/types.jl
```

- [ ] **Step 2: For each file, document one clear responsibility**

Required sections:

- file path
- purpose
- key exported or internal functions/types
- direct dependencies within `src/`

- [ ] **Step 3: Add subsystem-level overviews**

Document:

- `dtn/` as operator infrastructure
- `solver/` as timestepper and linear solve infrastructure

- [ ] **Step 4: Cross-check against `src/VerticalDynamics.jl` exports**

Run: `sed -n '1,200p' src/VerticalDynamics.jl`
Expected: every exported symbol belongs to a documented source file.

## Chunk 3: Test Guide

### Task 3: Document the test suite by intent

**Files:**
- Create: `docs/test-guide.md`

- [ ] **Step 1: Inventory the test files**

Use:

```text
test/test_api.jl
test/test_convergence.jl
test/test_fft_apply.jl
test/test_fourier_modes.jl
test/test_impulse_response.jl
test/test_kernel_structure.jl
test/test_nullspace.jl
test/test_solver_api.jl
test/test_solver_assembly.jl
test/test_solver_domain.jl
test/test_solver_operator.jl
test/test_solver_step.jl
```

- [ ] **Step 2: For each test file, state what behavior it protects**

Required columns or bullets:

- test file
- implementation files covered
- behavior verified
- whether it is a DtN test or solver test

- [ ] **Step 3: Document the canonical verification command**

Include:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

and explain what a green run means at the current stage of the project.

## Chunk 4: Plans And Requirements Guide

### Task 4: Document the role of planning artifacts

**Files:**
- Create: `docs/plans/README.md`

- [ ] **Step 1: Inventory current planning and requirements artifacts**

Current known files:

```text
docs/plans/2017-galeano-rios-paper-backbone.md
docs/plans/2026-03-17-julia-toeplitz-dtn-implementation-plan.md
docs/plans/2026-03-17-prescribed-contact-plane-solver.md
docs/plans/2026-03-17-solver-scalability-memo.md
docs/plans/dtn-fourier-memo.md
docs/plans/julia-toeplitz-dtn-design-brief.md
.requirements/20260317T202333Z_corrected_kernel_dtn/REQUIREMENTS.md
.requirements/20260317T231548Z_prescribed_contact_plane_solver/REQUIREMENTS.md
```

- [ ] **Step 2: Group them by purpose**

Required groups:

- literature/model notes
- implementation plans
- scalability/architecture memos
- per-task requirements snapshots

- [ ] **Step 3: Add guidance on which documents are authoritative**

Clarify:

- source code is authoritative for implementation
- tests are authoritative for verified behavior
- plans and memos capture rationale and intended direction

## Chunk 5: Final Integration

### Task 5: Stitch the docs into a coherent navigation path

**Files:**
- Modify: `README.md`
- Modify: `docs/project-map.md`
- Modify: `docs/source-tree.md`
- Modify: `docs/test-guide.md`
- Modify: `docs/plans/README.md`

- [ ] **Step 1: Add cross-links between the new docs**

Required links:

- `README.md` → `docs/project-map.md`
- `docs/project-map.md` → `docs/source-tree.md`
- `docs/project-map.md` → `docs/test-guide.md`
- `docs/project-map.md` → `docs/plans/README.md`

- [ ] **Step 2: Run a final existence check**

Run: `find docs -maxdepth 3 -type f | sort`
Expected: all new documentation files exist and are linked from the top-level map.

- [ ] **Step 3: Run the Julia test suite to ensure documentation work did not disturb code paths**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: PASS

- [ ] **Step 4: Commit**

```bash
git add README.md docs
git commit -m "docs: add project file documentation"
```

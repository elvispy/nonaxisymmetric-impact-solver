# Prescribed-Contact Plane Solver Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a first non-axisymmetric finite-difference impact solver in Julia with a prescribed contact mask, rigid-plane contact kinematics, unknown contact pressure, and implicit Euler time stepping.

**Architecture:** The solver mirrors the structure of `km-disk-vibrating`: a domain builder prepares spatial operators, a one-step assembly builds a square implicit linear system, and a top-level `solve_motion` loop records the solution history. The bath is fully 2D on a Cartesian grid, while the contacting body is reduced to rigid-plane generalized coordinates `q = (a, b, c)` with a generalized mass matrix induced by a prescribed grid density field.

**Tech Stack:** Julia 1.12, `FFTW`, `LinearAlgebra`, Julia `Test` stdlib.

---

## Chunk 1: Solver Scaffolding

### Task 1: Add solver module files and exports

**Files:**
- Create: `src/solver/types.jl`
- Create: `src/solver/domain.jl`
- Create: `src/solver/assembly.jl`
- Create: `src/solver/time_marching.jl`
- Modify: `src/VerticalDynamics.jl`
- Create: `test/test_solver_api.jl`

- [ ] **Step 1: Write the failing test**
- [ ] **Step 2: Run the solver test to verify the API is missing**
- [ ] **Step 3: Add the minimal solver types and exports**
- [ ] **Step 4: Run the test to verify construction succeeds**

## Chunk 2: Domain And Rigid-Body Operators

### Task 2: Build Cartesian domain operators

**Files:**
- Modify: `src/solver/domain.jl`
- Create: `test/test_solver_domain.jl`

- [ ] **Step 1: Add tests for Laplacian size, mask bookkeeping, and rigid-body mass moments**
- [ ] **Step 2: Implement the grid coordinates, 5-point Laplacian, and generalized mass matrix**
- [ ] **Step 3: Run the domain tests until they pass**

## Chunk 3: Implicit Step Assembly

### Task 3: Assemble the square linear system

**Files:**
- Modify: `src/solver/assembly.jl`
- Create: `test/test_solver_assembly.jl`

- [ ] **Step 1: Add a failing test for expected matrix dimensions**
- [ ] **Step 2: Implement dense DtN matrix assembly and the block implicit system**
- [ ] **Step 3: Run the assembly test and verify the system is square with size `2N + 6`**

## Chunk 4: Time Marching

### Task 4: Advance one step and solve multiple steps

**Files:**
- Modify: `src/solver/time_marching.jl`
- Create: `test/test_solver_step.jl`
- Modify: `test/runtests.jl`

- [ ] **Step 1: Add failing tests for zero-solution invariance and rigid-plane enforcement on the contact mask**
- [ ] **Step 2: Implement `advance_one_step` and `solve_motion`**
- [ ] **Step 3: Run the full Julia test suite and fix remaining issues**

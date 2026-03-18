# Case Video Script Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a small Julia script that runs one representative fixed-contact case and writes a top-view heatmap MP4 of the free-surface elevation.

**Architecture:** Keep the solver untouched and add a focused visualization helper under `src/` so the script stays thin. Use `Plots.jl` to render one heatmap per saved timestep with a fixed global color range, then call local `ffmpeg` to assemble those frames into a video under `outputs/`.

**Tech Stack:** Julia, `Plots.jl`, local `ffmpeg`

---

## File Map

- Create: `src/visualization.jl`
- Modify: `src/VerticalDynamics.jl`
- Modify: `Project.toml`
- Modify: `Manifest.toml`
- Create: `scripts/run_case_and_make_video.jl`
- Create: `test/test_visualization.jl`
- Modify: `test/runtests.jl`

## Chunk 1: Add Visualization Helpers

### Task 1: Add the plotting dependency and wire a visualization source file

**Files:**
- Modify: `Project.toml`
- Modify: `Manifest.toml`
- Create: `src/visualization.jl`
- Modify: `src/VerticalDynamics.jl`

- [ ] **Step 1: Add `Plots` to the Julia project**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.add("Plots")'
```

Expected: `Project.toml` and `Manifest.toml` gain `Plots` and its resolved dependencies.

- [ ] **Step 2: Include the visualization file in the package entrypoint**

Add:

```julia
include("visualization.jl")
export compute_eta_color_limits
export render_eta_video
```

to `src/VerticalDynamics.jl`.

- [ ] **Step 3: Create the visualization helper file with docstring-covered function shells**

Add:

```julia
"""
    compute_eta_color_limits(history)

Return a consistent `(zmin, zmax)` pair for heatmap rendering across all saved states.
"""

"""
    render_eta_video(history, x, y, output_path; fps=20, title_prefix="eta", frames_dir=nothing)

Render top-view heatmap frames of `η` and assemble them into an MP4 with `ffmpeg`.
"""
```

- [ ] **Step 4: Commit**

```bash
git add Project.toml Manifest.toml src/VerticalDynamics.jl src/visualization.jl
git commit -m "Add visualization helpers for eta video rendering"
```

## Chunk 2: Test a Pure Visualization Helper

### Task 2: Add and pass a focused visualization unit test

**Files:**
- Create: `test/test_visualization.jl`
- Modify: `test/runtests.jl`

- [ ] **Step 1: Write the failing test**

Add a test that builds a short fake history of `SolverState` values with known `η` extrema and checks:

```julia
@test compute_eta_color_limits(history) == (-2.0, 3.0)
```

- [ ] **Step 2: Run the focused test and verify failure**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test(test_args=["test_visualization"])'
```

Expected: fail because the helper is missing or incomplete.

- [ ] **Step 3: Implement the minimal helper**

Implement `compute_eta_color_limits(history)` in `src/visualization.jl`.

- [ ] **Step 4: Re-run the focused test**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test(test_args=["test_visualization"])'
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add test/test_visualization.jl test/runtests.jl src/visualization.jl
git commit -m "Test eta color-limit helper"
```

## Chunk 3: Implement Video Rendering

### Task 3: Implement frame rendering and ffmpeg assembly

**Files:**
- Modify: `src/visualization.jl`

- [ ] **Step 1: Write the minimal rendering implementation**

Implement:
- output directory creation
- global color-limit reuse
- one PNG frame per state using `heatmap`
- `ffmpeg` invocation to assemble `frame_%05d.png` into MP4
- explicit error if `ffmpeg` is missing

- [ ] **Step 2: Smoke-check the helper on a tiny synthetic history if useful**

Run a small inline Julia command that writes a short video to `/tmp`.

- [ ] **Step 3: Commit**

```bash
git add src/visualization.jl
git commit -m "Implement eta heatmap video rendering"
```

## Chunk 4: Add the Run Script

### Task 4: Add a script that runs one representative case and renders the video

**Files:**
- Create: `scripts/run_case_and_make_video.jl`

- [ ] **Step 1: Create the script**

The script should:
- activate the local project
- define a small fixed-contact case
- call `solve_motion`
- build output paths under `outputs/`
- call `render_eta_video`
- print the final video path

- [ ] **Step 2: Run the script as a smoke test**

Run:

```bash
julia --project=. scripts/run_case_and_make_video.jl
```

Expected: a video appears under `outputs/`.

- [ ] **Step 3: Commit**

```bash
git add scripts/run_case_and_make_video.jl
git commit -m "Add demo script for eta video output"
```

## Chunk 5: Final Verification

### Task 5: Run the full suite and summarize usage

**Files:**
- Modify: `README.md`

- [ ] **Step 1: Add a short usage note for the new script**

Add a brief command example to `README.md`.

- [ ] **Step 2: Run the full test suite**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS.

- [ ] **Step 3: Re-run the script if needed and confirm the MP4 path**

- [ ] **Step 4: Commit**

```bash
git add README.md
git commit -m "Document case video script"
```

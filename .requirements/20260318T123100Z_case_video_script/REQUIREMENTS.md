# As Is

The repository contains a working prescribed-contact solver and a matrix-free DtN path, but there is no script for running a representative case and no visualization pipeline for turning solver output into a video. The Julia project also does not yet depend on a plotting library.

# To Be

The repository should provide a small runnable script that builds one representative fixed-contact case, runs the solver, renders a top-view heatmap of the free-surface elevation `η(x,y,t)` for each saved timestep, and assembles those frames into an MP4 video. The implementation should stay simple: a plotting dependency in Julia for frame generation, and local `ffmpeg` for video assembly.

# Requirements

1. The project shall expose a reusable Julia helper that renders a sequence of solver states to top-view heatmap frames and assembles them into a video.
2. The repository shall include a script that runs one representative fixed-contact case and writes a video to an output directory.
3. The rendered heatmap shall use a consistent color scale across all frames.
4. The script shall create its output directories if they do not exist and shall fail clearly if `ffmpeg` is unavailable.
5. The implementation shall preserve the existing solver test suite.

# Acceptance Criteria

1.1 A source-level function exists that accepts solver history plus output settings and writes a video.
1.2 The helper is documented with at least a one-line docstring signature.
2.1 A script under `scripts/` runs a small case with the current solver API.
2.2 Running the script produces a video file under `outputs/`.
3.1 The plotting code computes a single `(zmin, zmax)` range from the full history and reuses it for every frame.
4.1 The script or helper creates missing frame/video directories automatically.
4.2 If `ffmpeg` is not found, the error message states that explicitly.
5.1 `julia --project=. -e 'using Pkg; Pkg.test()'` still passes after the change.

# Testing Plan

- Add a small unit test for any pure visualization helper that is practical to test without rendering a full video.
- Run the full Julia test suite after integrating the new source file and dependency.
- Run the video script once as a smoke test and confirm it writes an MP4 file.

# Implementation Plan

1. Add the plotting dependency and create a focused visualization source file.
   Test: package resolves and existing tests still compile.
2. Add a small pure helper test, likely for consistent color-limit calculation from a solver history.
   Test: run the new focused test first.
3. Implement the visualization helper with frame rendering and `ffmpeg` assembly.
   Test: run the focused test and a manual smoke invocation on a tiny synthetic history if needed.
4. Add `scripts/run_case_and_make_video.jl` to define and run one representative case.
   Test: execute the script and confirm frame/video output exists.
5. Run the full Julia test suite.
   Test: `julia --project=. -e 'using Pkg; Pkg.test()'`.

# vertical-dynamics

`vertical-dynamics` is a Julia testbed for non-wetting impact models on a 2D surface grid. The current codebase focuses on a non-axisymmetric bath with a prescribed contact mask, a rigid-plane contact model, and structured Dirichlet-to-Neumann operators for the fluid coupling.

## Current Scope

The repository currently includes:

- a corrected nonperiodic DtN operator on a uniform Cartesian grid
- a prescribed-contact rigid-plane solver for non-axisymmetric bath dynamics
- dense small-grid reference assembly and a matrix-free GMRES solve path
- a demo script that runs one case and renders a top-view heatmap video of `η(x,y,t)`

The current solver is intended as a research testbed, not a fully validated production code.

## Repository Map

Start with:

- [docs/project-map.md](/Users/eaguerov/Documents/Github/vertical-dynamics/docs/project-map.md)
- [docs/source-tree.md](/Users/eaguerov/Documents/Github/vertical-dynamics/docs/source-tree.md)
- [docs/test-guide.md](/Users/eaguerov/Documents/Github/vertical-dynamics/docs/test-guide.md)

## Requirements

- Julia `1.12`
- `ffmpeg` in `PATH` if you want to render MP4 videos

Julia dependencies are managed through the local project:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## Running Tests

Run the full test suite with:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Running the Demo Case

Run one demo case and produce an MP4 heatmap video with:

```bash
julia --project=. scripts/run_case_and_make_video.jl
```

The script writes output under `outputs/` and prints the final MP4 path.

## Current Model Assumptions

The current coupled model assumes:

- the contact mask is prescribed
- the contact region is constrained to remain coplanar
- the contact pressure is solved as an unknown field on the mask
- the body motion is reduced to rigid-plane generalized coordinates
- the bath is discretized on a uniform Cartesian grid

The code does not yet solve for an unknown contact region.

## Example Output

After running the demo script, a typical output file looks like:

- [outputs/case_eta_20260318_095813.mp4](/Users/eaguerov/Documents/Github/vertical-dynamics/outputs/case_eta_20260318_095813.mp4)

## Notes

- The periodic DtN oracle was intentionally removed from the active API to keep the codebase focused on the nonperiodic path.
- The visualization path uses `Plots.jl` for frame rendering and `ffmpeg` for MP4 assembly.

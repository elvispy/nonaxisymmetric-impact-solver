# As Is

The repository currently exposes both the nonperiodic corrected-kernel DtN path and the older periodic DtN oracle through the public API. `VerticalDynamics` exports `PeriodicDtN2D` and `build_periodic_dtn`, and the test suite still contains periodic-specific checks. The periodic path is useful as historical reference but is no longer the intended production direction and now adds noise to the codebase.

# To Be

The current state of the repository should be frozen in a git commit so it can be recovered later if needed. After that, the public codebase should expose only the corrected nonperiodic DtN path. Periodic-only entry points, tests, and top-level documentation references should be removed so the remaining API and docs describe the solver that is actually being pursued.

# Requirements

1. The current project state shall be captured in a git commit before removing the periodic DtN entry point.
2. The package shall no longer export or expose the periodic DtN entry point as part of the main source API.
3. Periodic-only tests shall be removed from the active test suite.
4. Source-level and project-level documentation shall no longer present the periodic DtN route as a current supported entry point.
5. The remaining nonperiodic codebase shall continue to pass the test suite after the cleanup.

# Acceptance Criteria

1.1 A git commit exists that records the state prior to periodic-entry-point removal.
1.2 The commit message clearly describes the frozen solver baseline.
2.1 `src/VerticalDynamics.jl` no longer exports `PeriodicDtN2D` or `build_periodic_dtn`.
2.2 The DtN source no longer exposes the periodic constructor or apply methods in the main codepath.
3.1 `test/runtests.jl` no longer includes periodic-only tests.
3.2 No active test depends on `build_periodic_dtn` or `PeriodicDtN2D`.
4.1 `README.md` and `docs/source-tree.md` no longer describe the periodic DtN path as part of the active API.
4.2 Any remaining mentions of the periodic route are limited to historical design notes under `docs/plans/`.
5.1 `julia --project=. -e 'using Pkg; Pkg.test()'` exits successfully after the cleanup.

# Testing Plan

- Freeze the current tree in a commit before cleanup.
- Remove the periodic entry point and its test references.
- Run the full Julia test suite after cleanup.
- Spot-check the exported API and repository grep results to ensure the periodic path is no longer exposed from the active source tree.

# Implementation Plan

1. Create this task requirements note and inspect the current public exports, tests, and periodic references.
   Test: repository inspection only.
2. Commit the current project state with a baseline message.
   Test: `git log -1 --oneline`.
3. Remove periodic DtN exports, constructors, apply methods, and periodic-only tests.
   Test: grep the source and test tree for `PeriodicDtN2D` and `build_periodic_dtn`.
4. Update active documentation files to describe only the nonperiodic path.
   Test: grep project docs outside `docs/plans/` for periodic API references.
5. Run the full Julia test suite and confirm the cleaned codebase is green.
   Test: `julia --project=. -e 'using Pkg; Pkg.test()'`.

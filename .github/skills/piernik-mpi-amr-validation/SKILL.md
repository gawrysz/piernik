---
name: piernik-mpi-amr-validation
description: Select and interpret Piernik validation tests for numerical, MPI, AMR, load-balancing, restart, or physics changes.
---

# Piernik numerical and distributed validation

Read the root and relevant subsystem `AGENTS.md` files, then use the
[test map](../../../jenkins/README.md) and current Makefile to select existing
`artifacts` and `gold` targets matching the affected physics and invariants.
Start with focused targets; use broader groups when the change warrants them
and resources permit. Run `make QA` for code changes as directed by the project
instructions.

Before running a target, inspect its commands and existing output paths. The
test map describes important cleanup and `.setuprc` side effects. Do not delete
or overwrite user data; if a target would do so, ask before proceeding or use
an isolated environment. Do not run `make gold-clean` merely for convenience.

Distinguish a test failure from an informative `Differs!` accuracy metric.
Report the measured metrics and reference values, and investigate unexpected
changes without inventing numerical tolerances. Gold tests require zero total
`gdf_distance`, not byte-identical HDF5 files; never update a gold reference
solely to make a test pass.

If existing targets do not exercise the changed behavior, explain the gap and
add or run a focused regression test when appropriate. For Hilbert ordering,
face-adjacency can be checked on a complete uniform grid of power-of-two
square or cubic blocks, but is not a universal invariant for sparse samples
or non-power-of-two and rectangular blocks.

Report selected commands, results, meaningful metric differences, any skipped
tests and their reasons. Preserve logs needed to investigate failures.

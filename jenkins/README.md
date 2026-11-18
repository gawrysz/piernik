# Validation tests

Run the targets below from the repository root. See the [Makefile](../Makefile)
for the commands and current test lists; `make help` summarizes the groups.
`make CI` runs `QA`, `artifacts` and `gold`.

The project skill [`piernik-mpi-amr-validation`](../.github/skills/piernik-mpi-amr-validation/SKILL.md)
uses this map to select and interpret tests for numerical and distributed changes.

## Artifact tests

`make artifacts` runs the following targets. Individual targets can also be run
with `make <target>`.

| Target | Purpose |
| --- | --- |
| `dep` | Generate a module dependency graph; the default problem is `testing_and_debuging/chimaera`. |
| `py3` | Check that setup runs under Python 3. |
| `noHDF5` | Build a configuration without HDF5. |
| `I64` | Build with default 64-bit integers and conversion warnings treated as errors. |
| `IOv2` | Check restart reproducibility, including a restart with a different number of MPI processes. |
| `Jeans` | Evolve an oscillating self-gravitating system: periodic multigrid gravity, numerical accuracy and restart reproducibility. Sensitive to changes in precision and timestep. |
| `Maclaurin` | Compare a self-gravitating ellipsoid with the analytical potential on an AMR grid, using open boundaries and the multipole solver. |
| `3body` | Integrate a periodic figure-eight orbit: particle integration accuracy, self-gravity, N-body, open boundaries and particle restart reproducibility. |
| `CPAW` | Compare MHD evolution with the analytical circularly polarized Alfvén wave solution. |

`Jeans`, `Maclaurin`, `3body` and `CPAW` report measured error metrics and mark
deviations from stored values as `Differs!`. These deviations are informative:
they can reveal unwanted drift or an accuracy improvement and do not by
themselves make the target fail. The restart comparisons in `Jeans` and
`3body`, and the comparison in `IOv2`, require zero total difference according
to `gdf_distance`. Do not infer an acceptance tolerance from the displayed
reference values.

## Gold tests

`make gold` runs the named targets and `custom_gold`; each named target can
also be run separately, for example `make gold_MHDsedovAMR`. For an individual
target, ensure `jenkins/workspace/` exists first.

| Target | Purpose |
| --- | --- |
| `gold_MHDsedovAMR` | Detect changes in MHD results on an AMR grid, including MPI communication problems. |
| `gold_mcrtest` | Older, relatively simple cosmic-ray test. It could potentially become an analytical test of anisotropic diffusion of a Gaussian blob. |
| `gold_mcrwind` | Older cosmic-ray wind test exercising more modules than `mcrtest`. |
| `gold_CRESP` | Check spectral cosmic-ray diffusion over energy bins, for which this test has no analytical reference. |
| `gold_resist` | Exercise non-ideal MHD with electrical resistivity; this is the only existing test of that module. |
| `gold_streaming_instability` | Check a simplified streaming-instability model related to published work; its detailed coverage needs further investigation. |

Each gold test compares the current result with one generated from a specified
reference commit. It passes only if the total difference reported by
[`gdf_distance`](../bin/gdf_distance) is zero; it is not a byte-for-byte
comparison of the entire HDF5 files. A difference needs investigation, whether
or not the numerical change was intended. Never update `GOLD_COMMIT` merely
to make a test pass.

## Run-time effects and prerequisites

- `make artifacts` removes `jenkins/artifacts/` before running the group;
  individual targets also write artifacts there.
- Gold tests use and modify `jenkins/workspace/`, setup-generated object and
  run directories, and temporarily edit `.setuprc` when it exists. Inspect
  existing data and configuration before running them; do not assume all paths
  are fresh or safe to remove. A failed or interrupted run may need inspection.
- Gold tests may clone the reference repository and build both versions.
  Ensure adequate disk space and network access. Existing reference outputs
  can be reused.
- `make gold-clean` removes the gold workspace contents **and** artifacts.
  Use it only when their removal is intended.
- `make gold-serial` uses the older `gold_test_list.sh` runner; GNU Parallel is
  needed by that runner's parallel mode.

## Gold configuration

[`gold_test.sh`](gold_test.sh) reads a trusted file from
[`gold_configs/`](gold_configs/):

- `GOLD_COMMIT`: reference commit, updated only for a justified change.
- `PROBLEM_NAME`: problem under `problems/`.
- `SETUP_PARAMS`: setup options for both builds.
- `GOLD_PARAMS`: runtime arguments for the reference executable.
- `OUTPUT`: output file to compare.
- `NTHR`: number of MPI processes; defaults to 1.

The root Makefile maps named targets to configs and `custom_gold` discovers
other configs. The older `gold_test_list.sh` runner also discovers configs.

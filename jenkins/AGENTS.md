# jenkins/ – gold (regression) tests and CI glue

Parent: [../AGENTS.md](../AGENTS.md). User-level description: [README.md](README.md). `make gold` (part of `make CI`) calls `gold_test_list.sh`.

## How a gold test works
`gold_test.sh <config>` reads these variables from `gold_configs/<name>.config`: `GOLD_COMMIT` (SHA of a commit considered correct), `PROBLEM_NAME`,
`SETUP_PARAMS` (current code), `GOLD_PARAMS` (setup options for the old commit), `OUTPUT` (file to compare), optional `NTHR`.
It clones the repository (plus all `refs/pull/*`) into a workspace (`jenkins/workspace/`, git-ignored), checks out `GOLD_COMMIT`, builds and runs it, builds and runs the current tree
(using the **current** `setup` and `compilers/` – rsynced over the old checkout), and compares the HDF5 outputs with `bin/gdf_distance` (result mapped via log10 to a pass/fail threshold).
`gold_test_list.sh` runs all `gold_configs/*.config` in parallel (GNU `parallel`) or serially (`SERIAL=1`, `make gold-serial`).
Configs: `MHDsedovAMR`, `mcrtest_new`, `mcrtest_CRESP`, `mcrwind`, `resist`, `streaming_instability`.

## Updating a gold
A failing gold test means that results changed. If the change is intended (bug fix, new physics default), move `GOLD_COMMIT` to a commit
that contains the new behaviour, **explaining it in the commit message**. Never update it just to make CI green.
Adding a test: a config file + a `gold_<name>` target and an entry in `GOLD_TESTS` in the `Makefile` (`custom_gold` automatically runs every config whose name is *not* hard-coded in the Makefile case list, so a new config is picked up by `make gold` even without a target).

## Other
- `extract_jobs.sh` – dumps shell commands from Jenkins job XMLs (`$JENKINS_HOME`), for migrating the server configuration into the repo.

## Observations
- `bin/regenerate_golds.sh` documents the older approach (stored files on a private disk); the commit-based approach here replaced it.
- `gold_test.sh` parses configs with `eval $(sed ...)` – configs must be trusted files, quoting of values with spaces matters.
- `GOLD_COMMIT` must be pushed to GitHub (the script clones `$PIERNIK_REPO`), and shallow clones break checkout of old commits.

# Piernik Codebase Guide for AI Agents

## Project Overview

**Piernik** is a hydrodynamic simulation code written primarily in Fortran 90, designed for astrophysical MHD (magnetohydrodynamics) and fluid dynamics problems. It supports:
- Structured AMR (Adaptive Mesh Refinement) grids
- Parallel execution via MPI and OpenMP
- HDF5 I/O for checkpoints and outputs
- Multiple physics modules (fluids, magnetic fields, gravity, particles, dust)
- Problem-specific implementations in dedicated directories

**Key fact**: The build system uses a Python-based setup script that generates per-problem Makefiles from templates and dependency graphs.

## Critical Build & Setup Workflow

### Setup and Compilation

1. **Configuration**: Edit `.setuprc` files (compiler choice, flags, optimization) or create new compiler definitions in `compilers/`
   - Default: `compilers/h5pfc_gnu.in` (OpenMPI+HDF5+gfortran)
   - Each `.in` file defines: `F90`, `F90FLAGS`, `LIBS`, `LDFLAGS`, debug modes

2. **Setup phase** (`./setup` Python script):
   - Scans `src/` and problem-specific subdirectories for Fortran files
   - Builds dependency graphs using regex on module statements and `#include` directives
   - Auto-generates `obj_<PROBLEM>/Makefile` for each problem
   - Compiles if `--nocompile` is not set

3. **Build commands**:
   ```bash
   make                    # Rebuild all obj*/piernik executables
   make obj_ABC            # Build specific problem (obj_ABC/)
   RS=1 make obj_ABC       # Setup only, skip compile
   make resetup            # Re-run setup for all obj* using stored .setuprc
   make allsetup           # Create obj directories for all valid problems
   ```

4. **Important**: Each `obj_<PROBLEM>/` directory is independent. Files in `src/` are **symlinked** or copied into the obj directory during setup.

### Critical Compiler Flags

- **Default precision**: `-fdefault-real-8` (real*8 = 64-bit floats) — changing this breaks numerics
- **Debug mode**: Set `PIERNIK_DEBUG=1` in compiler config for stack traces, runtime checks
- **Gfortran 14.x fix**: `-fno-recursive -fcheck=recursion` to prevent CRESP initialization errors
- **MPI compatibility**: `-fallow-argument-mismatch` for gfortran ≥ 10 (old MPI interface)

## Architecture: Source Tree Organization

```
src/
├── base/              # Core framework: grid, time stepping, profiling
├── scheme/            # Numerical schemes (MUSCL, PPM, RK integration)
├── grid/              # Grid management, AMR refinement, boundaries
├── fluids/            # Equation of state, density, velocity
├── magnetic/          # MHD equations, divergence cleaning
├── gravity/           # Gravity solvers (multigrid, Poisson)
├── IO/                # I/O abstractions (HDF5, restart, particles)
├── particles/         # Particle tracer support
├── scheme/            # Flux/reconstruction schemes
└── utils/             # Utilities (constants, physical units)

problems/
├── <problem_name>/    # Each problem: problem.par, initproblem.F90, pierce.h
```

**Module dependency order matters**: Files in `base/` are foundational; later directories depend on them. The setup script uses `use` statements and `#include` directives to build dependency graphs.

## Key Patterns & Conventions

### Problem Definition Structure

Each problem in `problems/<name>/` requires:
- **`problem.par`**: Problem parameters (domain size, timesteps, physics flags)
- **`initproblem.F90`**: Initial conditions and problem-specific subroutines
- **`pierce.h`** (optional): C-style includes for build-time flags

### Fortran Conventions

- **Free-form**: `-ffree-form`, no fixed-column restrictions
- **Type declarations**: Must use `implicit none` everywhere
- **Module structure**:
  ```fortran
  module module_name
    implicit none
    ! declarations
  contains
    subroutine/function definitions
  end module module_name
  ```
- **File naming**: `<module_name>.F90` (always uppercase `.F90` for preprocessing)
- **Includes**: Use `#include "file.h"` for compile-time definitions

### Dependency Resolution

The setup script (`python/DirWalk.py`) scans for:
- `use <module>` statements → module dependencies
- `#include "<file.h>"` → header includes
- Modules pulled in from other files override local ones (checked by setup logic)

**Critical**: Module names must match file stems (e.g., `mod_global.f90` → `module mod_global`)

### Testing & Quality Assurance

- **`make qa`**: Run static analysis (checks error messages, license headers)
- **`make pycodestyle`**: Python code style checks
- **`make gold`**: Run gold tests (reference runs in `jenkins/workspace/`) with structured outputs
- **`make CI`**: Local equivalent of Jenkins CI pipeline

### Debugging

- Set `PIERNIK_DEBUG=1` in compiler config for gdb-friendly builds
- Check `bin/piernik` output format and error codes
- Stack traces appear with `-fbacktrace`

## Cross-Component Communication

### Global State & Parameters

- **`cg_list`**: Central grid structure (all blocks, refinement levels)
- **`mpi_setup`**: MPI communicators and rank/size
- **Global modules** in `src/base/`: Define shared constants and state

### Boundary Conditions

- Grid boundaries handle periodic/reflective/extrapolation
- Handled in `src/grid/` via standard interfaces
- Problem-specific BCs can override in `initproblem.F90`

### Data Flow for Initialization

1. Read `problem.par` parameters
2. Call `initproblem()` to set initial conditions
3. Grid and MPI topology initialized before physics
4. I/O system prepares checkpoint format (HDF5 by default)

## External Dependencies

- **Fortran compiler**: gfortran, ifort (Intel Fortran), or oneAPI compiler
- **MPI**: OpenMPI or MPICH (via wrapper compilers: `mpif90`, `h5pfc`)
- **HDF5**: Parallel HDF5 library with Fortran bindings
- **FFTW3**: For Fourier transforms (optional physics modules)
- **Python 3**: For setup script and post-processing utilities

## When Adding Features

1. **New physics module**: Create `src/<phys>/module_name.F90` with proper `use` dependencies on `base/` modules
2. **New problem**: Create `problems/<name>/` with `problem.par`, `initproblem.F90`
3. **Modify core**: Changes in `src/base/` or `src/scheme/` require re-setup of all problems
4. **Test**: Run `make allsetup` first, then selective builds with specific obj directories

## Commands for Common Tasks

```bash
# Quick build and test for one problem
./setup
make obj_tearing           # or any obj_*

# Check setup for all problems without building
make allsetup

# Find module dependencies for a problem
make dep P=problems/tearing/initproblem.F90

# Run full CI checks locally
make CI

# Rebuild with different compiler
cat > .setuprc << EOF
COMPILER=compilers/ifx.in
EOF
make resetup
```

## Known Gotchas

- **Real precision**: `-fdefault-real-8` is mandatory; 32-bit floats break physics
- **Gfortran recursion bug (14.x)**: Use `-fno-recursive -fcheck=recursion` to avoid CRESP crashes
- **Module name case**: Fortran is case-insensitive; setup script treats `MODULE MOD_X` and `module mod_x` as equivalent
- **Parallel I/O**: HDF5 setup must be compiled with MPI; check with `h5pfc -showversion`
- **obj_<PROBLEM> independence**: Changing src files requires re-setup; `make resetup` uses saved .setuprc

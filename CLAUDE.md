# Claude Instructions for MorFit.jl

## Overview

MorFit.jl is the core Julia package for molecular fitting using Monte Carlo methods and morphometric energy functions.

## Permissions & Skills

See [../CLAUDE.md](../CLAUDE.md) for allowed commands, permissions, and available skills.

## Key Files

- `ARCHITECTURE.md` - Repository structure reference
- `src/MorFit.jl` - Package entry point
- `tests/Tests.jl` - Test runner

## Data Format

**Simulation state**: `Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}` - one (R, t) per molecule

**Templates** (unified list format, supports heterogeneous molecules):
- `centers::Vector{Matrix{Float64}}` - 3×N matrix per molecule
- `radii::Vector{Vector{Float64}}` - radius vector per molecule

**MOLECULE_DATA** - Centered building blocks (JLD2):
```julia
MorFit.MOLECULE_DATA["4ty7:protein"].centers  # 3×N Matrix
MorFit.MOLECULE_DATA["4ty7:protein"].radii    # Vector
```

**EXPERIMENTAL_ASSEMBLIES** - Realized coordinates for RMSD (JLD2):
```julia
key = ["4ty7:protein", "4ty7:ligand"]
assembly = MorFit.EXPERIMENTAL_ASSEMBLIES[key][1]  # First equivalent assembly
assembly.centers[1]  # 3×N Matrix - REALIZED coords for molecule 1
assembly.centers[2]  # 3×N Matrix - REALIZED coords for molecule 2
assembly.radii[1]    # Vector{Float64} - radii for molecule 1
```

## Running Tests

```julia
using MorFit
MorFit.Tests.run()  # All tests
MorFit.Tests.run_morphometric_approach_tests()  # Working tests only
```

## Related Repositories

- `HPC_MorFit` - HPC simulation runner (uses MorFit.jl)
- `Notebooks_MorFit` - Visualization and analysis notebooks

## Recent Changes (2026-02-02)

### Realization Functions Refactored
The 6 separate realization functions were consolidated into a single `get_realization` function:

```julia
# Old API (removed):
get_flat_realization(x, templates)
get_point_vector_realization(x, templates)
get_matrix_realization_per_mol(x, templates)
# ... etc

# New API:
get_realization(x, templates)                      # default: format=:matrix
get_realization(x, templates, format=:flat)        # Vector{Float64}
get_realization(x, templates, format=:points)      # Vector{Vector{Float64}}
get_realization(x, templates, format=:points_per_mol)
get_realization(x, templates, format=:point3f)
get_realization(x, templates, format=:point3f_per_mol)
```

### Algorithm Module Improvements
- Consolidated 3 type-specific `add_to_output` functions into 1 generic version
- Added timestamps to resume simulation (cumulative across sessions)
- Added comprehensive docstrings to `RandomWalkMetropolis`, `ConnectedComponentRandomWalkMetropolis`, and `simulate!` methods

### Removed Unused Functions
- `get_moment_of_inertia` (had undefined variable bug)
- `perturb_single_randomly_chosen_only_translations` (had off-by-one bug)
- `get_initial_state_only_translations` (inconsistent return type)

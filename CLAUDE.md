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

**EXPERIMENTAL_ASSEMBLIES** - Reference states for RMSD (JLD2):
```julia
key = ["4ty7:protein", "4ty7:ligand"]
state = MorFit.EXPERIMENTAL_ASSEMBLIES[key][1]  # First equivalent state
R, t = state[1]  # First molecule's rotation and translation
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

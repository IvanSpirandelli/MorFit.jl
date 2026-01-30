# Claude Instructions for MorFit.jl

## Overview

MorFit.jl is the core Julia package for molecular fitting using Monte Carlo methods and morphometric energy functions.

## Permissions

See [../CLAUDE.md](../CLAUDE.md) for allowed commands and permissions.

## Key Files

- `ARCHITECTURE.md` - Repository structure reference
- `src/MorFit.jl` - Package entry point
- `tests/Tests.jl` - Test runner

## Data Format

Templates use unified list format (supports heterogeneous molecules):
- `template_centers::Vector{Matrix{Float64}}`
- `template_radii::Vector{Vector{Float64}}`

## Running Tests

```julia
using MorFit
MorFit.Tests.run()  # All tests
MorFit.Tests.run_morphometric_approach_tests()  # Working tests only
```

## Related Repositories

- `HPC_MorFit` - HPC simulation runner (uses MorFit.jl)
- `Notebooks_MorFit` - Visualization and analysis notebooks

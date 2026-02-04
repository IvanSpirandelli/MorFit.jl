# Claude Instructions for MorFit.jl

## Overview

MorFit.jl is the core Julia package for molecular fitting using Monte Carlo methods and morphometric energy functions.

## Permissions & Skills

See [../CLAUDE.md](../CLAUDE.md) for allowed commands, permissions, and available skills.

**Important:** Before starting any task, check if a relevant skill applies (e.g., `/morfit-refactor` for refactoring, `/julia-review` for code review).

## Key Files

- `ARCHITECTURE.md` - Repository structure and data formats
- `src/MorFit.jl` - Package entry point
- `src/modules/Energies/src/types.jl` - Parameter struct definitions
- `tests/Tests.jl` - Test runner

## Data Format

**Simulation state**: `Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}` - one (R, t) per molecule

**Templates** (unified list format):
- `centers::Vector{Matrix{Float64}}` - 3×N matrix per molecule
- `radii::Vector{Vector{Float64}}` - radius vector per molecule

**Global data** (loaded from JLD2):
- `MorFit.MOLECULE_DATA[key]` - Centered building blocks (see `ARCHITECTURE.md:67-76`)
- `MorFit.EXPERIMENTAL_ASSEMBLIES[key]` - Realized coordinates for RMSD (see `ARCHITECTURE.md:78-95`)

## Running Tests

```bash
julia --project=. -e 'using MorFit; MorFit.Tests.run()'
```

Test subsets: `run_morphometric_approach_tests()`, `run_energy_call_tests()`, `run_configuration_distance_tests()`

## Energy Model

The total energy uses an additive model: `E = θ_G·F_sol + θ_O·O + θ_T·Σ(λᵢ·Pᵢ)`

### Parameter Types

Defined in `src/modules/Energies/src/types.jl`:
- `EnergyScales(θ_G, θ_O, θ_T)` - Scaling factors (default 1.0)
- `SolvationParams(rs, η)` - Probe radius, packing fraction
- `LinearOverlapParams(jump, slope)` - Overlap penalty
- `TopologyParams(λ)` - Persistence weights [λ₀, λ₁, λ₂]
- `NumericalParams(delaunay_eps, exact_delaunay)` - Tolerances
- `MolecularSystem(centers, radii, bounds)` - System definition

### Energy Functions

- `calculate_combined_energy(...)` - New additive model (recommended)
- `calculate_combined_potential(...)` - Legacy μ-interpolation (deprecated)

See `src/modules/Energies/src/combined.jl` for implementation.

## Realization API

Single function with format parameter (see `src/modules/Utilities/src/realizations.jl`):

```julia
get_realization(x, templates)                    # default: format=:matrix
get_realization(x, templates, format=:flat)      # Vector{Float64}
get_realization(x, templates, format=:points)    # Vector{Vector{Float64}}
```

## Related Repositories

- `mor-fit-hpc` - HPC simulation runner (uses MorFit.jl)
- `mor-fit-analysis` - Visualization and analysis notebooks

## Recent Changes (2026-02-03)

- New additive energy model: `E = θ_G·F_sol + θ_O·O + θ_T·T`
- Added `types.jl` with parameter structs
- `calculate_combined_energy` replaces μ-interpolation

## Recent Changes (2026-02-02)

- Consolidated 6 realization functions into single `get_realization` with format parameter
- Consolidated 3 `add_to_output` functions into 1 generic version
- Added timestamps to resume simulation
- Removed unused functions with bugs

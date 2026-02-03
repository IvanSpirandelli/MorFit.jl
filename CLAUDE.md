# Claude Instructions for MorFit.jl

## Overview

MorFit.jl is the core Julia package for molecular fitting using Monte Carlo methods and morphometric energy functions.

## Permissions & Skills

See [../CLAUDE.md](../CLAUDE.md) for allowed commands, permissions, and available skills.

**Important:** Before starting any task, check if a relevant skill applies (e.g., `/morfit-refactor` for refactoring, `/julia-review` for code review). Use `/` to see available skills.

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

- `mor-fit-hpc` - HPC simulation runner (uses MorFit.jl)
- `mor-fit-analysis` - Visualization and analysis notebooks

## Energy Model

The total energy uses an additive model with independent scaling factors:

```
E = θ_G · F_sol + θ_O · O + θ_T · Σ(λᵢ · Pᵢ)
```

### Parameter Types (MorFit.Energies)

```julia
using MorFit.Energies

# Energy scaling (defaults to 1.0)
scales = EnergyScales(θ_G=1.0, θ_O=1.2, θ_T=0.5)

# Solvation parameters (prefactors computed via White Bear)
sol = SolvationParams(rs=1.4, η=0.3665)

# Overlap penalty
ol = LinearOverlapParams(jump=1.0, slope=100.0)

# Topology (persistence weights for H₀, H₁, H₂)
topo = TopologyParams([1.0, 0.0, 0.0])

# Numerical tolerances
num = NumericalParams(delaunay_eps=1.0, exact_delaunay=false)

# Molecular system
system = MolecularSystem(centers, radii, bounds)
```

### Energy Functions

```julia
# New additive model (recommended)
energy, measures = calculate_combined_energy(
    x, system, sol_params, ol_params, topo_params,
    num_params, scales, precomputed, bol_check
)

# Legacy μ-interpolation (deprecated)
energy, measures = calculate_combined_potential(x, ..., μ, ...)
```

## Recent Changes (2026-02-03)

### Energy System Refactored
- New additive energy model: `E = θ_G·F_sol + θ_O·O + θ_T·T`
- Added `types.jl` with parameter structs
- `calculate_combined_energy` replaces μ-interpolation
- Legacy functions preserved for backward compatibility

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

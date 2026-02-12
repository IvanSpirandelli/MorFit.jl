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

## Simulation I/O

Defined in `src/modules/Algorithms/src/simulation_io.jl`:

### SimulationOutput

Structured output container replacing manual `Dict{String, Vector}` initialization:

```julia
output = Algorithms.SimulationOutput()
simulate!(rwm, x_init, 60.0, 10000, output)
output.total_step_attempts  # scalar Int
output.E_total              # Vector{Float64}
output.measures["E_G"]      # Vector{Float64} (auto-created)
```

- `record!(output, E, state, step, measures)` — auto-expanding, no pre-initialization
- `to_dict(output)` → `Dict{String, Vector}` for JLD2 serialization
- `SimulationOutput(d::Dict{String, Vector})` — reconstruct from saved Dict (legacy compatible)

### build_input_dict

Canonical function to assemble the input/config dictionary for saving:

```julia
saved_config = Algorithms.build_input_dict(
    molecule_ids=..., system=..., sol_params=..., ol_params=...,
    topo_params=..., num_params=..., scales=..., pert_params=...,
    perturbation=..., T_sim=..., wall_clock_limit_minutes=...,
    target_iterations=..., x_init=..., seed=...,
)
```

### Resume API

```julia
# Fresh start
simulate!(rwm, x_init, wall_clock_limit, target_iters, output)
# Resume from previous output
simulate!(rwm, wall_clock_limit, target_iters, output)
```

## Related Repositories

- `mor-fit-hpc` - HPC simulation runner (uses MorFit.jl)
- `mor-fit-analysis` - Visualization and analysis notebooks

## Recent Changes (2026-02-12)

- Replaced `Dict{String, Vector}` output with `SimulationOutput` struct
- Removed `add_to_output` — replaced by `record!` with auto-expanding measures
- Added `build_input_dict` to centralize input dict assembly
- Resume API: `simulate!(rwm, time, iters, output)` instead of passing Dict
- Backward compatible: `to_dict`/`SimulationOutput(dict)` round-trip preserves JLD2 format

## Recent Changes (2026-02-03)

- New additive energy model: `E = θ_G·F_sol + θ_O·O + θ_T·T`
- Added `types.jl` with parameter structs
- `calculate_combined_energy` replaces μ-interpolation

## Recent Changes (2026-02-02)

- Consolidated 6 realization functions into single `get_realization` with format parameter
- Consolidated 3 `add_to_output` functions into 1 generic version
- Added timestamps to resume simulation
- Removed unused functions with bugs

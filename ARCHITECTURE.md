# MorFit.jl Architecture

## Overview

Core Julia package for molecular fitting using Monte Carlo methods and morphometric energy functions.

## Directory Structure

```
MorFit.jl/
├── src/
│   ├── MorFit.jl                  # Package entry point
│   ├── modules/
│   │   ├── Utilities/             # State handling, RMSD, initialization
│   │   ├── Algorithms/            # Monte Carlo samplers (RWM, CC-RWM)
│   │   └── Energies/              # Energy functions (solvation, morphometric)
│   └── templates/                 # Protein structure data (JLD2 binary)
│       ├── molecule_data.jld2         # MOLECULE_DATA: centered templates
│       └── experimental_assemblies.jld2  # EXPERIMENTAL_ASSEMBLIES: ground truth
├── tests/
│   ├── Tests.jl                   # Test runner
│   ├── test_morphometric_approach.jl  # ✓ Working
│   ├── test_energy_calls.jl           # Stubbed - needs implementation
│   ├── test_configuration_distances.jl # ✓ Working (RMSD tests)
│   └── test_simulation.jld2           # Test simulation data for RMSD tests
├── Project.toml
└── LICENSE (MIT)
```

---

## Core Modules

### Utilities
| File | Purpose |
|------|---------|
| `center_of_mass_computation.jl` | Center-of-mass calculations |
| `configuration_distances.jl` | RMSD, alignment, metrics to ground truth |
| `initialization.jl` | Initial state generation |
| `perturbation.jl` | State perturbation for MC moves |
| `realizations.jl` | Realize configurations from state vectors |
| `state_conversion.jl` | State representation conversions |
| `temperature.jl` | Temperature calculations |

### Algorithms
| File | Purpose |
|------|---------|
| `random_walk_metropolis.jl` | Standard RWM with time/iteration limits + resume |
| `connected_component_random_walk_metropolis.jl` | RWM for n_mol > 2 |

### Energies
| File | Purpose |
|------|---------|
| `solvation_free_energy.jl` | Solvation energy calculations |
| `connected_component_calculations.jl` | CC-wise energy (supports heterogeneous molecules) |
| `combined.jl` | Combined fsol + twasp potentials |
| `persistence_computations.jl` | Persistent homology metrics |
| `alpha_shapes/alpha_shape_diagrams.jl` | Alpha-complex via Python |
| `morphometric_approach/ball_union_measures.jl` | Union of balls calculations |
| `morphometric_approach/prefactors.jl` | Pre-computed energy factors |

---

## Global Data Structures

### MOLECULE_DATA
Centered molecular building blocks (templates) loaded from JLD2:
```julia
MorFit.MOLECULE_DATA["4ty7:protein"]  # Returns NamedTuple
  .centers  # 3×N Matrix{Float64} - atom coordinates CENTERED AT ORIGIN (COM = [0,0,0])
  .radii    # Vector{Float64} - atomic radii
```
- Keys follow format: `"pdb_id:component"` (e.g., `"4ty7:protein"`, `"4ty7:ligand"`)
- Each entry is ONE molecule template, centered at the origin
- Used as building blocks for simulations (apply rotation + translation to position)

### EXPERIMENTAL_ASSEMBLIES
Ground truth configurations for RMSD calculation. Stores **realized coordinates** (already positioned in space):
```julia
key = ["4ty7:protein", "4ty7:ligand"]
assemblies = MorFit.EXPERIMENTAL_ASSEMBLIES[key]  # Returns Vector of equivalent assemblies

# Each assembly is a NamedTuple:
assembly = assemblies[1]
assembly.centers  # Vector{Matrix{Float64}} - one 3×N matrix per molecule (REALIZED positions)
assembly.radii    # Vector{Vector{Float64}} - one radius vector per molecule

# Access individual molecules:
assembly.centers[1]  # 3×N Matrix - protein coordinates (already in final position)
assembly.centers[2]  # 3×N Matrix - ligand coordinates (already in final position)
```
- Keys are `Vector{String}` of molecule IDs specifying which templates are used
- Values are `Vector` of equivalent assemblies (for permutation symmetry in homodimers/trimers)
- Coordinates are **pre-realized** (not centered) - ready for direct RMSD comparison

### Data Flow
```
MOLECULE_DATA (centered templates)     EXPERIMENTAL_ASSEMBLIES (realized)
         │                                        │
         │ Apply R * template + t                 │ Direct comparison
         ▼                                        ▼
   Simulation state  ──────────────────────►  RMSD calculation
```

---

## Data Format

Templates use unified list format:
- `template_centers::Vector{Matrix{Float64}}` - one 3×N matrix per molecule
- `template_radii::Vector{Vector{Float64}}` - one radius vector per molecule

This supports both same-type assemblies (sta) and heterogeneous assemblies (ppii).

---

## Dependencies

| Package | Purpose |
|---------|---------|
| `AlphaMolWrapper_jll` | Alpha shape computation |
| `GeometryBasics` | 3D geometry primitives |
| `PyCall` | Python interop for alpha_shape_diagrams.jl |
| `Rotations` | Rotation matrices |
| `StaticArrays` | Performance |
| `Distributions` | Probability distributions |
| `JLD2` | Binary serialization |
| `Combinatorics` | Permutations for RMSD |
| `Graphs` | Connected component detection |
| `PDBTools` | PDB file handling |

---

## HPC_MorFit Integration

HPC_MorFit depends on MorFit.jl for core functionality. Key exports used:

**Algorithms:**
- `RandomWalkMetropolis`, `ConnectedComponentRandomWalkMetropolis`
- `simulate!` (both fresh start and resume variants)

**Energies:**
- `get_prefactors`, `get_bounding_radii`, `get_single_subunit_energy_and_measures`
- `solvation_free_energy_and_measures_*` functions
- `connected_component_wise_solvation_free_energy_and_measures`
- `calculate_combined_potential`

**Utilities:**
- `get_initial_state`, `perturb_single_randomly_chosen`
- `get_center_of_mass`, `get_realization`
- `get_rmsd_for_fixed_target_inhibitor_pair`

**Data:**
- `MOLECULE_DATA` - Centered molecular templates
- `EXPERIMENTAL_ASSEMBLIES` - Ground truth configurations for RMSD
- `get_reference_templates(molecule_ids)` - Helper to extract templates

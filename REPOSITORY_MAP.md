# MorFit.jl Repository Map

**Purpose:** Documentation for repository cleanup
**Total Size:** ~150 MB (96.7% is accumulated output data)

---

## Quick Summary

| Directory | Size | Status |
|-----------|------|--------|
| `assets/output/` | 145 MB | **DELETE** - Old simulation outputs |
| `src/templates/` | 2.3 MB | Review - Large data embedded in .jl files |
| `examples/` | 804 KB | Update - Uses old module name "MorphoMol" |
| `src/modules/` | ~100 KB | **KEEP** - Core code |
| `tests/` | 28 KB | **KEEP** - Has some commented code |

---

## Directory Structure

```
MorFit.jl/
├── src/                           # Main source code (2.8 MB)
│   ├── MorFit.jl                  # Package entry point
│   ├── simulation_setup.jl        # Factory functions for simulations
│   ├── modules/                   # Core functionality
│   │   ├── Utilities/             # State handling, RMSD, initialization
│   │   ├── Algorithms/            # Monte Carlo samplers
│   │   └── Energies/              # Energy functions (solvation, morphometric)
│   └── templates/                 # Protein structure data (2.3 MB)
│       ├── target_and_inhibitor_templates.jl    (44 KB)
│       ├── experimental_assembly.jl             (389 KB)
│       ├── asymmetric_unit_templates.jl         (756 KB)
│       └── protein_ligand_data.jld2             (1.5 MB binary)
├── tests/                         # Test suite (28 KB)
├── examples/                      # Jupyter notebooks (804 KB)
├── assets/                        # OUTPUT DATA (145 MB)
│   └── output/                    # 474 .poly files - OLD SIMULATION RESULTS
├── Project.toml                   # Julia dependencies
├── Manifest.toml                  # Dependency lock (shouldn't be committed)
├── LICENSE                        # MIT
└── Readme.md                      # Minimal (just header)
```

---

## Detailed Breakdown

### 1. `src/modules/` - Core Code (KEEP)

#### Utilities (~500 lines)
| File | Lines | Purpose |
|------|-------|---------|
| `center_of_mass_computation.jl` | 72 | Center-of-mass calculations |
| `configuration_distances.jl` | 211 | RMSD, alignment, metrics to ground truth |
| `initialization.jl` | 12 | Initial state generation |
| `perturbation.jl` | 30 | State perturbation |
| `realizations.jl` | 61 | Realize configurations from state vectors |
| `state_conversion.jl` | 12 | State representation conversions |
| `temperature.jl` | 51 | Temperature calculations |

#### Algorithms (~600 lines)
| File | Lines | Purpose |
|------|-------|---------|
| `random_walk_metropolis.jl` | 121 | Standard RWM algorithm |
| `simulated_annealing.jl` | 159 | Simulated annealing |
| `connected_component_random_walk_metropolis.jl` | ? | RWM for connected components |
| `hamiltonian_monte_carlo.jl` | ? | HMC sampling |

#### Energies (~1600 lines)
| File | Lines | Purpose |
|------|-------|---------|
| `solvation_free_energy.jl` | 218 | Solvation energy calculations |
| `connected_component_calculations.jl` | 235 | Connected component energy |
| `combined.jl` | 109 | Combined energy potentials |
| `persistence_computations.jl` | 100 | Persistent homology metrics |
| `alpha_shapes/alpha_shape_diagrams.jl` | 61 | Alpha-complex computations |
| `morphometric_approach/ball_union_measures.jl` | 159 | Union of balls calculations |
| `morphometric_approach/prefactors.jl` | 97 | Pre-computed factors |

---

### 2. `src/templates/` - Data Files (REVIEW)

These are large data files embedded as Julia code:

| File | Size | Content |
|------|------|---------|
| `target_and_inhibitor_templates.jl` | 44 KB | `TARGETS` dict with protein structures |
| `experimental_assembly.jl` | 389 KB | `EXPERIMENTAL_ASSEMBLY` dict - ground truth conformations |
| `asymmetric_unit_templates.jl` | 756 KB | Asymmetric unit definitions |
| `protein_ligand_data.jld2` | 1.5 MB | Binary pre-computed data |

**Consideration:** Could extract to JSON/JLD2 instead of .jl files for cleaner separation.

---

### 3. `tests/` - Test Suite (KEEP, minor cleanup)

| File | Size | Status |
|------|------|--------|
| `Tests.jl` | 667 B | Master test runner |
| `test_morphometric_approach.jl` | ? | Active |
| `test_interface.jl` | 1.4 KB | **Has commented/legacy code** |
| `test_energy_calls.jl` | 5.0 KB | Active |
| `test_configuration_distances.jl` | 1.3 KB | Active |

---

### 4. `examples/` - Notebooks (UPDATE OR DELETE)

| Notebook | Size | Issue |
|----------|------|-------|
| `rotating_needle.ipynb` | 236 KB | Uses "MorphoMol" (old name) |
| `proteins_in_solvation.ipynb` | 508 KB | Uses "MorphoMol" (old name) |
| `moving_point.ipynb` | 52 KB | Uses "MorphoMol" (old name) |
| `connected_components.ipynb` | 8 KB | Uses "MorphoMol" (old name) |

**Issue:** All notebooks reference old module name "MorphoMol" instead of "MorFit"

---

### 5. `assets/output/` - OLD DATA (DELETE)

- **Size:** 145 MB
- **Content:** 474 .poly files (polytope data)
- **Format:** ASCII with 3D coordinates, radii, RGBA colors
- **Purpose:** Old simulation outputs from PhD experiments
- **Status:** Not used by any current code

**Recommendation:** Delete entirely or move to external archive.

---

## Issues Found

### Critical
1. **assets/output/** - 145 MB of unused simulation outputs bloating repo

### Medium
2. **.gitignore mismatch** - Ignores `assets/`, `examples/`, `Manifest.toml` but they're committed
3. **Old module name** - Examples use "MorphoMol" instead of "MorFit"
4. **Manifest.toml committed** - Generally shouldn't be in version control for packages

### Minor
5. **Readme.md empty** - Just contains header
6. **.DS_Store files** - macOS metadata scattered around
7. **test_interface.jl** - Contains commented/legacy tests

---

## Cleanup Recommendations

### Phase 1: Immediate Space Recovery
```bash
# Remove 145 MB of old simulation outputs
rm -rf assets/

# Remove macOS metadata
find . -name ".DS_Store" -delete

# Remove lock file (regenerates on install)
rm Manifest.toml
```

### Phase 2: Fix Consistency
- Update `.gitignore` to match what should actually be ignored
- Either update examples to use "MorFit" or delete them
- Decide: should `examples/` be in the repo?

### Phase 3: Code Cleanup
- Review `test_interface.jl` for dead code
- Consider extracting template data from .jl to .jld2 format
- Add proper Readme.md documentation

---

## Dependency Analysis

**Project.toml dependencies:**
- `AlphaMolWrapper_jll` - Alpha shape computation
- `GeometryBasics` - 3D geometry
- `PyCall` - Python interop (purpose unclear)
- `Rotations` - Rotation matrices
- `StaticArrays` - Performance
- `Distributions` - Probability
- `JLD2` - Binary serialization
- `Combinatorics` - Combinatorial algorithms
- `Graphs` - Graph algorithms
- `PDBTools` - PDB file handling

**Note:** PyCall is required for `alpha_shape_diagrams.jl` - computes alpha shape diagrams via Python.

---

## File Counts by Type

| Extension | Count | Total Size |
|-----------|-------|------------|
| `.poly` | 474 | 145 MB |
| `.jl` | 24 | ~200 KB |
| `.ipynb` | 4 | 804 KB |
| `.jld2` | 1 | 1.5 MB |
| `.toml` | 2 | 25 KB |
| `.md` | 1 | 12 B |

---

## Decision Checklist

Before cleanup, decide:

- [x] Delete `assets/output/` entirely? (145 MB saved) - DONE
- [x] Keep or delete `examples/`? - DELETED
- [x] If keeping examples, update to use "MorFit"? - N/A (deleted)
- [x] Remove `Manifest.toml` from git? - DONE
- [ ] Keep template data as .jl files or convert to .jld2?
- [x] What is PyCall used for - keep or remove? - KEEP (alpha_shape_diagrams.jl)
- [ ] Add proper documentation to Readme.md?

---

## HPC_MorFit Dependency Analysis

The `../HPC_MorFit` repository is an HPC framework for running molecular simulations and Bayesian optimization. It depends on MorFit.jl for core functionality.

### Files Using MorFit

| File | Usage Level |
|------|-------------|
| `julia_scripts/start_simulations.jl` | Heavy - core simulation runner |
| `gaussian_process/objective_function.jl` | Light - RMSD calculations |
| `gaussian_process/batch_launcher.jl` | Import only |

### Required MorFit Functions

**From `MorFit` (top-level):**
- `get_initialization(input, flag)` - Initialize molecular assembly
- `get_energy(input)` - Get energy calculation function
- `get_perturbation(input)` - Get perturbation function for MC moves
- `are_bounding_spheres_overlapping(x, id1, id2, radius)` - Collision detection
- `get_bounding_radius(centers, radii, rs)` - Bounding sphere calculation
- `get_initial_connected_component_energies(...)` - Multi-molecule energy init

**From `MorFit.Algorithms`:**
- `RandomWalkMetropolis(energy, perturbation, β)` - Standard RWM sampler
- `ConnectedComponentRandomWalkMetropolis(...)` - RWM for n_mol > 2
- `simulate!(rwm, x_init, time_minutes, target_iters, output)` - Run simulation

**From `MorFit.Energies`:**
- `get_prefactors(rs, η)` - Energy prefactors from radii and packing fraction

**From `MorFit.Utilities`:**
- `get_center_of_mass(points, elements)` - Center of mass calculation
- `get_rmsd_for_fixed_target_inhibitor_pair(...)` - RMSD with fixed alignment

**Global Data Accessed:**
- `MorFit.TEMPLATES[mol_id]["template_centers"]` - Molecular template coordinates
- `MorFit.TEMPLATES[mol_id]["template_radii"]` - Molecular template radii
- `MorFit.PROTEIN_LIGAND_DATA[id]["protein"]["centers"]` - Protein coordinates
- `MorFit.PROTEIN_LIGAND_DATA[id]["protein"]["radii"]` - Protein radii
- `MorFit.PROTEIN_LIGAND_DATA[id]["protein"]["elements"]` - Protein elements
- `MorFit.PROTEIN_LIGAND_DATA[id]["ligand"]["centers"]` - Ligand coordinates
- `MorFit.PROTEIN_LIGAND_DATA[id]["ligand"]["radii"]` - Ligand radii
- `MorFit.PROTEIN_LIGAND_DATA[id]["ligand"]["elements"]` - Ligand elements

### Cleanup Implications

**MUST KEEP** for HPC_MorFit compatibility:
- `src/simulation_setup.jl` - provides `get_initialization`, `get_energy`, `get_perturbation`
- `src/modules/Algorithms/` - RWM and CC-RWM samplers
- `src/modules/Energies/src/morphometric_approach/prefactors.jl` - `get_prefactors`
- `src/modules/Utilities/src/center_of_mass_computation.jl` - `get_center_of_mass`
- `src/modules/Utilities/src/configuration_distances.jl` - RMSD functions
- `src/templates/target_and_inhibitor_templates.jl` - `TEMPLATES` dict
- `src/templates/protein_ligand_data.jld2` - `PROTEIN_LIGAND_DATA`

**Potentially removable** (not used by HPC_MorFit):
- `src/modules/Energies/src/persistence_computations.jl`
- `src/modules/Energies/src/alpha_shapes/` (unless used by energy functions)
- `src/templates/experimental_assembly.jl`
- `src/templates/asymmetric_unit_templates.jl`
- `src/modules/Algorithms/src/simulated_annealing.jl`
- `src/modules/Algorithms/src/hamiltonian_monte_carlo.jl`

---

## Changes Made

_Track all cleanup actions here:_

| Date | Action | Details |
|------|--------|---------|
| 2026-01-30 | Created | Initial repository mapping |
| 2026-01-30 | Deleted | `assets/` - 145 MB of old .poly simulation outputs |
| 2026-01-30 | Deleted | All `.DS_Store` files |
| 2026-01-30 | Deleted | `Manifest.toml` |
| 2026-01-30 | Deleted | `examples/` - 4 outdated notebooks using old module name |
| 2026-01-30 | Documented | HPC_MorFit dependency analysis - identified required functions |
| 2026-01-30 | Refactored | **Unified sta/ppii data formats** - Both assembly types now use list format: `template_centers::Vector{Matrix{Float64}}`, `template_radii::Vector{Vector{Float64}}`. For sta, this is n_mol copies of the same template. Changes: HPC `_get_template_centers_and_radii`, MorFit `_get_realization_radii_and_sizes`, added Vector overloads in `connected_component_calculations.jl` |
| 2026-01-30 | Deleted | `hamiltonian_monte_carlo.jl` and `simulated_annealing.jl` - Unused algorithms |
| 2026-01-30 | Deleted | All old Matrix{Float64} format overloads from `connected_component_calculations.jl`, `realizations.jl`, `solvation_free_energy.jl` - Only unified Vector{Matrix} format remains |

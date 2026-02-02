# MorFit.jl Potential Improvements

## Medium Priority

### 1. Parameter Explosion in Energy Functions
**File:** `src/modules/Energies/src/combined.jl`

Functions with 13-18 parameters are hard to use correctly. Consider:
- Grouping into structs (already have `SolvationParams`, etc. in HPC_MorFit)
- Using keyword arguments with defaults

### 2. Hardcoded Default in `perturb_single_specified`
**File:** `src/modules/Utilities/src/perturbation.jl:12`

```julia
function perturb_single_specified(x, σ_r, σ_t; specified_index = 2)
```

The default `specified_index = 2` assumes protein-ligand format. Document this assumption or make it required.

### 3. Persistence Calculation Duplication
**File:** `src/modules/Energies/src/persistence_computations.jl`

`_calculate_persistence_energy` (line 90) and `_calculate_persistence_energy_and_measures` (line 60) have nearly identical logic. Refactor:

```julia
function _calculate_persistence_energy_and_measures(pdgm, persistence_weights)
    # ... calculation ...
    return energy, measures
end

function _calculate_persistence_energy(pdgm, persistence_weights)
    energy, _ = _calculate_persistence_energy_and_measures(pdgm, persistence_weights)
    return energy
end
```

---

## Low Priority

### 1. Silent Failures with `@warn`
**Files:**
- `center_of_mass_computation.jl:28-30` (unknown elements)
- `MorFit.jl:52-54` (missing assemblies)

Consider making these errors or at least logging more context.

### 2. Assertions in Production Code
**Files:**
- `connected_component_calculations.jl:37`

Replace with proper error handling:

```julia
# Instead of:
@assert haskey(input, "simulation_time_minutes")

# Use:
if !haskey(input, "simulation_time_minutes")
    throw(ArgumentError("input must contain 'simulation_time_minutes'"))
end
```

### 3. No Input Validation in Prefactor Functions
**File:** `src/modules/Energies/src/morphometric_approach/prefactors.jl`

Functions like `pressure(rs, η)` don't validate that `η ∈ (0, 1)`. Could produce NaN/Inf.

### 4. Performance: Repeated Time Computation
**File:** `src/modules/Algorithms/src/random_walk_metropolis.jl:31`

```julia
# Computed every iteration:
Dates.value(now() - start_time) / 60000.0
```

Could cache `start_time` as milliseconds and use simpler arithmetic.

---

## Test Coverage

### Stubbed Tests Need Implementation
**File:** `tests/test_energy_calls.jl`

The following tests throw `@assert false`:
- `test_solvation_free_energy()`
- `test_cc_energy_two_molecules()`
- `test_cc_energy_multi_molecules()`
- `test_combined_potential()`
- `test_bounding_sphere_overlap()`

---

## Documentation Gaps

### Missing Module Docstring
Main `MorFit.jl` lacks a module-level docstring explaining:
- What the package does
- Main data structures (MOLECULE_DATA, EXPERIMENTAL_ASSEMBLIES)
- How to run tests

### Undocumented Functions
Priority functions needing docstrings:
1. `_find_superposition_transform` - Kabsch algorithm implementation
2. `_calculate_driven_rmsd` - Fixed target alignment logic
3. All perturbation functions

### Data Format Documentation
The `MOLECULE_DATA` and `EXPERIMENTAL_ASSEMBLIES` structures are only documented in comments. Add proper docstrings.

---

## Summary

| Priority | Count | Effort |
|----------|-------|--------|
| Medium | 3 | Medium |
| Low | 4 | Low |
| Tests | 5 | Medium |
| Docs | 8+ | Low |

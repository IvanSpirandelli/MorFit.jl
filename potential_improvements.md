# MorFit.jl Potential Improvements

## Critical Bugs

### 1. Undefined Variable in `get_moment_of_inertia`
**File:** `src/modules/Utilities/src/center_of_mass_computation.jl:79`

The function uses `num_atoms` which is never defined, causing a runtime crash.

```julia
# Current (broken):
for i in 1:num_atoms  # num_atoms undefined!

# Fix:
for i in 1:size(centered_coordinates, 2)
```

### 2. Off-by-One Indexing in `perturb_single_randomly_chosen_only_translations`
**File:** `src/modules/Utilities/src/perturbation.jl:27`

Uses 0-based indexing with Julia's 1-based arrays:

```julia
# Current (broken):
i = rand(0:(length(x)÷3)-1)  # Returns 0, 1, 2, ... but Julia is 1-indexed

# Fix:
i = rand(1:(length(x)÷3))
```

### 3. Potential Issue in `get_initial_state_only_translations`
**File:** `src/modules/Utilities/src/initialization.jl:11`

Returns flat vector instead of `Vector{Vector}`, inconsistent with other initialization functions.

---

## High Priority

### 1. Code Duplication in `add_to_output`
**File:** `src/modules/Algorithms/src/Algorithms.jl:6-22`

Three nearly identical implementations for different value types. Consolidate:

```julia
# Replace all three with:
function add_to_output(measures::Dict, output::Dict{String, Vector})
    for (key, value) in measures
        if haskey(output, key)
            push!(output[key], value)
        else
            output[key] = [value]
        end
    end
end
```

### 2. Inconsistent Timestamp Handling in `simulate!`
**File:** `src/modules/Algorithms/src/random_walk_metropolis.jl`

- Fresh simulation (line 17): Adds "timestamps" to output
- Resume simulation (line 63): Does NOT add timestamps

Both should be consistent.

### 3. Missing Documentation Throughout
No docstrings for:
- All Utilities functions (`perturbation.jl`, `realizations.jl`, `initialization.jl`)
- Algorithm structs and `simulate!` methods
- Energy function parameters and return values

---

## Medium Priority

### 1. Repeated Realization Code
**File:** `src/modules/Utilities/src/realizations.jl`

Five functions doing essentially the same thing with different output types:
- `get_flat_realization`
- `get_point_vector_realization`
- `get_point_vector_realization_per_mol`
- `get_point3f_realization`
- `get_point3f_realization_per_mol`

Consider a generic function with type parameter.

### 2. Complex One-Liners
**File:** `src/modules/Utilities/src/realizations.jl:5, 10`

Dense `hvcat` expressions are hard to read:

```julia
# Current:
[(hvcat((n_mol), [R * tc .+ t for ((R,t), tc) in zip(x, template_centers)]...)...)...]

# Clearer:
function get_flat_realization(x, template_centers)
    realized = [R * tc .+ t for ((R, t), tc) in zip(x, template_centers)]
    return vec(hcat(realized...))
end
```

### 3. Parameter Explosion in Energy Functions
**File:** `src/modules/Energies/src/combined.jl`

Functions with 13-18 parameters are hard to use correctly. Consider:
- Grouping into structs (already have `SolvationParams`, etc. in HPC_MorFit)
- Using keyword arguments with defaults

### 4. Hardcoded Default in `perturb_single_specified`
**File:** `src/modules/Utilities/src/perturbation.jl:12`

```julia
function perturb_single_specified(x, σ_r, σ_t; specified_index = 2)
```

The default `specified_index = 2` assumes protein-ligand format. Document this assumption or make it required.

### 5. Persistence Calculation Duplication
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
- `random_walk_metropolis.jl:39`
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
4. All realization functions
5. `simulate!` methods - algorithm details, output format

### Data Format Documentation
The `MOLECULE_DATA` and `EXPERIMENTAL_ASSEMBLIES` structures are only documented in comments. Add proper docstrings.

---

## Summary

| Priority | Count | Effort |
|----------|-------|--------|
| Critical | 3 | Low |
| High | 3 | Medium |
| Medium | 5 | Medium |
| Low | 4 | Low |
| Tests | 5 | Medium |
| Docs | 10+ | Low |

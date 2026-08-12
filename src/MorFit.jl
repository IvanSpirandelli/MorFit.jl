"""
    MorFit

Molecular fitting via Monte Carlo sampling of rigid-body configurations with
morphometric energy functions.

A configuration ("state") is a `Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}`
holding one (rotation, translation) pair per molecule. Energies combine solvation
free energy (morphometric approach), an overlap penalty, and persistent-homology
terms; see `MorFit.Energies`. Sampling algorithms (Random Walk Metropolis,
Simulated Annealing) live in `MorFit.Algorithms`.
"""
module MorFit

using JLD2, Rotations

include("modules/Utilities/src/Utilities.jl")
include("modules/Algorithms/src/Algorithms.jl")
include("modules/Energies/src/Energies.jl")

const _templates_dir = joinpath(@__DIR__, "templates")

@load joinpath(_templates_dir, "molecule_data.jld2") molecule_data

"""
    MOLECULE_DATA

Centered molecular building blocks (templates).

Keys are `"pdb_id:component"` strings (e.g. `"4ty7:protein"`, `"4ty7:ligand"`).
Each value is a `NamedTuple` with:
- `centers::Matrix{Float64}`: 3×N atom coordinates, centered at the origin
- `radii::Vector{Float64}`: atomic radii

Templates are positioned in space by applying a state's rotation and translation.
"""
const MOLECULE_DATA = molecule_data

@load joinpath(_templates_dir, "experimental_assemblies.jld2") experimental_assemblies

"""
    EXPERIMENTAL_ASSEMBLIES

Experimentally determined reference configurations for RMSD calculations.

Keys are `Vector{String}` of molecule IDs (e.g. `["4ty7:protein", "4ty7:ligand"]`).
Each value is a `Vector` of equivalent assemblies (accounting for permutation
symmetry); an assembly is a `NamedTuple` with:
- `centers::Vector{Matrix{Float64}}`: one 3×N matrix of *realized* coordinates
  per molecule (already positioned in space, unlike the centered templates)
- `radii::Vector{Vector{Float64}}`: atomic radii per molecule
"""
const EXPERIMENTAL_ASSEMBLIES = experimental_assemblies

#=============================================================================
# RMSD convenience functions (use MOLECULE_DATA and EXPERIMENTAL_ASSEMBLIES)
=============================================================================#

"""
    get_reference_templates(molecule_ids::Vector{String})

Get reference templates from MOLECULE_DATA for the given molecule IDs.
"""
function get_reference_templates(molecule_ids::Vector{String})
    [MOLECULE_DATA[mol_id].centers for mol_id in molecule_ids]
end

"""
    get_rmsd_to_ground_truth(molecule_ids, sim_templates, state)

Compute RMSD between simulation state and experimental ground truth.
Returns minimum RMSD across all equivalent reference assemblies.

EXPERIMENTAL_ASSEMBLIES stores realized coordinates (not states), so we:
1. Realize the simulation state to get coordinates
2. Compare with reference coordinates directly using identity states
"""
function get_rmsd_to_ground_truth(
    molecule_ids::Vector{String},
    sim_templates::Vector{Matrix{Float64}},
    state::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}
)
    if !haskey(EXPERIMENTAL_ASSEMBLIES, molecule_ids)
        @warn "molecule_ids $molecule_ids not found in EXPERIMENTAL_ASSEMBLIES"
        return Inf
    end

    # The RMSD calculation expects unrealized coordinates and rigid transformations,
    # but our reference assemblies are already realized coordinates.
    # Therefore we pass identity states.
    ref_assemblies = EXPERIMENTAL_ASSEMBLIES[molecule_ids]
    identity_state = [(one(QuatRotation), [0.0, 0.0, 0.0]) for _ in 1:length(molecule_ids)]

    min_rmsd = Inf
    for ref_assembly in ref_assemblies
        rmsd = Utilities.get_rmsd_for_fixed_target_inhibitor_pair(
            sim_templates, ref_assembly.centers, state, identity_state
        )
        min_rmsd = min(min_rmsd, rmsd)
    end

    return min_rmsd
end

_sim_templates(input) = haskey(input, "centers") ? input["centers"] : input["template_centers"]

"""
    get_min_rmsd(input, output)

Get RMSD for the minimum energy state of a simulation.
Accepts both a saved `Dict` and a `SimulationOutput` for `output`.
"""
function get_min_rmsd(input, output::Algorithms.SimulationOutput)
    molecule_ids = Vector{String}(input["molecule_ids"])
    min_energy_state = output.states[argmin(output.E_total)]
    get_rmsd_to_ground_truth(molecule_ids, _sim_templates(input), min_energy_state)
end

get_min_rmsd(input, output::AbstractDict) =
    get_min_rmsd(input, Algorithms.SimulationOutput(output))

"""
    get_min_rmsd_cutoff(input, output, cutoff_index)

Get RMSD for minimum energy state up to a given iteration cutoff.
Accepts both a saved `Dict` and a `SimulationOutput` for `output`.
"""
function get_min_rmsd_cutoff(input, output::Algorithms.SimulationOutput, cutoff_index::Int)
    molecule_ids = Vector{String}(input["molecule_ids"])

    max_index = findlast(x -> x <= cutoff_index, output.αs)
    if isnothing(max_index)
        return Inf
    end

    min_energy_state = output.states[1:max_index][argmin(output.E_total[1:max_index])]
    get_rmsd_to_ground_truth(molecule_ids, _sim_templates(input), min_energy_state)
end

get_min_rmsd_cutoff(input, output::AbstractDict, cutoff_index::Int) =
    get_min_rmsd_cutoff(input, Algorithms.SimulationOutput(output), cutoff_index)

end # module MorFit

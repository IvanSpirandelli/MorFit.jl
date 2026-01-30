module MorFit
using JLD2, Rotations

# Modules
include("modules/Utilities/src/Utilities.jl")
include("modules/Algorithms/src/Algorithms.jl")
include("modules/Energies/src/Energies.jl")

# Load template data from JLD2 files
const _templates_dir = joinpath(@__DIR__, "templates")

# MOLECULE_DATA: centered building blocks
# Keys: "pdb_id:component" (e.g., "4ty7:protein", "4ty7:ligand")
# Values: NamedTuple (centers=Matrix{Float64}, radii=Vector{Float64})
@load joinpath(_templates_dir, "molecule_data.jld2") molecule_data
const MOLECULE_DATA = molecule_data

# EXPERIMENTAL_ASSEMBLIES: reference configurations for RMSD
# Keys: Vector{String} (e.g., ["4ty7:protein", "4ty7:ligand"])
# Values: Vector of states, each state is Vector{Tuple{QuatRotation, Vector{Float64}}}
@load joinpath(_templates_dir, "experimental_assemblies.jld2") experimental_assemblies
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
        return Inf
    end

    ref_assemblies = EXPERIMENTAL_ASSEMBLIES[molecule_ids]
    identity_state = [(QuatRotation(1.0, 0.0, 0.0, 0.0), [0.0, 0.0, 0.0]) for _ in 1:length(molecule_ids)]

    min_rmsd = Inf
    for ref_assembly in ref_assemblies
        # ref_assembly.centers are already realized coordinates
        rmsd = Utilities.get_rmsd_for_fixed_target_inhibitor_pair(
            sim_templates, ref_assembly.centers, state, identity_state
        )
        min_rmsd = min(min_rmsd, rmsd)
    end

    return min_rmsd
end

"""
    get_min_rmsd(input, output)

Get RMSD for the minimum energy state of a simulation.
"""
function get_min_rmsd(input, output)
    molecule_ids = Vector{String}(input["molecule_ids"])
    sim_templates = haskey(input, "centers") ? input["centers"] : input["template_centers"]
    min_energy_state = output["states"][argmin(output["Es"])]
    get_rmsd_to_ground_truth(molecule_ids, sim_templates, min_energy_state)
end

"""
    get_min_rmsd_cutoff(input, output, cutoff_index)

Get RMSD for minimum energy state up to a given iteration cutoff.
"""
function get_min_rmsd_cutoff(input, output, cutoff_index::Int)
    molecule_ids = Vector{String}(input["molecule_ids"])
    sim_templates = haskey(input, "centers") ? input["centers"] : input["template_centers"]

    max_index = findlast(x -> x <= cutoff_index, output["αs"])
    if isnothing(max_index)
        return Inf
    end

    min_energy_state = output["states"][1:max_index][argmin(output["Es"][1:max_index])]
    get_rmsd_to_ground_truth(molecule_ids, sim_templates, min_energy_state)
end

# Tests
include("../tests/Tests.jl")

end #module MorFit
#=============================================================================
# Combined Energy Functions
#
# Computes total energy using the additive model:
#   E = θ_G · F_sol + θ_O · O + θ_T · Σ(λᵢ · Pᵢ)
#
# Routes computation based on θ values - skips components when θ = 0.
=============================================================================#

"""
    calculate_combined_energy(x, system, sol_params, ol_params, topo_params,
                              num_params, scales, precomputed, bol_check)

Calculate combined energy for a 2-molecule system using the additive model:

    E = θ_G · F_sol + θ_O · O + θ_T · T

Routes computation based on θ values - skips components when corresponding θ = 0.

# Arguments
- `x`: State vector of (rotation, translation) tuples
- `system`: MolecularSystem with templates and bounds
- `sol_params`: SolvationParams (rs, η)
- `ol_params`: LinearOverlapParams (jump, slope)
- `topo_params`: TopologyParams (λ weights)
- `num_params`: NumericalParams (delaunay tolerances)
- `scales`: EnergyScales (θ_G, θ_O, θ_T)
- `precomputed`: Precomputed single-subunit values
- `bol_check`: Bounding overlap check function

# Returns
- `(energy, measures_dict)` tuple
"""
function calculate_combined_energy(
    x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}},
    system::MolecularSystem,
    sol_params::SolvationParams,
    ol_params::LinearOverlapParams,
    topo_params::TopologyParams,
    num_params::NumericalParams,
    scales::EnergyScales,
    precomputed::Precomputed,
    bol_check::Function
)
    !in_bounds(x, system.bounds) && return Inf, Dict{String, Any}()

    energy = 0.0
    measures = Dict{String, Any}()

    # Solvation + Overlap (currently coupled in C++ lib, skip if both θ = 0)
    if scales.θ_G != 0.0 || scales.θ_O != 0.0
        prefactors = get_wb_prefactors(sol_params.rs, sol_params.η)
        fsol, overlap = solvation_free_energy_and_separate_overlap_with_bounding_container_check(
            x, system.centers, system.radii,
            sol_params.rs, prefactors, ol_params.jump, ol_params.slope,
            num_params.delaunay_eps, precomputed.single_energies, bol_check
        )
        energy += scales.θ_G * fsol + scales.θ_O * overlap
        measures["E_G"] = fsol
        measures["E_O"] = overlap
    else
        measures["E_G"] = 0.0
        measures["E_O"] = 0.0
    end

    # Topology (independent, skip if θ_T = 0)
    if scales.θ_T != 0.0
        topo_energy, topo_measures = compute_total_alpha_shape_persistence(
            x, system.centers, system.radii,
            topo_params.λ, num_params.exact_delaunay, true
        )
        energy += scales.θ_T * topo_energy
        merge!(measures, topo_measures)
        measures["E_T"] = topo_energy
    else
        measures["E_T"] = 0.0
    end

    return energy, measures
end

"""
    calculate_combined_energy(ccs, p_id, x, system, sol_params, ol_params,
                              topo_params, num_params, scales, precomputed, bol_check)

Calculate combined energy for n>2 molecule systems with connected component tracking.

Routes computation based on θ values - skips components when corresponding θ = 0.

# Additional Arguments
- `ccs`: Connected component energy cache from previous iteration
- `p_id`: Index of the perturbed molecule

# Returns
- `(energy, measures_dict, updated_ccs)` tuple
"""
function calculate_combined_energy(
    ccs::Dict{Vector{Int64}, Tuple{Float64, Dict{String, Any}}},
    p_id::Int,
    x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}},
    system::MolecularSystem,
    sol_params::SolvationParams,
    ol_params::LinearOverlapParams,
    topo_params::TopologyParams,
    num_params::NumericalParams,
    scales::EnergyScales,
    precomputed::Precomputed,
    bol_check::Function
)
    !in_bounds(x, system.bounds) && return Inf, Dict{String, Any}(), ccs

    energy = 0.0
    measures = Dict{String, Any}()
    updated_ccs = ccs

    # Solvation + Overlap (currently coupled, skip if both θ = 0)
    if scales.θ_G != 0.0 || scales.θ_O != 0.0
        prefactors = get_wb_prefactors(sol_params.rs, sol_params.η)
        fsol, fsol_measures, updated_ccs = connected_component_wise_solvation_free_energy_and_measures(
            ccs, p_id, x, system.centers, system.radii,
            sol_params.rs, prefactors, ol_params.jump, ol_params.slope,
            num_params.delaunay_eps, precomputed.single_energies, precomputed.single_measures, bol_check
        )
        overlap = get(fsol_measures, "E_O", 0.0)
        fsol_pure = fsol - overlap
        energy += scales.θ_G * fsol_pure + scales.θ_O * overlap
        merge!(measures, fsol_measures)
        measures["E_G"] = fsol_pure
        measures["E_O"] = overlap
    else
        measures["E_G"] = 0.0
        measures["E_O"] = 0.0
    end

    # Topology (independent, skip if θ_T = 0)
    if scales.θ_T != 0.0
        topo_energy, topo_measures = compute_total_alpha_shape_persistence(
            x, system.centers, system.radii,
            topo_params.λ, num_params.exact_delaunay, true
        )
        energy += scales.θ_T * topo_energy
        merge!(measures, topo_measures)
        measures["E_T"] = topo_energy
    else
        measures["E_T"] = 0.0
    end

    return energy, measures, updated_ccs
end

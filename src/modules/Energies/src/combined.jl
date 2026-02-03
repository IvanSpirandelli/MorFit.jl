#=============================================================================
# Combined Energy Functions
#
# Computes total energy using the additive model:
#   E = θ_G · F_sol + θ_O · O + θ_T · Σ(λᵢ · Pᵢ)
=============================================================================#

"""
    calculate_combined_energy(x, system, sol_params, ol_params, topo_params,
                              num_params, scales, precomputed, bol_check)

Calculate combined energy for a 2-molecule system using the additive model:

    E = θ_G · F_sol + θ_O · O + θ_T · T

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
    # Guard clause: Exit early if out of bounds
    !in_bounds(x, system.bounds) && return Inf, Dict{String, Any}()

    # Compute prefactors from solvation parameters (White Bear model)
    prefactors = get_wb_prefactors(sol_params.rs, sol_params.η)

    # 1. Compute topology energy: T = Σ(λᵢ · Pᵢ)
    topo_energy, topo_measures = compute_total_alpha_shape_persistence(
        x, system.template_centers, vcat(system.template_radii...),
        topo_params.λ, num_params.exact_delaunay, true  # compute_weighted=true
    )

    # 2. Compute solvation free energy and overlap separately
    fsol, overlap = solvation_free_energy_and_separate_overlap_with_bounding_container_check(
        x, system.template_centers, system.template_radii,
        sol_params.rs, prefactors, ol_params.jump, ol_params.slope,
        num_params.delaunay_eps, precomputed.single_energies, bol_check
    )

    # 3. Combine with scaling factors: E = θ_G · F_sol + θ_O · O + θ_T · T
    energy = scales.θ_G * fsol + scales.θ_O * overlap + scales.θ_T * topo_energy

    # 4. Build measures dictionary
    measures = merge(topo_measures, Dict{String, Any}(
        "Es_fsol" => fsol,
        "Es_overlap" => overlap,
        "Es_topo" => topo_energy,
        "OLs" => overlap
    ))

    return energy, measures
end

"""
    calculate_combined_energy(ccs, p_id, x, system, sol_params, ol_params,
                              topo_params, num_params, scales, precomputed, bol_check)

Calculate combined energy for n>2 molecule systems with connected component tracking.

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
    # Guard clause: Exit early if out of bounds
    !in_bounds(x, system.bounds) && return Inf, Dict{String, Any}(), ccs

    # Compute prefactors from solvation parameters (White Bear model)
    prefactors = get_wb_prefactors(sol_params.rs, sol_params.η)

    # 1. Compute topology energy (whole system): T = Σ(λᵢ · Pᵢ)
    radii = vcat(system.template_radii...)
    topo_energy, topo_measures = compute_total_alpha_shape_persistence(
        x, system.template_centers, radii,
        topo_params.λ, num_params.exact_delaunay, true  # compute_weighted=true
    )

    # 2. Compute solvation free energy for connected components
    # Note: CC calculation currently returns combined fsol+overlap
    # TODO: Separate overlap tracking for CC case
    fsol, fsol_measures, updated_ccs = connected_component_wise_solvation_free_energy_and_measures(
        ccs, p_id, x, system.template_centers, system.template_radii,
        sol_params.rs, prefactors, ol_params.jump, ol_params.slope,
        num_params.delaunay_eps, precomputed.single_energies, precomputed.single_measures, bol_check
    )

    # Extract overlap from measures if available
    overlap = get(fsol_measures, "OLs", 0.0)
    fsol_pure = fsol - overlap  # Remove overlap from combined value

    # 3. Combine with scaling factors: E = θ_G · F_sol + θ_O · O + θ_T · T
    energy = scales.θ_G * fsol_pure + scales.θ_O * overlap + scales.θ_T * topo_energy

    # 4. Build measures dictionary
    measures = merge(fsol_measures, topo_measures, Dict{String, Any}(
        "Es_fsol" => fsol_pure,
        "Es_overlap" => overlap,
        "Es_topo" => topo_energy
    ))

    return energy, measures, updated_ccs
end

#=============================================================================
# Legacy API (deprecated, for backward compatibility)
=============================================================================#

"""
    calculate_combined_potential(...)

DEPRECATED: Use `calculate_combined_energy` with the new parameter structs.

This function maintains backward compatibility with the old μ-interpolation model:
    E = μ · (F_sol + O) + (1-μ) · T
"""
function calculate_combined_potential(
    x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}},
    template_centers,
    radii,
    rs::Float64,
    prefactors::AbstractVector,
    overlap_jump::Float64,
    overlap_slope::Float64,
    bounds::Float64,
    persistence_weights::Vector{Float64},
    delaunay_eps::Float64,
    exact_delaunay::Bool,
    ssu_energy,
    ssu_measures,
    bol_nmol::Function,
    μ::Float64,
    compute_weighted::Bool
)
    # Guard clause: Exit early if out of bounds
    !in_bounds(x, bounds) && return Inf, Dict{String, Any}()

    # 1. Compute total alpha shape persistence
    tasp, tasp_measures = compute_total_alpha_shape_persistence(
        x, template_centers, radii, persistence_weights, exact_delaunay, compute_weighted
    )

    # 2. Compute solvation free energy
    fsol, fsol_measures = solvation_free_energy_and_measures_with_bounding_container_check(
        x, template_centers, radii, rs, prefactors,
        overlap_jump, overlap_slope, delaunay_eps, ssu_energy, ssu_measures, bol_nmol
    )

    # 3. Combine results (legacy μ-interpolation)
    energy = μ * fsol + (1 - μ) * tasp
    measures = merge!(fsol_measures, tasp_measures)
    measures = merge!(measures, Dict("Es_fsol" => fsol, "Es_topo" => tasp))
    return energy, measures
end

function calculate_combined_potential_with_separate_overlap(
    x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}},
    template_centers,
    radii,
    rs::Float64,
    prefactors::AbstractVector,
    overlap_jump::Float64,
    overlap_slope::Float64,
    bounds::Float64,
    persistence_weights::Vector{Float64},
    delaunay_eps::Float64,
    exact_delaunay::Bool,
    ssu_energy,
    bol_nmol::Function,
    μ::Float64
)
    # Guard clause: Exit early if out of bounds
    !in_bounds(x, bounds) && return Inf, Dict{String, Any}()

    # 1. Compute total alpha shape persistence
    twasp = total_weighted_alpha_shape_persistence(
        x, template_centers, radii, persistence_weights, exact_delaunay
    )

    # 2. Compute solvation free energy
    fsol, ol = solvation_free_energy_and_separate_overlap_with_bounding_container_check(
        x, template_centers, radii, rs, prefactors,
        overlap_jump, overlap_slope, delaunay_eps, ssu_energy, bol_nmol
    )

    # 3. Combine results (legacy model with separate overlap)
    energy = μ * fsol + (1 - μ) * twasp + ol
    return energy, Dict{String, Any}("Es_fsol" => fsol, "Es_topo" => twasp, "OLs" => ol)
end

function calculate_combined_potential(
    ccs,
    p_id,
    x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}},
    template_centers,
    template_radii,
    rs::Float64,
    prefactors::AbstractVector,
    overlap_jump::Float64,
    overlap_slope::Float64,
    bounds::Float64,
    persistence_weights::Vector{Float64},
    delaunay_eps::Float64,
    exact_delaunay::Bool,
    ssu_energy,
    ssu_measures,
    bol_nmol::Function,
    μ::Float64,
    compute_weighted::Bool
)
    # Guard clause: Exit early if out of bounds
    !in_bounds(x, bounds) && return Inf, Dict{String, Any}(), ccs

    # 1. Compute total alpha shape persistence
    radii = vcat(template_radii...)
    tasp, tasp_measures = compute_total_alpha_shape_persistence(
        x, template_centers, radii, persistence_weights, exact_delaunay, compute_weighted
    )

    # 2. Compute solvation free energy for the connected components
    fsol, fsol_measures, updated_ccs = connected_component_wise_solvation_free_energy_and_measures(
        ccs, p_id, x, template_centers, template_radii, rs, prefactors,
        overlap_jump, overlap_slope, delaunay_eps, ssu_energy, ssu_measures, bol_nmol
    )

    # 3. Combine results (legacy μ-interpolation)
    energy = μ * fsol + (1 - μ) * tasp
    measures = merge!(fsol_measures, tasp_measures)
    measures = merge!(measures, Dict("Es_fsol" => fsol, "Es_topo" => tasp))
    return energy, measures, updated_ccs
end

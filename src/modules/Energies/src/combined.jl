# This function is called for systems with exactly two molecules, where the full connected component construction is unnecessary
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
    tasp, tasp_measures = compute_total_alpha_shape_persistence(x, template_centers, radii, persistence_weights, exact_delaunay, compute_weighted)

    # 2. Compute solvation free energy
    fsol, fsol_measures = solvation_free_energy_and_measures_with_bounding_container_check(
        x, template_centers, radii, rs, prefactors, 
        overlap_jump, overlap_slope, delaunay_eps, ssu_energy, ssu_measures, bol_nmol
    )

    # 3. Combine results and return
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
    twasp = total_weighted_alpha_shape_persistence(x, template_centers, radii, persistence_weights, exact_delaunay)

    # 2. Compute solvation free energy
    fsol, ol = solvation_free_energy_and_separate_overlap_with_bounding_container_check(
        x, template_centers, radii, rs, prefactors, 
        overlap_jump, overlap_slope, delaunay_eps, ssu_energy, bol_nmol
    )
    # 3. Combine results and return
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
    radii = vcat([template_radii for _ in 1:length(x)]...)
    tasp, tasp_measures = compute_total_alpha_shape_persistence(x, template_centers, radii, persistence_weights, exact_delaunay, compute_weighted)

    # 2. Compute solvation free energy for the connected components
    fsol, fsol_measures, updated_ccs = connected_component_wise_solvation_free_energy_and_measures(
        ccs, p_id, x, template_centers, template_radii, rs, prefactors,
        overlap_jump, overlap_slope, delaunay_eps, ssu_energy, ssu_measures, bol_nmol
    )

    # 3. Combine results and return
    energy = μ * fsol + (1 - μ) * tasp
    measures = merge!(fsol_measures, tasp_measures)
    measures = merge!(measures, Dict("Es_fsol" => fsol, "Es_topo" => tasp))
    return energy, measures, updated_ccs
end
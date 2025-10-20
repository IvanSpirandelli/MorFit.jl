function get_initialization(input, non_overlapping = false)
    n_mol = input["n_mol"]
    bounds = input["bounds"]
    if input["initialization"] == "random"
        if non_overlapping  
            return () -> get_non_overlapping_initial_state(input)
        else
            return () -> Utilities.get_initial_state(n_mol, bounds)
        end
    elseif input["initialization"] == "random_only_translations"
        return () -> Utilities.get_initial_state_only_translations(n_mol, bounds)
    end
end

function get_non_overlapping_initial_state(input, max_count = 25)
    n_mol = input["n_mol"]
    bounds = input["bounds"]
    x = Utilities.get_initial_state(n_mol, bounds)
    attempt = 0
    while _is_overlapping_state(x, input)
        x = Utilities.get_initial_state(n_mol, bounds)
        attempt += 1
        if attempt >= max_count
            println("Maximum attempts to find non overlapping initial state reached.")
            @assert false 
        end
    end
    return x
end    

function _get_realization_radii_and_sizes(x, input)
    assembly_type = split(input["assembly_spec"],"_")[1]
    if assembly_type == "sta"
        return Utilities.get_flat_realization(x, input["template_centers"]), fill(length(input["template_radii"]), input["n_mol"]), vcat([input["template_radii"] for i in 1:input["n_mol"]]...)
    elseif assembly_type == "ppii"
        return Utilities.get_flat_realization(x, input["template_centers"]), [length(r) for r in input["template_radii"]], vcat(input["template_radii"]...)
    end
end 
function _is_overlapping_state(x, input)
    realz, sizes, radii = _get_realization_radii_and_sizes(x, input)
    measures = Energies.get_geometric_measures_and_overlap_value(
        realz, 
        sizes, 
        radii, 
        0.0, 
        1.0,
        0.0, 
        100.0
    )
    return !isapprox(0.0, measures[5])
end

function get_perturbation(input)
    σ_r = input["σ_r"]
    σ_t = input["σ_t"]
    if input["perturbation"] == "single_random"
        return (x) -> Utilities.perturb_single_randomly_chosen(x, σ_r, σ_t)
    elseif input["perturbation"] == "single_specified"
        return (x) -> Utilities.perturb_single_specified(x, σ_r, σ_t)
    elseif input["perturbation"] == "single_random_get_index"
        return (x) -> Utilities.get_index_and_perturb_single_randomly_chosen(x, σ_r, σ_t)    
    elseif input["perturbation"] == "single_random_only_translations"
        return (x) -> Utilities.perturb_single_randomly_chosen_only_translations(x, σ_t)
    elseif input["perturbation"] == "all"
        return (x) -> Utilities.perturb_all(x, σ_r, σ_t)
    elseif input["perturbation"] == "hs_all"
        return (x) -> Utilities.perturb_all(x, σ_t, σ_t)
    end
end

function get_energy(input)
    if input["energy"] == "tasp"
        return (x) -> Energies.total_alpha_shape_persistence_and_measures(x, input["template_centers"], input["persistence_weights"], input["exact_delaunay"])
    elseif input["energy"] == "twasp"
        radii = vcat([input["template_radii"] for _ in 1:input["n_mol"]]...)
        return (x) -> Energies.total_weighted_alpha_shape_persistence_and_measures(x, input["template_centers"], radii, input["persistence_weights"], input["exact_delaunay"])
    elseif input["energy"] == "hstasp"
        return (x) -> Energies.hs_total_alpha_shape_persistence(x, input["persistence_weights"], input["exact_delaunay"])
    elseif input["energy"] == "fsol"
        radii = vcat([input["template_radii"] for _ in 1:input["n_mol"]]...)
        return (x) -> Energies.solvation_free_energy_and_measures_in_bounds(x, input["template_centers"], radii, input["rs"], input["prefactors"], input["overlap_jump"], input["overlap_slope"], input["bounds"], input["delaunay_eps"])
    elseif input["energy"] == "cc_fsol"
        return get_connected_component_solvation_free_energy_in_bounds_energy_call(input)
    elseif input["energy"] == "cc_fsol_twasp_interpolated"
        return get_connected_component_solvation_free_energy_with_total_alpha_shape_persistence_interpolated_in_bounds_energy_call(input, true)
    elseif input["energy"] == "cc_fsol_twasp_interpolated_ol"
        return get_connected_component_solvation_free_energy_with_total_alpha_shape_persistence_interpolated_in_bounds_seperate_overlap_energy_call(input, true)
    else    
        return (x) -> 0.0
    end
end

function get_connected_component_solvation_free_energy_with_total_alpha_shape_persistence_interpolated_in_bounds_energy_call(input, weighted = true)
    rs = input["rs"]
    prefactors = input["prefactors"]
    overlap_jump = input["overlap_jump"]
    overlap_slope = input["overlap_slope"]
    delaunay_eps = input["delaunay_eps"]
    template_centers = input["template_centers"]
    template_radii = input["template_radii"]
    bounds = input["bounds"]
    persistence_weights = input["persistence_weights"]
    exact_delaunay = input["exact_delaunay"]
    μ = input["μ"]

    ssu_energy, ssu_measures = Energies.get_single_subunit_energy_and_measures(template_centers, template_radii, rs, prefactors, overlap_jump, overlap_slope, delaunay_eps)

    if input["n_mol"] == 2
        radii = vcat([input["template_radii"] for _ in 1:input["n_mol"]]...)
        bol_nmol = (x) -> Energies.are_bounding_spheres_overlapping(x, 1, 2, Energies.get_bounding_radii(template_centers, template_radii, rs))
        return (x) -> Energies.calculate_combined_potential(
            x, template_centers, radii, 
            rs, prefactors, overlap_jump, overlap_slope, bounds, persistence_weights, 
            delaunay_eps, exact_delaunay, 
            ssu_energy, ssu_measures, bol_nmol, 
            μ, weighted
        )
    else
        bol_nmol = (x, id1, id2) -> Energies.are_bounding_spheres_overlapping(x, id1, id2, Energies.get_bounding_radii(template_centers, template_radii, rs))
        return (ccs, p_id, x) -> Energies.calculate_combined_potential(
            ccs, p_id, 
            x, template_centers, template_radii, 
            rs, prefactors, overlap_jump, overlap_slope, bounds, persistence_weights, 
            delaunay_eps, exact_delaunay, 
            ssu_energy, ssu_measures, bol_nmol, 
            μ, weighted
        )
    end
end

function get_connected_component_solvation_free_energy_with_total_alpha_shape_persistence_interpolated_in_bounds_seperate_overlap_energy_call(input, weighted = true)
    rs = input["rs"]
    prefactors = input["prefactors"]
    overlap_jump = input["overlap_jump"]
    overlap_slope = input["overlap_slope"]
    delaunay_eps = input["delaunay_eps"]
    template_centers = input["template_centers"]
    template_radii = input["template_radii"]
    bounds = input["bounds"]
    persistence_weights = input["persistence_weights"]
    exact_delaunay = input["exact_delaunay"]
    μ = input["μ"]

    ssu_energy, _ = Energies.get_single_subunit_energy_and_measures(template_centers, template_radii, rs, prefactors, overlap_jump, overlap_slope, delaunay_eps)
    radii = vcat([input["template_radii"] for _ in 1:input["n_mol"]]...)
    bol_nmol = (x) -> Energies.are_bounding_spheres_overlapping(x, 1, 2, Energies.get_bounding_radii(template_centers, template_radii, rs))
    return (x) -> Energies.calculate_combined_potential_with_separate_overlap(
        x, template_centers, radii, 
        rs, prefactors, overlap_jump, overlap_slope, bounds, persistence_weights, 
        delaunay_eps, exact_delaunay, 
        ssu_energy, bol_nmol, 
        μ
    )
end

function get_connected_component_solvation_free_energy_in_bounds_energy_call(input)
    rs = input["rs"]
    prefactors = input["prefactors"]
    overlap_jump = input["overlap_jump"]
    overlap_slope = input["overlap_slope"]
    delaunay_eps = input["delaunay_eps"]
    template_centers = input["template_centers"]
    template_radii = input["template_radii"]
    bounds = input["bounds"]

    ssu_energy, ssu_measures = Energies.get_single_subunit_energy_and_measures(template_centers, template_radii, rs, prefactors, overlap_jump, overlap_slope, delaunay_eps)

    if input["n_mol"] == 2
        radii = vcat([input["template_radii"] for _ in 1:input["n_mol"]]...)
        bol_nmol = (x) -> Energies.are_bounding_spheres_overlapping(x, 1, 2, Energies.get_bounding_radii(template_centers, template_radii, rs))
        return (x) -> Energies.solvation_free_energy_and_measures__with_bounding_container_check_in_bounds(
            x, template_centers, radii, 
            rs, prefactors, overlap_jump, overlap_slope, bounds, 
            delaunay_eps, ssu_energy, ssu_measures, bol_nmol
        )
    else
        bol_nmol = (x, id1, id2) -> Energies.are_bounding_spheres_overlapping(x, id1, id2, get_bounding_radii(template_centers, template_radii, rs))
        return (ccs, p_id, x) -> Energies.connected_component_wise_solvation_free_energy_and_measures_in_bounds(
            ccs, p_id, 
            x, template_centers, template_radii, 
            rs, prefactors, overlap_jump, overlap_slope, bounds, 
            delaunay_eps, ssu_energy, ssu_measures, bol_nmol
        )
    end
end
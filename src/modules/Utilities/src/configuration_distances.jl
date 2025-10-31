using Combinatorics
using LinearAlgebra
using Statistics

function get_metric_to_ground_truth(metric_func::Function, n_mol::Int, mol_type::String, sim_templates, state; metric_kwargs...)
    if !(mol_type in keys(EXPERIMENTAL_ASSEMBLY))
        return Inf
    end
    
    sim_state = typeof(state) == typeof(Vector{Float64}()) ? convert_flat_state_to_tuples(state) : state

    # Prepare reference states from experimental data
    if n_mol == 2
        ref_templates = EXPERIMENTAL_ASSEMBLY[mol_type]["template_centers"]
        flat_exp = EXPERIMENTAL_ASSEMBLY[mol_type]["state"]
        ref_state = convert_flat_state_to_tuples(flat_exp)
        return get_configuration_distance(sim_state, ref_state, sim_templates, ref_templates, metric_func; metric_kwargs...)
    
    elseif n_mol == 3 && mol_type == "6r7m" || mol_type == "1aly"
        ref_templates = EXPERIMENTAL_ASSEMBLY[mol_type]["template_centers"]
        exp_states_flat = get_assembled_three_states(mol_type) 
        
        min_metric_val = Inf
        for flat in exp_states_flat
            ref_state = convert_flat_state_to_tuples(flat)
            # Call the master distance function for each possible ground truth state
            cand_metric = get_configuration_distance(sim_state, ref_state, sim_templates, ref_templates, metric_func; metric_kwargs...)
            min_metric_val = min(min_metric_val, cand_metric)
        end
        return min_metric_val
    end
end

get_theta_to_ground_truth(n_mol::Int, mol_type::String, sim_templates, state) = get_metric_to_ground_truth(get_theta_of_pair, n_mol, mol_type, sim_templates, state)
get_screw_axis_distance_to_ground_truth(n_mol::Int, mol_type::String, sim_templates, state; λ_dist=0.1) = get_metric_to_ground_truth(get_screw_axis_distance_of_pair, n_mol, mol_type, sim_templates, state; λ_dist=λ_dist)
get_rmsd_align_one_to_ground_truth(n_mol::Int, mol_type::String, sim_templates, state) = get_metric_to_ground_truth(get_rmsd_align_one_of_pair, n_mol, mol_type, sim_templates, state)

"""
    get_min_metric(metric_func, input, output; metric_kwargs...)

A generic top-level function to get the specified metric for the minimum energy state
of a simulation, comparing it against the experimental ground truth.
"""
function get_min_metric(metric_func::Function, input::Dict{String, Any}, output::Dict{String, Vector}; metric_kwargs...)
    mol_type = input["mol_type"]
    n_mol = input["n_mol"]
    sim_templates = input["template_centers"] # Assuming Matrix{Float64}
    
    # Get the minimum energy state from the simulation
    min_energy_state = output["states"][argmin(output["Es"])]
    min_energy_state = typeof(min_energy_state) == typeof(Vector{Float64}()) ? convert_flat_state_to_tuples(min_energy_state) : min_energy_state

    get_metric_to_ground_truth(metric_func, n_mol, mol_type, sim_templates, min_energy_state; metric_kwargs...)
end

get_min_theta(input, output) = get_min_metric(get_theta_of_pair, input, output)
get_min_screw_axis_distance(input, output; λ_dist=0.1) = get_min_metric(get_screw_axis_distance_of_pair, input, output; λ_dist=λ_dist)
get_min_rmsd_align_one(input, output) = get_min_metric(get_rmsd_align_one_of_pair, input, output)

function get_min_metric_cut_off(metric_func::Function, input::Dict{String, Any}, output::Dict{String, Vector}, cutoff_index::Int; metric_kwargs...)
    mol_type = input["mol_type"]
    n_mol = input["n_mol"]
    sim_templates = input["template_centers"] # Assuming Matrix{Float64}

    max_index = findlast(x -> x <= cutoff_index, output["αs"])

    # Get the minimum energy state from the simulation
    min_energy_state = output["states"][1:max_index][argmin(output["Es"][1:max_index])]
    min_energy_state = typeof(min_energy_state) == typeof(Vector{Float64}()) ? convert_flat_state_to_tuples(min_energy_state) : min_energy_state

    get_metric_to_ground_truth(metric_func, n_mol, mol_type, sim_templates, min_energy_state; metric_kwargs...)
end

get_min_theta_cut_off(input, output, cutoff_index) = get_min_metric_cut_off(get_theta_of_pair, input, output, cutoff_index)
get_min_screw_axis_distance_cut_off(input, output, cutoff_index; λ_dist=0.1) = get_min_metric_cut_off(get_screw_axis_distance_of_pair, input, output, cutoff_index; λ_dist=λ_dist)
get_min_rmsd_align_one_cut_off(input, output, cutoff_index) = get_min_metric_cut_off(get_rmsd_align_one_of_pair, input, output, cutoff_index)


"""
Calculates the distance between two multi-molecule configurations using a specified pairwise metric.
It finds the optimal permutation of molecules in `state_a` to match `state_b` by
exhaustively checking all permutations and minimizing the sum of pairwise metric values.
"""
function get_configuration_distance(
    state_a::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}},
    state_b::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}},
    templates_a::Matrix{Float64},
    templates_b::Matrix{Float64},
    pair_metric::Function;
    metric_kwargs... # To pass things like λ_dist
)
    n_mol = length(state_a)
    @assert n_mol == length(state_b) "States must have the same number of molecules."
    if n_mol < 2
        return 0.0
    end
    
    identity_permutation = 1:n_mol
    pairings = collect(combinations(identity_permutation, 2))
    
    if isempty(pairings)
        return 0.0
    end

    min_sum = Inf
    # Find best permutation of state_a to match state_b by fixing state_b's permutation
    fixed_perm_b = collect(identity_permutation)
    for perm_a in permutations(collect(identity_permutation))
        # Calculate the sum of the metric over all pairs for the current permutation
        current_sum = sum(
            mean( # Note: `mean` assumes the metric returns a vector of per-atom errors.
                pair_metric(
                    templates_a, templates_b,
                    state_a, state_b,
                    [perm_a[p[1]], perm_a[p[2]]], # Get molecules from permuted state_a
                    [fixed_perm_b[p[1]], fixed_perm_b[p[2]]], # Get molecules from fixed state_b
                    ;metric_kwargs...
                )
            ) for p in pairings
        )
        min_sum = min(min_sum, current_sum)
    end

    return min_sum / length(pairings)
end

function get_theta_of_pair(
    template_centers_a::Matrix{Float64}, template_centers_b::Matrix{Float64},
    state_a::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}},
    state_b::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}},
    subset_a::Vector{Int}, subset_b::Vector{Int}
)
    # This nested helper correctly uses a specific template with its corresponding state
    function transform_dist(template_centers, state, ss, idx)
        i, j = ss[1], ss[2]
        R1, T1 = state[i]
        R2, T2 = state[j]
        return euclidean(T1 + R1 * template_centers[:,idx], T2 + R2 * template_centers[:,idx])
    end
    
    return [abs(transform_dist(template_centers_a, state_a, subset_a, i) - transform_dist(template_centers_b, state_b, subset_b, i)) for i in 1:size(template_centers_a, 2)]
end

"""
Calculates the distance between two molecular pairs based on their relative screw motion.

The process is:
1.  Calculate the relative transformation (rotation and translation) for the pair in `state_a`.
2.  Calculate the relative transformation for the corresponding pair in `state_b`.
3.  Decompose each transformation into its screw axis parameters: rotation angle `θ` and
    translation along the axis `d`.
4.  Compute the final distance as the Euclidean distance in the weighted parameter space:
    `sqrt((θ_sim - θ_ref)^2 + (λ_dist * (d_sim - d_ref))^2)`.
    The `λ_dist` parameter weights the importance of translational vs. rotational differences.
"""
function get_screw_axis_distance_of_pair(
    templates_a::Matrix{Float64}, templates_b::Matrix{Float64}, # Unused, for signature compatibility
    state_a::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}},
    state_b::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}},
    subset_a::Vector{Int}, subset_b::Vector{Int};
    λ_dist::Float64 = 0.1
)
    # --- Helper to get relative motion ---
    function get_diff_motion(state)
        R1, t1 = state[1]; R2, t2 = state[2]
        return R1' * R2, R1' * (t2 - t1)
    end
    
    # --- helper to calculate screw parameters from a single relative motion ---
    function get_screw_params(R_diff, t_diff)
        # Angle of rotation
        # Clamp argument to acos to prevent domain errors from floating point inaccuracies
        θ = acos(clamp((tr(R_diff) - 1) / 2, -1.0, 1.0))

        # Handle the pure translation case (no rotation)
        if isapprox(θ, 0.0, atol=1e-8)
            # For pure translation, the "axis" is undefined, and the translation
            # along the axis is simply the magnitude of the translation vector.
            return 0.0, norm(t_diff)
        end

        # For rotation, find the screw axis and the translation along it
        axis = [R_diff[3, 2] - R_diff[2, 3], R_diff[1, 3] - R_diff[3, 1], R_diff[2, 1] - R_diff[1, 2]] / (2 * sin(θ))
        d = dot(axis, t_diff) # Translation component parallel to the axis
        
        return θ, d
    end

    # --- Main logic ---

    # 1. Use the subsets to extract the specific pair from the full state
    sim_pair_state = [state_a[subset_a[1]], state_a[subset_a[2]]]
    ref_pair_state = [state_b[subset_b[1]], state_b[subset_b[2]]]

    # 2. Get the relative motion for both the simulation and reference pairs
    R_sim, t_sim = get_diff_motion(sim_pair_state)
    R_ref, t_ref = get_diff_motion(ref_pair_state)
    
    # 3. Calculate the screw parameters (θ, d) for each motion
    θ_sim, d_sim = get_screw_params(R_sim, t_sim)
    θ_ref, d_ref = get_screw_params(R_ref, t_ref)

    # 4. Calculate the weighted Euclidean distance between the (θ, d) parameter pairs
    #    This is the core fix: actually comparing sim vs. ref.
    distance = sqrt((θ_sim - θ_ref)^2 + (λ_dist * (d_sim - d_ref))^2)

    # Return as a single-element vector to work with the `mean` call in the master function
    return [distance]
end

function _find_superposition_transform(mobile, ref)
    c_mob = mean(mobile); c_ref = mean(ref)
    X = hcat([p .- c_mob for p in mobile]...); Y = hcat([p .- c_ref for p in ref]...)
    H = X * Y'; svd_res = svd(H)
    d = sign(det(svd_res.V * svd_res.U'))
    R = svd_res.V * diagm([1.0, 1.0, d]) * svd_res.U'
    t = c_ref - R * c_mob
    return R, t
end

function _calculate_driven_rmsd(driver_mob, follower_mob, driver_ref, follower_ref)
    R, t = _find_superposition_transform(driver_mob, driver_ref)

    k = length(follower_mob)
    if k == 0
        return 0.0
    end
    
    sq_err = sum(norm((R * follower_mob[i] + t) - follower_ref[i])^2 for i in 1:k)

    return sqrt(sq_err / k)
end

function get_rmsd_align_one_of_pair(
    templates_a::Matrix{Float64}, templates_b::Matrix{Float64},
    state_a::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}},
    state_b::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}},
    subset_a::Vector{Int}, subset_b::Vector{Int};
)
    sim_mol_1 = state_a[subset_a[1]]
    sim_mol_2 = state_a[subset_a[2]]
    ref_mol_1 = state_b[subset_b[1]]
    ref_mol_2 = state_b[subset_b[2]]

    coords_W = get_point_vector_realization([sim_mol_1], templates_a)
    coords_X = get_point_vector_realization([sim_mol_2], templates_a)
    coords_Y = get_point_vector_realization([ref_mol_1], templates_b)
    coords_Z = get_point_vector_realization([ref_mol_2], templates_b)

    rmsd_A = _calculate_driven_rmsd(coords_W, coords_X, coords_Y, coords_Z)
    rmsd_B = _calculate_driven_rmsd(coords_W, coords_X, coords_Z, coords_Y)

    return [min(rmsd_A, rmsd_B)]
end



function get_global_rmsd_of_target_inhibitor_pair(
    templates_sim::Vector{Matrix{Float64}}, templates_ref::Vector{Matrix{Float64}},
    state_sim::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}},
    state_ref::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}
)
    #@assert false "Still needs to be implemented. Handleing the Vector data input"

    sim_target_state = state_sim[1]
    sim_inhib_state = state_sim[2]
    ref_target_state = state_ref[1]
    ref_inhib_state = state_ref[2]

    coords_W = get_point_vector_realization([sim_target_state], templates_sim[1])
    coords_X = get_point_vector_realization([sim_inhib_state], templates_sim[2])
    coords_Y = get_point_vector_realization([ref_target_state], templates_ref[1])
    coords_Z = get_point_vector_realization([ref_inhib_state], templates_ref[2])

    mobile = [coords_W; coords_X]
    reference = [coords_Y; coords_Z]

    (R,t) = _find_superposition_transform(mobile, reference)

    k = length(mobile)
    sq_err = sum(norm((R * mobile[i] + t) - reference[i])^2 for i in 1:k)
    return sqrt(sq_err / k)
end


function get_rmsd_for_fixed_target_inhibitor_pair(
    templates_sim::Vector{Matrix{Float64}}, templates_ref::Vector{Matrix{Float64}},
    state_sim::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}},
    state_ref::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}
)
    #@assert false "Still needs to be implemented. Handleing the Vector data input"

    sim_target_state = state_sim[1]
    sim_inhib_state = state_sim[2]
    ref_target_state = state_ref[1]
    ref_inhib_state = state_ref[2]

    coords_W = get_point_vector_realization([sim_target_state], templates_sim[1])
    coords_X = get_point_vector_realization([sim_inhib_state], templates_sim[2])
    coords_Y = get_point_vector_realization([ref_target_state], templates_ref[1])
    coords_Z = get_point_vector_realization([ref_inhib_state], templates_ref[2])

    return _calculate_driven_rmsd(coords_W, coords_X, coords_Y, coords_Z)
end

# These methods simply convert the flat vectors to tuples

function get_theta_of_pair(
    template_centers_a::Matrix{Float64}, template_centers_b::Matrix{Float64},
    state_a::Vector{Float64}, state_b::Vector{Float64},
    subset_a::Vector{Int}, subset_b::Vector{Int}
)
    state_a_tuples = convert_flat_state_to_tuples(state_a)
    state_b_tuples = convert_flat_state_to_tuples(state_b)
    return get_theta_of_pair(template_centers_a, template_centers_b, state_a_tuples, state_b_tuples, subset_a, subset_b)
end
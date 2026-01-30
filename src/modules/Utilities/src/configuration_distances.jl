using LinearAlgebra
using Statistics

#=============================================================================
# Core RMSD calculations
=============================================================================#

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

"""
    get_rmsd_for_fixed_target_inhibitor_pair(templates_sim, templates_ref, state_sim, state_ref)

Compute RMSD for a target-inhibitor pair by aligning on the target (molecule 1)
and measuring the ligand (molecule 2) displacement.
"""
function get_rmsd_for_fixed_target_inhibitor_pair(
    templates_sim::Vector{Matrix{Float64}}, templates_ref::Vector{Matrix{Float64}},
    state_sim::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}},
    state_ref::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}
)
    sim_target_state = state_sim[1]
    sim_ligand_state = state_sim[2]
    ref_target_state = state_ref[1]
    ref_ligand_state = state_ref[2]

    coords_W = get_point_vector_realization([sim_target_state], [templates_sim[1]])
    coords_X = get_point_vector_realization([sim_ligand_state], [templates_sim[2]])
    coords_Y = get_point_vector_realization([ref_target_state], [templates_ref[1]])
    coords_Z = get_point_vector_realization([ref_ligand_state], [templates_ref[2]])

    return _calculate_driven_rmsd(coords_W, coords_X, coords_Y, coords_Z)
end
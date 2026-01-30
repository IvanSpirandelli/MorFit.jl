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
"""
    random_rotation(σ_r) -> QuatRotation{Float64}

Random rotation whose rotation-vector components are drawn from `N(0, σ_r²)`.
"""
random_rotation(σ_r) = QuatRotation(exp(Rotations.RotationVecGenerator((randn(3) .* σ_r)...)))

_perturbed_entry((R, t), σ_r, σ_t) = (R * random_rotation(σ_r), t .+ randn(3) .* σ_t)

function perturb_all(x, σ_r, σ_t)
    [_perturbed_entry(e, σ_r, σ_t) for e in x]
end

function perturb_single_specified(x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}, σ_r, σ_t; specified_index = 2)
    x_cand = copy(x)
    x_cand[specified_index] = _perturbed_entry(x[specified_index], σ_r, σ_t)
    x_cand
end

function get_index_and_perturb_single_randomly_chosen(x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}, σ_r, σ_t)
    i = rand(1:length(x))
    i, perturb_single_specified(x, σ_r, σ_t; specified_index = i)
end

function perturb_single_randomly_chosen(x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}, σ_r, σ_t)
    get_index_and_perturb_single_randomly_chosen(x, σ_r, σ_t)[2]
end

function perturb_single_point(x::Vector{Vector{Float64}}, σ_t::Float64)
    x_cand = copy(x)
    i = rand(1:length(x))
    x_cand[i] = x_cand[i] .+ randn(3) .* σ_t
    return x_cand
end

function perturb_all(x, σ_r, σ_t)
    [(R * exp(Rotations.RotationVecGenerator(randn(3) .* σ_r...)), t .+ (randn(3) .* σ_t)) for (R, t) in x]
end

function perturb_single_randomly_chosen(x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}, σ_r, σ_t)
    x_cand = deepcopy(x)
    i  = rand(1:length(x))
    x_cand[i] = (x_cand[i][1] * exp(Rotations.RotationVecGenerator(randn(3) .* σ_r...)), x_cand[i][2] .+ (randn(3) .* σ_t))
    x_cand
end

function perturb_single_specified(x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}, σ_r, σ_t; specified_index = 2)
    x_cand = deepcopy(x)
    x_cand[specified_index] = (x_cand[specified_index][1] * exp(Rotations.RotationVecGenerator(randn(3) .* σ_r...)), x_cand[specified_index][2] .+ (randn(3) .* σ_t))
    x_cand
end

function get_index_and_perturb_single_randomly_chosen(x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}, σ_r, σ_t)
    x_cand = deepcopy(x)
    i  = rand(1:length(x))
    x_cand[i] = (x_cand[i][1] * exp(Rotations.RotationVecGenerator(randn(3) .* σ_r...)), x_cand[i][2] .+ (randn(3) .* σ_t))
    i, x_cand
end

function perturb_single_randomly_chosen_only_translations(x, σ_t)
    x_cand = deepcopy(x)
    i  = rand(0:(length(x)÷3)-1)
    x_cand[(i*3)+1:(i*3)+3] = x_cand[(i*3)+1:(i*3)+3] .+ (randn(3) .* σ_t)
    x_cand
end

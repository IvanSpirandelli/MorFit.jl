function convert_flat_state_to_tuples(state::Vector{Float64})
    mod(length(state), 6) == 0 || throw(ArgumentError("flat state vector length must be a multiple of 6, got $(length(state))"))
    n_mol = length(state) ÷ 6
    state_tuples = Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}(undef, n_mol)
    for i in 1:n_mol
        idx_start = (i - 1) * 6 + 1
        R = QuatRotation(exp(Rotations.RotationVecGenerator(state[idx_start:idx_start+2]...)))
        T = state[idx_start+3:idx_start+5]
        state_tuples[i] = (R, T)
    end
    return state_tuples
end

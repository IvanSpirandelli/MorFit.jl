using Rotations

function get_flat_realization(x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}, template_centers::Vector{Matrix{Float64}})
    n_mol = length(x)
    [(hvcat((n_mol), [R * tc .+ t for ((R,t), tc) in zip(x, template_centers)]...)...)...]
end

function get_point_vector_realization(x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}, template_centers::Vector{Matrix{Float64}})
    n_mol = length(x)
    [Vector{Float64}(e) for e in eachcol(hvcat((n_mol), [R * tc .+ t for ((R,t), tc) in zip(x, template_centers)]...))]
end

function get_point_vector_realization_per_mol(x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}, template_centers::Vector{Matrix{Float64}})
    [[Vector{Float64}(e) for e in eachcol(R * tc .+ t)] for ((R, t), tc) in zip(x, template_centers)]
end

function get_point3f_realization(x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}, template_centers::Vector{Matrix{Float64}})
    n_mol = length(x)
    [Point3f(e) for e in eachcol(hvcat((n_mol), [R * tc .+ t for ((R,t), tc) in zip(x, template_centers)]...))]
end

function get_point3f_realization_per_mol(x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}, template_centers::Vector{Matrix{Float64}})
    [[Point3f(e) for e in eachcol(R * tc .+ t)] for ((R, t), tc) in zip(x, template_centers)]
end

function get_matrix_realization_per_mol(x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}, template_centers::Vector{Matrix{Float64}})
    [R * tc .+ t for ((R, t), tc) in zip(x, template_centers)]
end

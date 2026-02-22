function in_bounds(x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}, bounds::Float64)
    all(all(0.0 <= e <= bounds for e in t) for (_,t) in x)
end

function in_bounds(x::Vector{Vector{Float64}}, bounds::Float64)
    all(all(0.0 <= c <= bounds for c in point) for point in x)
end

function in_bounds_sphere(x::Vector{Vector{Float64}}, radius::Float64)
    all(sum(point.^2) <= radius^2 for point in x)
end

function _prepare_args(centers::Vector{Matrix{Float64}}, radii::Vector{Vector{Float64}})
    n_atoms_per_mol = [size(tc, 2) for tc in centers]
    flat_radii = vcat(radii...)
    return n_atoms_per_mol, flat_radii
end

function _format_measures_output(measures::AbstractVector, prefactors::AbstractVector)
    energy = sum(measures .* [prefactors; 1.0])
    measures_dict = Dict{String, Any}(
        "Vs" => measures[1], "As" => measures[2], "Cs" => measures[3],
        "Xs" => measures[4], "E_O" => measures[5]
    )
    return energy, measures_dict
end

function solvation_free_energy_and_measures(
    x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}},
    centers::Vector{Matrix{Float64}},
    radii::Vector{Vector{Float64}},
    rs::Float64,
    prefactors::AbstractVector,
    overlap_jump::Float64,
    overlap_slope::Float64,
    delaunay_eps::Float64
    )
    n_atoms_per_mol, flat_radii = _prepare_args(centers, radii)

    flat_realization = Utilities.get_realization(x, centers, format=:flat)

    measures = get_geometric_measures_and_overlap_value(
        flat_realization,
        n_atoms_per_mol,
        flat_radii,
        rs,
        overlap_jump,
        overlap_slope,
        delaunay_eps
    )

    return _format_measures_output(measures, prefactors)
end

function solvation_free_energy_and_overlap(
    x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}},
    centers::Vector{Matrix{Float64}},
    radii::Vector{Vector{Float64}},
    rs::Float64,
    prefactors::AbstractVector,
    overlap_jump::Float64,
    overlap_slope::Float64,
    delaunay_eps::Float64
    )
    n_atoms_per_mol, flat_radii = _prepare_args(centers, radii)

    flat_realization = Utilities.get_realization(x, centers, format=:flat)

    measures = get_geometric_measures_and_overlap_value(
        flat_realization,
        n_atoms_per_mol,
        flat_radii,
        rs,
        overlap_jump,
        overlap_slope,
        delaunay_eps
    )
    fsol = sum(measures[1:4] .* prefactors)
    ol = measures[5]
    return fsol, ol
end

function solvation_free_energy(
    atom_coordinates::Vector{Float64},
    atom_radii::Vector{Float64},
    probe_radius::Float64,
    prefactors::AbstractVector,
    delaunay_eps::Float64 = 1.0)
    measures = get_geometric_measures(
        atom_coordinates,
        atom_radii,
        probe_radius,
        delaunay_eps
    )
    sum(measures .* prefactors)
end

function solvation_free_energy(
    atom_coordinates::Vector{Float64},
    molecule_sizes::Vector{Int},
    atom_radii::Vector{Float64},
    probe_radius::Float64,
    prefactors::AbstractVector,
    overlap_existence_penalty::Float64,
    overlap_penalty_slope::Float64,
    delaunay_eps::Float64 = 1.0)
    measures = get_geometric_measures_and_overlap_value(
        atom_coordinates,
        molecule_sizes,
        atom_radii,
        probe_radius,
        overlap_existence_penalty,
        overlap_penalty_slope,
        delaunay_eps
    )
    sum(measures .* [prefactors; 1.0])
end

function in_bounds(x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}, bounds::Float64)
    all(all(0.0 <= e <= bounds for e in t) for (_,t) in x)
end

function _prepare_args(template_centers::Vector{Matrix{Float64}}, radii::Vector{Vector{Float64}})
    n_atoms_per_mol = [size(tc, 2) for tc in template_centers]
    flat_radii = vcat(radii...)
    return n_atoms_per_mol, flat_radii
end

function solvation_free_energy(
    x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}},
    template_centers::Vector{Matrix{Float64}},
    radii::Vector{Vector{Float64}},
    rs::Float64,
    prefactors::AbstractVector,
    overlap_jump::Float64,
    overlap_slope::Float64,
    delaunay_eps::Float64
    )
    n_atoms_per_mol, flat_radii = _prepare_args(template_centers, radii)
    atom_coordinates = Utilities.get_realization(x, template_centers, format=:flat)
    solvation_free_energy(atom_coordinates, n_atoms_per_mol, flat_radii, rs, prefactors, overlap_jump, overlap_slope, delaunay_eps)
end

function solvation_free_energy_in_bounds(
    x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}},
    template_centers::Vector{Matrix{Float64}},
    radii::Vector{Vector{Float64}},
    rs::Float64,
    prefactors::AbstractVector,
    overlap_jump::Float64,
    overlap_slope::Float64,
    bounds::Float64,
    delaunay_eps::Float64
    )
    !in_bounds(x, bounds) && return Inf
    solvation_free_energy(x, template_centers, radii, rs, prefactors, overlap_jump, overlap_slope, delaunay_eps)
end

function _format_measures_output(measures::AbstractVector, prefactors::AbstractVector)
    energy = sum(measures .* [prefactors; 1.0])
    measures_dict = Dict{String, Any}(
        "Vs" => measures[1], "As" => measures[2], "Cs" => measures[3],
        "Xs" => measures[4], "OLs" => measures[5]
    )
    return energy, measures_dict
end

function solvation_free_energy_and_measures(
    x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}},
    template_centers::Vector{Matrix{Float64}},
    radii::Vector{Vector{Float64}},
    rs::Float64,
    prefactors::AbstractVector,
    overlap_jump::Float64,
    overlap_slope::Float64,
    delaunay_eps::Float64
    )
    n_atoms_per_mol, flat_radii = _prepare_args(template_centers, radii)

    flat_realization = Utilities.get_realization(x, template_centers, format=:flat)

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
    template_centers::Vector{Matrix{Float64}},
    radii::Vector{Vector{Float64}},
    rs::Float64,
    prefactors::AbstractVector,
    overlap_jump::Float64,
    overlap_slope::Float64,
    delaunay_eps::Float64
    )
    n_atoms_per_mol, flat_radii = _prepare_args(template_centers, radii)

    flat_realization = Utilities.get_realization(x, template_centers, format=:flat)

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

function solvation_free_energy_and_measures_in_bounds(
    x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}},
    template_centers::Vector{Matrix{Float64}},
    radii::Vector{Vector{Float64}},
    rs::Float64,
    prefactors::AbstractVector,
    overlap_jump::Float64,
    overlap_slope::Float64,
    bounds::Float64,
    delaunay_eps::Float64
    )
    !in_bounds(x, bounds) && return (Inf, Dict{String, Any}())
    return solvation_free_energy_and_measures(
        x, template_centers, radii, rs, prefactors, overlap_jump, overlap_slope, delaunay_eps
    )
end

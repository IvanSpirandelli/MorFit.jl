module AlphaMolWrap
    using CxxWrap
    using AlphaMolWrapper_jll
    
    @wrapmodule(() -> libalphamolwrapper)

    #@wrapmodule(() -> joinpath("../../AlphaMolWrapper/build","libalphamolwrapper"))

    function __init__()
        @initcxx
    end
end

function get_geometric_measures_and_overlap_value(
    atom_coordinates::Vector{Float64}, 
    molecule_size::Int, 
    atom_radii::Vector{Float64}, 
    probe_radius::Float64, 
    overlap_existence_penalty::Float64,
    overlap_penalty_slope::Float64, 
    delaunay_eps::Float64 = 1.0)

    n_mol = length(atom_radii) ÷ molecule_size
    get_geometric_measures_and_overlap_value(
        atom_coordinates,
        fill(molecule_size, n_mol),
        atom_radii,
        probe_radius,
        overlap_existence_penalty,
        overlap_penalty_slope,
        delaunay_eps
    )
end

function get_geometric_measures_and_overlap_value(
    atom_coordinates::Vector{Float64}, 
    molecule_sizes::Vector{Int}, 
    atom_radii::Vector{Float64}, 
    probe_radius::Float64, 
    overlap_existence_penalty::Float64,
    overlap_penalty_slope::Float64, 
    delaunay_eps::Float64 = 1.0)
    outs = [0.0, 0.0, 0.0, 0.0, 0.0]
    AlphaMolWrap.get_geometric_measures_and_overlap_value(
        outs,
        atom_coordinates,
        molecule_sizes,
        atom_radii,
        probe_radius,
        overlap_existence_penalty,
        overlap_penalty_slope,
        delaunay_eps
    )
    outs
end


function get_geometric_measures_and_overlap_value_with_derivatives(
    atom_coordinates::Vector{Float64}, 
    molecule_size::Int, 
    atom_radii::Vector{Float64}, 
    probe_radius::Float64, 
    overlap_existence_penalty::Float64,
    overlap_penalty_slope::Float64, 
    delaunay_eps::Float64 = 1.0)

    n_mol = length(atom_radii) ÷ molecule_size
    get_geometric_measures_and_overlap_value_with_derivatives(
        atom_coordinates,
        fill(molecule_size, n_mol),
        atom_radii,
        probe_radius,
        overlap_existence_penalty,
        overlap_penalty_slope,
        delaunay_eps
    )
end

function get_geometric_measures_and_overlap_value_with_derivatives(
    atom_coordinates::Vector{Float64}, 
    molecule_sizes::Vector{Int},
    atom_radii::Vector{Float64}, 
    probe_radius::Float64, 
    overlap_existence_penalty::Float64,
    overlap_penalty_slope::Float64, 
    delaunay_eps::Float64 = 1.0)
    
    measure_outs = [0.0, 0.0, 0.0, 0.0, 0.0]

    n = size(atom_coordinates)[1]

    dvol_outs = [0.0 for _ in 1:n]
    dsurf_outs = [0.0 for _ in 1:n]
    dmean_outs = [0.0 for _ in 1:n]
    dgauss_outs = [0.0 for _ in 1:n]
    dlol_outs = [0.0 for _ in 1:n]

    AlphaMolWrap.get_geometric_measures_and_overlap_value_with_derivatives(
        measure_outs,
        dvol_outs,
        dsurf_outs, 
        dmean_outs,
        dgauss_outs,
        dlol_outs,
        atom_coordinates,
        molecule_sizes,
        atom_radii,
        probe_radius,
        overlap_existence_penalty,
        overlap_penalty_slope,
        delaunay_eps
    )
    measure_outs, dvol_outs, dsurf_outs, dmean_outs, dgauss_outs, dlol_outs
end

function get_geometric_measures(
    atom_coordinates::Vector{Float64}, 
    atom_radii::Vector{Float64}, 
    probe_radius::Float64,
    delaunay_eps::Float64 = 1.0)
    outs = [0.0, 0.0, 0.0, 0.0]
    AlphaMolWrap.get_geometric_measures(
        outs,
        atom_coordinates,
        atom_radii,
        probe_radius,
        delaunay_eps
    )
    outs
end

function get_geometric_measures_with_derivatives(
    atom_coordinates::Vector{Float64}, 
    atom_radii::Vector{Float64}, 
    probe_radius::Float64, 
    delaunay_eps::Float64 = 1.0)
    
    measure_outs = [0.0, 0.0, 0.0, 0.0, 0.0]

    n = size(atom_coordinates)[1]

    dvol_outs = [0.0 for _ in 1:n]
    dsurf_outs = [0.0 for _ in 1:n]
    dmean_outs = [0.0 for _ in 1:n]
    dgauss_outs = [0.0 for _ in 1:n]

    AlphaMolWrap.get_geometric_measures_with_derivatives(
        measure_outs,
        dvol_outs,
        dsurf_outs, 
        dmean_outs,
        dgauss_outs,
        atom_coordinates,
        atom_radii,
        probe_radius,
        delaunay_eps
    )
    measure_outs, dvol_outs, dsurf_outs, dmean_outs, dgauss_outs
end

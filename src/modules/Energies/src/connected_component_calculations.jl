using Graphs

function get_single_subunit_energy_and_measures(centers::Vector{Matrix{Float64}}, radii::Vector{Vector{Float64}}, rs::Float64, prefactors::AbstractVector, overlap_jump, overlap_slope, delaunay_eps)
    combined = [solvation_free_energy_and_measures([(QuatRotation(exp(Rotations.RotationVecGenerator([0.0, 0.0, 0.0]...))), [0.0, 0.0, 0.0])], [tc], [tr], rs, prefactors, overlap_jump, overlap_slope, delaunay_eps) for (tc,tr) in zip(centers, radii)]
    [c[1] for c in combined], [c[2] for c in combined]
end

function get_bounding_radii(centers::Vector{Matrix{Float64}}, radii::Vector{Vector{Float64}}, rs::Float64)
    max_ds = [maximum([euclidean([0.0, 0.0, 0.0], e) for e in eachcol(tc)]) for tc in centers]
    max_rs = [maximum(r) for r in radii]
    return [d+r+rs for (d,r) in zip(max_ds, max_rs)]
end

function are_bounding_spheres_overlapping(x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}, id_one::Int, id_two::Int, bounding_radii::Vector{Float64})
    t_one = x[id_one][2]
    t_two = x[id_two][2]
    r_one = bounding_radii[id_one]
    r_two = bounding_radii[id_two]
    euclidean(t_one, t_two) < r_one + r_two
end

function solvation_free_energy_and_separate_overlap_with_bounding_container_check(
    x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}},
    centers::Vector{Matrix{Float64}},
    radii::Vector{Vector{Float64}},
    rs::Float64,
    prefactors::AbstractVector,
    overlap_jump::Float64,
    overlap_slope::Float64,
    delaunay_eps::Float64,
    single_subunit_energies::Vector{Float64},
    molecule_boundary_overlap_check::Function
    )
    @assert length(x) == 2 "This function cannot handle more than two molecules"
    if !molecule_boundary_overlap_check(x)
        return single_subunit_energies[1] + single_subunit_energies[2], 0.0
    end
    solvation_free_energy_and_overlap(x, centers, radii, rs, prefactors, overlap_jump, overlap_slope, delaunay_eps)
end

# Stuff for 3 or more subunits
get_sub_state(x, indices) = x[indices]

function construct_overlap_graph(x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}, molecule_boundary_overlap_check)
    n = Int(length(x))
    graph = SimpleGraph(n)
    for i in 1:n
        for j = i+1:n
            if molecule_boundary_overlap_check(x, i, j)
                add_edge!(graph, i, j)
            end
        end
    end
    graph
end

function connected_component_wise_solvation_free_energy_and_measures(
    last_iteration_ccs_energies_and_measures::Dict{Vector{Int64}, Tuple{Float64, Dict{String, Any}}},
    transformed_index::Int,
    x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}},
    centers::Vector{Matrix{Float64}},
    radii::Vector{Vector{Float64}},
    rs::Float64,
    prefactors::AbstractVector,
    overlap_jump::Float64,
    overlap_slope::Float64,
    delaunay_eps::Float64,
    single_subunit_energies::Vector{Float64},
    single_subunit_measures::Vector{Dict{String, Any}},
    molecule_boundary_overlap_check::Function,
    )
    graph = construct_overlap_graph(x, molecule_boundary_overlap_check)
    ccs = connected_components(graph)
    indexed_cc = [e for e in connected_components(graph) if transformed_index in e][1]
    ccs_energies_and_measures = Dict{Vector{Int64}, Tuple{Float64, Dict{String, Any}}}()

    for cc in ccs
        if length(cc) == 1
            mol_idx = cc[1]
            ccs_energies_and_measures[cc] = single_subunit_energies[mol_idx], single_subunit_measures[mol_idx]
        elseif !(cc in keys(last_iteration_ccs_energies_and_measures)) || cc == indexed_cc
            sub_state = get_sub_state(x, cc)
            cc_centers = centers[cc]
            cc_radii = radii[cc]
            ccs_energies_and_measures[cc] = solvation_free_energy_and_measures(
                sub_state, cc_centers, cc_radii,
                rs, prefactors, overlap_jump, overlap_slope, delaunay_eps
            )
        else
            ccs_energies_and_measures[cc] = last_iteration_ccs_energies_and_measures[cc]
        end
    end
    f_sol = 0.0
    measures = Dict{String, Any}(k => 0.0 for k in keys(first(ccs_energies_and_measures)[2][2]))
    for comb_val in values(ccs_energies_and_measures)
        f_sol += comb_val[1]
        for (k,v) in comb_val[2]
            measures[k] += v
        end
    end
    f_sol, measures, ccs_energies_and_measures
end

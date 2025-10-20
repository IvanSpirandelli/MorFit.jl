using Graphs

function get_single_subunit_energy_and_measures(template_centers::Matrix{Float64}, template_radii::Vector{Float64}, rs::Float64, prefactors::AbstractVector, overlap_jump, overlap_slope, delaunay_eps)
    solvation_free_energy_and_measures([(QuatRotation(exp(Rotations.RotationVecGenerator([0.0, 0.0, 0.0]...))), [0.0, 0.0, 0.0])], template_centers, template_radii, rs, prefactors, overlap_jump, overlap_slope, delaunay_eps)
end

function get_single_subunit_energy_and_measures(template_centers::Vector{Matrix{Float64}}, template_radii::Vector{Vector{Float64}}, rs::Float64, prefactors::AbstractVector, overlap_jump, overlap_slope, delaunay_eps)
    combined = [solvation_free_energy_and_measures([(QuatRotation(exp(Rotations.RotationVecGenerator([0.0, 0.0, 0.0]...))), [0.0, 0.0, 0.0])], tc, tr, rs, prefactors, overlap_jump, overlap_slope, delaunay_eps) for (tc,tr) in zip(template_centers, template_radii)]
    [c[1] for c in combined], [c[2] for c in combined]
end

# Molecule of single type
function get_bounding_radii(template_centers::Matrix{Float64}, template_radii::Vector{Float64}, rs::Float64)
    d = maximum([euclidean([0.0, 0.0, 0.0], e) for e in eachcol(template_centers)])
    r = maximum(template_radii)
    return [d+r+rs]
end

# Multiple types of molecules
function get_bounding_radii(template_centers::Vector{Matrix{Float64}}, template_radii::Vector{Vector{Float64}}, rs::Float64)
    max_ds = [maximum([euclidean([0.0, 0.0, 0.0], e) for e in eachcol(tc)]) for tc in template_centers]
    max_rs = [maximum(rs) for rs in template_radii]
    return [d+r+rs for (d,r) in zip(max_ds, max_rs)]
end

function are_bounding_spheres_overlapping(x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}, id_one::Int, id_two::Int, bounding_radii::Vector{Float64})
    t_one = x[id_one][2]
    t_two = x[id_two][2]
    r_one = length(bounding_radii) == 1 ? bounding_radii[1] : bounding_radii[id_one]
    r_two = length(bounding_radii) == 1 ? bounding_radii[1] : bounding_radii[id_two]
    euclidean(t_one, t_two) < r_one + r_two
end

#Use with standard RWM and 2 subunits only
function solvation_free_energy_and_measures_with_bounding_container_check_in_bounds(
    x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}},
    template_centers::Vector{Matrix{Float64}}, 
    radii::Vector{Vector{Float64}},
    rs::Float64, 
    prefactors::AbstractVector, 
    overlap_jump::Float64,
    overlap_slope::Float64,
    bounds::Float64,
    delaunay_eps::Float64,
    single_subunit_energies::Vector{Float64},
    single_subunit_measures::Vector{Dict{String, Any}},
    molecule_boundary_overlap_check::Function
    )
    @assert length(x) == 2 "This function cannot handle more than two molecules"
    if in_bounds(x, bounds)
        solvation_free_energy_and_measures_with_bounding_container_check(x, template_centers, radii, rs, prefactors, overlap_jump, overlap_slope, delaunay_eps, single_subunit_energies, single_subunit_measures, molecule_boundary_overlap_check)
    else
        Inf, Dict{String, Any}()
    end
end

function solvation_free_energy_and_measures_with_bounding_container_check(
    x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}},
    template_centers::Vector{Matrix{Float64}}, 
    radii::Vector{Vector{Float64}},
    rs::Float64, 
    prefactors::AbstractVector, 
    overlap_jump::Float64,
    overlap_slope::Float64,
    delaunay_eps::Float64,
    single_subunit_energies::Vector{Float64},
    single_subunit_measures::Vector{Dict{String, Any}},
    molecule_boundary_overlap_check::Function
    )
    @assert length(x) == 2 "This function cannot handle more than two molecules"
    if !molecule_boundary_overlap_check(x)
        return single_subunit_energies[1] + single_subunit_energies[2], Dict{String, Any}(k => single_subunit_measures[1][k] + single_subunit_measures[2][k] for k in keys(single_subunit_measures[1]))
    end
    solvation_free_energy_and_measures(x, template_centers, radii, rs, prefactors, overlap_jump, overlap_slope, delaunay_eps)
end

function solvation_free_energy_and_separate_overlap_with_bounding_container_check(
    x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}},
    template_centers::Vector{Matrix{Float64}}, 
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
    solvation_free_energy_and_overlap(x, template_centers, radii, rs, prefactors, overlap_jump, overlap_slope, delaunay_eps)
end

function solvation_free_energy_and_measures_with_bounding_container_check_in_bounds(
    x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}},
    template_centers::Matrix{Float64}, 
    radii::Vector{Float64},
    rs::Float64, 
    prefactors::AbstractVector, 
    overlap_jump::Float64,
    overlap_slope::Float64,
    bounds::Float64,
    delaunay_eps::Float64,
    single_subunit_energy::Float64,
    single_subunit_measures::Dict{String, Any},
    molecule_boundary_overlap_check::Function
    )
    @assert length(x) == 2 "This function cannot handle more than two molecules"
    if in_bounds(x, bounds)
        solvation_free_energy_and_measures_with_bounding_container_check(x, template_centers, radii, rs, prefactors, overlap_jump, overlap_slope, delaunay_eps, single_subunit_energy, single_subunit_measures, molecule_boundary_overlap_check)
    else
        Inf, Dict{String, Any}()
    end
end

function solvation_free_energy_and_measures_with_bounding_container_check(
    x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}},
    template_centers::Matrix{Float64}, 
    radii::Vector{Float64},
    rs::Float64, 
    prefactors::AbstractVector, 
    overlap_jump::Float64,
    overlap_slope::Float64,
    delaunay_eps::Float64,
    single_subunit_energy::Float64,
    single_subunit_measures::Dict{String, Any},
    molecule_boundary_overlap_check::Function
    )
    @assert length(x) == 2 "This function cannot handle more than two molecules"
    if !molecule_boundary_overlap_check(x)
        return 2*single_subunit_energy, Dict{String, Any}(k => v.*2 for (k,v) in single_subunit_measures)
    end
    solvation_free_energy_and_measures(x, template_centers, radii, rs, prefactors, overlap_jump, overlap_slope, delaunay_eps)
end

#Stuff for 3 or more subunits
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

function connected_component_wise_solvation_free_energy_and_measures_in_bounds(last_iteration_ccs_energies_and_measures::Dict{Vector{Int64}, Tuple{Float64, Dict{String, Any}}},
    transformed_index::Int,
    x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}},
    template_centers::Matrix{Float64}, 
    template_radii::Vector{Float64},
    rs::Float64, 
    prefactors::AbstractVector, 
    overlap_jump::Float64,
    overlap_slope::Float64,
    bounds::Float64,
    delaunay_eps::Float64,
    single_subunit_energy::Float64,
    single_subunit_measures::Dict{String, Any},
    molecule_boundary_overlap_check::Function,
    )
    if in_bounds(x, bounds)
        connected_component_wise_solvation_free_energy_and_measures(last_iteration_ccs_energies_and_measures, transformed_index, x, template_centers, template_radii, rs, prefactors, overlap_jump, overlap_slope, delaunay_eps, single_subunit_energy, single_subunit_measures, molecule_boundary_overlap_check)
    else
        Inf, Dict{String, Any}(), Dict{Vector{Int64}, Tuple{Float64, Dict{String, Any}}}()
    end
end

function connected_component_wise_solvation_free_energy_and_measures(
    last_iteration_ccs_energies_and_measures::Dict{Vector{Int64}, Tuple{Float64, Dict{String, Any}}},
    transformed_index::Int,
    x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}},
    template_centers::Matrix{Float64}, 
    template_radii::Vector{Float64},
    rs::Float64, 
    prefactors::AbstractVector, 
    overlap_jump::Float64,
    overlap_slope::Float64,
    delaunay_eps::Float64,
    single_subunit_energy::Float64,
    single_subunit_measures::Dict{String, Any},
    molecule_boundary_overlap_check::Function,
    )
    graph = construct_overlap_graph(x, molecule_boundary_overlap_check)
    ccs = connected_components(graph)
    indexed_cc = [e for e in connected_components(graph) if transformed_index in e][1]
    ccs_energies_and_measures = Dict{Vector{Int64}, Tuple{Float64, Dict{String, Any}}}()

    for cc in ccs
        if length(cc) == 1
            ccs_energies_and_measures[cc] = single_subunit_energy, single_subunit_measures
        elseif !(cc in keys(last_iteration_ccs_energies_and_measures)) || cc == indexed_cc
            ccs_energies_and_measures[cc] = solvation_free_energy_and_measures(get_sub_state(x, cc), template_centers, vcat([template_radii for _ in 1:length(cc)]...), rs, prefactors, overlap_jump, overlap_slope, delaunay_eps) 
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

function get_initial_connected_component_energies(
    x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}},
    template_centers::Matrix{Float64}, 
    template_radii::Vector{Float64},
    rs::Float64, 
    prefactors::AbstractVector, 
    overlap_jump::Float64,
    overlap_slope::Float64,
    delaunay_eps::Float64,
    molecule_boundary_overlap_check::Function
    )

    graph = construct_overlap_graph(x, molecule_boundary_overlap_check)
    ccs = connected_components(graph)
    ccs_energies_and_measures = Dict{Vector{Int}, Tuple{Float64, Dict{String, Any}}}()

    for cc in ccs
        ccs_energies_and_measures[cc] = solvation_free_energy_and_measures(get_sub_state(x, cc), template_centers, vcat([template_radii for _ in 1:length(cc)]...), rs, prefactors, overlap_jump, overlap_slope, delaunay_eps)    
    end
    ccs_energies_and_measures
end
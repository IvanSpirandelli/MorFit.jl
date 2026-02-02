function compute_total_alpha_shape_persistence(
    x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}},
    template_centers::Vector{Matrix{Float64}},
    radii::Vector{Vector{Float64}},
    persistence_weights::Vector{Float64},
    exact_delaunay::Bool,
    compute_weighted::Bool
    )
    if compute_weighted
        return total_weighted_alpha_shape_persistence_and_measures(x, template_centers, radii, persistence_weights, exact_delaunay)
    else
        return total_alpha_shape_persistence_and_measures(x, template_centers, persistence_weights, exact_delaunay)
    end
end

function hs_total_alpha_shape_persistence(x::Vector{Float64}, persistence_weights::Vector{Float64}, exact_delaunay = false)
    points = collect(eachcol(reshape(x, (3,length(x)÷3))))
    pdgm = get_alpha_shape_persistence_diagram(points, exact_delaunay)
    
    return _calculate_persistence_energy_and_measures(pdgm, persistence_weights)
end

function total_alpha_shape_persistence_and_measures(
    x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}},
    template_centers::Vector{Matrix{Float64}},
    persistence_weights::Vector{Float64},
    exact_delaunay = false
    )
    points = Utilities.get_realization(x, template_centers, format=:points)
    pdgm = get_alpha_shape_persistence_diagram(points, exact_delaunay)

    return _calculate_persistence_energy_and_measures(pdgm, persistence_weights)
end

_flatten_radii(radii::Vector{Vector{Float64}}) = vcat(radii...)

function total_weighted_alpha_shape_persistence_and_measures(
    x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}},
    template_centers::Vector{Matrix{Float64}},
    radii::Vector{Vector{Float64}},
    persistence_weights::Vector{Float64},
    exact_delaunay = false
    )
    
    points = Utilities.get_realization(x, template_centers, format=:points)
    flat_radii = _flatten_radii(radii) 
    
    pdgm = get_weighted_alpha_shape_persistence_diagram(points, flat_radii, exact_delaunay)
    
    return _calculate_persistence_energy_and_measures(pdgm, persistence_weights)
end

function get_total_persistence(dgm)
    if length(dgm) == 0
        return 0.0
    end
    sum((dgm[:,2] - dgm[:,1]))
end

function _calculate_persistence_energy_and_measures(pdgm, persistence_weights::Vector{Float64})
    p0 = get_total_persistence(pdgm[1])
    p1 = get_total_persistence(pdgm[2])
    p2 = get_total_persistence(pdgm[3])
    
    λ0, λ1, λ2 = persistence_weights
    persistence_energy = λ0 * p0 + λ1 * p1 + λ2 * p2
    
    persistence_dict = Dict{String, Any}("P0s" => p0, "P1s" => p1, "P2s" => p2)
    
    return persistence_energy, persistence_dict
end

function total_weighted_alpha_shape_persistence(
    x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}},
    template_centers::Vector{Matrix{Float64}},
    radii::Vector{Vector{Float64}},
    persistence_weights::Vector{Float64},
    exact_delaunay = false
    )
    
    points = Utilities.get_realization(x, template_centers, format=:points)
    flat_radii = _flatten_radii(radii) 
    
    pdgm = get_weighted_alpha_shape_persistence_diagram(points, flat_radii, exact_delaunay)
    
    return _calculate_persistence_energy(pdgm, persistence_weights)
end


function _calculate_persistence_energy(pdgm, persistence_weights::Vector{Float64})
    p0 = get_total_persistence(pdgm[1])
    p1 = get_total_persistence(pdgm[2])
    p2 = get_total_persistence(pdgm[3])
    
    λ0, λ1, λ2 = persistence_weights
    persistence_energy = λ0 * p0 + λ1 * p1 + λ2 * p2

    return persistence_energy
end

function compute_total_alpha_shape_persistence(
    x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}},
    centers::Vector{Matrix{Float64}},
    radii::Vector{Vector{Float64}},
    persistence_weights::Vector{Float64},
    exact_delaunay::Bool,
    compute_weighted::Bool
    )
    if compute_weighted
        return total_weighted_alpha_shape_persistence_and_measures(x, centers, radii, persistence_weights, exact_delaunay)
    else
        return total_alpha_shape_persistence_and_measures(x, centers, persistence_weights, exact_delaunay)
    end
end

function total_alpha_shape_persistence_and_measures(
    x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}},
    centers::Vector{Matrix{Float64}},
    persistence_weights::Vector{Float64},
    exact_delaunay = false
    )
    points = Utilities.get_realization(x, centers, format=:points)
    pdgm = get_alpha_shape_persistence_diagram(points, exact_delaunay)

    return _calculate_persistence_energy_and_measures(pdgm, persistence_weights)
end

_flatten_radii(radii::Vector{Vector{Float64}}) = vcat(radii...)

function total_weighted_alpha_shape_persistence_and_measures(
    x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}},
    centers::Vector{Matrix{Float64}},
    radii::Vector{Vector{Float64}},
    persistence_weights::Vector{Float64},
    exact_delaunay = false
    )
    points = Utilities.get_realization(x, centers, format=:points)
    flat_radii = _flatten_radii(radii)

    pdgm = get_weighted_alpha_shape_persistence_diagram(points, flat_radii, exact_delaunay)

    return _calculate_persistence_energy_and_measures(pdgm, persistence_weights)
end

function get_total_persistence(dgm)
    isempty(dgm) && return 0.0
    sum(dgm[:,2] .- dgm[:,1])
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

function point_cloud_persistence_energy(
    points::Vector{Vector{Float64}},
    persistence_weights::Vector{Float64},
    exact_delaunay::Bool=false
)
    pdgm = get_alpha_shape_persistence_diagram(points, exact_delaunay)
    return _calculate_persistence_energy_and_measures(pdgm, persistence_weights)
end

function point_cloud_periodic_persistence_energy(
    points::Vector{Vector{Float64}},
    persistence_weights::Vector{Float64},
    box_lower::Vector{Float64},
    box_upper::Vector{Float64},
    exact_delaunay::Bool=false
)
    pdgm = get_periodic_alpha_shape_persistence_diagram(points, box_lower, box_upper, exact_delaunay)
    return _calculate_persistence_energy_and_measures(pdgm, persistence_weights)
end

function point_cloud_persistence_energy_with_edge_penalty(
    points::Vector{Vector{Float64}},
    persistence_weights::Vector{Float64},
    edge_penalty_fn::Function,
    exact_delaunay::Bool=false
)
    pdgm, edges = get_alpha_shape_persistence_diagram_and_edges(points, exact_delaunay)
    E_topo, measures = _calculate_persistence_energy_and_measures(pdgm, persistence_weights)

    E_penalty = 0.0
    for (i, j) in edges
        # diode returns 0-indexed vertex ids; shift to Julia's 1-based indexing
        d = euclidean(points[i+1], points[j+1])
        E_penalty += edge_penalty_fn(d)
    end

    measures["E_penalty"] = E_penalty
    measures["E_T"] = E_topo
    return E_topo + E_penalty, measures
end

# ── Variance penalty ─────────────────────────────────────────────────────────

"""
    get_persistence_variance(dgm) → Float64

Variance of bar lengths (death - birth) for one persistence diagram.
Returns 0.0 for empty or single-bar diagrams.
"""
function get_persistence_variance(dgm)
    size(dgm, 1) < 2 && return 0.0
    bars = dgm[:, 2] .- dgm[:, 1]
    n    = length(bars)
    μ    = sum(bars) / n
    return sum((b - μ)^2 for b in bars) / n
end

"""
    point_cloud_persistence_energy_with_variance(points, λ, γ, exact_delaunay)
        → (E_topo, E_var, measures)

Compute persistence energy and variance penalty from a **single** diagram call.

- `E_topo = λ₀P₀ + λ₁P₁ + λ₂P₂`  (linear; squaring is left to the caller)
- `E_var  = γ₀·Var₀ + γ₁·Var₁ + γ₂·Var₂`

The caller combines them as:
    E_total = (use_squared ? E_topo^2 : E_topo) + E_var

`measures` contains P0s/P1s/P2s plus Var0/Var1/Var2/E_var.
"""
function point_cloud_persistence_energy_with_variance(
    points::Vector{Vector{Float64}},
    persistence_weights::Vector{Float64},
    variance_weights::Vector{Float64},
    exact_delaunay::Bool=false
)
    pdgm   = get_alpha_shape_persistence_diagram(points, exact_delaunay)
    E_topo, measures = _calculate_persistence_energy_and_measures(pdgm, persistence_weights)

    v0 = get_persistence_variance(pdgm[1])
    v1 = get_persistence_variance(pdgm[2])
    v2 = get_persistence_variance(pdgm[3])
    γ0, γ1, γ2 = variance_weights
    E_var = γ0 * v0 + γ1 * v1 + γ2 * v2

    measures["Var0"]  = v0
    measures["Var1"]  = v1
    measures["Var2"]  = v2
    measures["E_var"] = E_var

    return E_topo, E_var, measures
end

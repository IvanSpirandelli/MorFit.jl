"""
    Algorithms

Monte Carlo simulation algorithms for molecular fitting.

Provides Random Walk Metropolis (RWM) implementations for exploring molecular
configuration space using energy functions and perturbation strategies.

# Exported Types
- `RandomWalkMetropolis`: Standard RWM algorithm
- `ConnectedComponentRandomWalkMetropolis`: RWM with connected component tracking

# Main Function
- `simulate!`: Run simulation with given algorithm and initial state
"""
module Algorithms
    using Dates
    using LinearAlgebra
    using Rotations

    """
        add_to_output(measures::Dict, output::Dict{String, Vector})

    Append key-value pairs from `measures` to the corresponding vectors in `output`.

    Only keys that exist in both `measures` and `output` are added.
    """
    function add_to_output(measures::Dict, output::Dict{String, Vector})
        for (k, v) in measures
            if haskey(output, k)
                push!(output[k], v)
            end
        end
    end

    include("random_walk_metropolis.jl")
    include("connected_component_random_walk_metropolis.jl")
end
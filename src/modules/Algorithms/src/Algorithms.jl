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

    Append all key-value pairs from `measures` to the corresponding vectors in `output`.

    Each key in `measures` should exist in `output`, and its value will be pushed
    to the vector at that key.
    """
    function add_to_output(measures::Dict, output::Dict{String, Vector})
        for (k, v) in measures
            push!(output[k], v)
        end
    end

    include("random_walk_metropolis.jl")
    include("connected_component_random_walk_metropolis.jl")
end
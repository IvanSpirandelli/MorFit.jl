"""
    Algorithms

Monte Carlo simulation algorithms for molecular fitting.

Provides Random Walk Metropolis (RWM) implementations for exploring molecular
configuration space using energy functions and perturbation strategies.

# Exported Types
- `RandomWalkMetropolis`: Standard RWM algorithm
- `ConnectedComponentRandomWalkMetropolis`: RWM with connected component tracking
- `SimulationOutput`: Structured output container for simulation trajectories

# Main Function
- `simulate!`: Run simulation with given algorithm and initial state
"""
module Algorithms
    using Dates
    using LinearAlgebra
    using Rotations

    include("simulation_io.jl")
    include("random_walk_metropolis.jl")
    include("connected_component_random_walk_metropolis.jl")
end

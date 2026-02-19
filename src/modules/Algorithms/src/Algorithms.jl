"""
    Algorithms

Monte Carlo simulation algorithms for molecular fitting.

Provides Random Walk Metropolis (RWM) and Simulated Annealing (SA) implementations
for exploring molecular configuration space using energy functions and perturbation strategies.

# Exported Types
- `RandomWalkMetropolis`: Standard RWM algorithm
- `ConnectedComponentRandomWalkMetropolis`: RWM with connected component tracking
- `SimulatedAnnealing`: SA with cooling schedule and adaptive step size
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
    include("simulated_annealing.jl")
    include("calibrate_temperature.jl")
end

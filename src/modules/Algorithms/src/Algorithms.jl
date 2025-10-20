module Algorithms
    using Dates
    using LinearAlgebra
    using Rotations

    function add_to_output(measures::Dict{String, Any}, output::Dict{String, Vector})
        for (k, v) in measures
            push!(output[k], v)
        end
    end

    function add_to_output(measures::Dict{String, Float64}, output::Dict{String, Vector})
        for (k, v) in measures
            push!(output[k], v)
        end
    end    
    
    function add_to_output(measures::Dict{String, Int64}, output::Dict{String, Vector})
        for (k, v) in measures
            push!(output[k], v)
        end
    end

    include("hamiltonian_monte_carlo.jl")
    include("connected_component_random_walk_metropolis.jl")
    include("random_walk_metropolis.jl")
    include("simulated_annealing.jl")
end
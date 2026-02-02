"""
    RandomWalkMetropolis{E, P}

Random Walk Metropolis algorithm for molecular configuration sampling.

# Fields
- `energy::E`: Energy function `(state) -> (energy, measures::Dict)`
- `perturbation::P`: Perturbation function `(state) -> new_state`
- `β::Float64`: Inverse temperature (higher = more selective acceptance)

# Example
```julia
algorithm = RandomWalkMetropolis(
    x -> my_energy(x),
    x -> perturb_single_randomly_chosen(x, σ_r, σ_t),
    1.0
)
output = simulate!(algorithm, initial_state, 60.0, 100000, Dict{String,Vector}())
```
"""
struct RandomWalkMetropolis{E, P}
    energy::E
    perturbation::P
    β::Float64
end

"""
    simulate!(algorithm::RandomWalkMetropolis, x, simulation_time_minutes, target_iterations, output) -> output

Run a fresh RWM simulation starting from state `x`.

# Arguments
- `algorithm`: The RWM algorithm with energy, perturbation, and β
- `x`: Initial state as `Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}`
- `simulation_time_minutes`: Maximum wall-clock time to run
- `target_iterations`: Maximum number of MCMC steps
- `output`: Dictionary to store results (will be populated with trajectory data)

# Output Dictionary Keys
- `"Es"`: Accepted energy values
- `"states"`: Accepted configurations
- `"αs"`: Step number at each acceptance
- `"timestamps"`: Wall-clock time (minutes) at each acceptance
- `"total_step_attempts"`: Total MCMC steps attempted
- Plus any keys from the energy function's measures dict
"""
function simulate!(algorithm::RandomWalkMetropolis, x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}, simulation_time_minutes::Float64, target_iterations::Int, output::Dict{String, Vector})
    start_time = now()
    energy = algorithm.energy
    perturbation = algorithm.perturbation
    β = algorithm.β

    E, measures = energy(x)

    total_step_attempts = 1

    add_to_output(merge!(measures, Dict("Es" => E, "states" => x, "αs" => total_step_attempts, "timestamps" => 0.0)), output)
    
    current_running_time = Dates.value(now() - start_time) / 60000.0
    while current_running_time < simulation_time_minutes && total_step_attempts < target_iterations
        total_step_attempts += 1
        x_cand = perturbation(x)
        E_cand, measures = energy(x_cand)

        if rand() < exp(-β*(E_cand - E))
            # The idea is that at entry i of the array it says at which number of steps m it was accepted. Giving i/m acceptance rate
            E = E_cand
            x = x_cand
            add_to_output(merge!(measures,Dict("Es" => E, "states" => x, "αs" => total_step_attempts, "timestamps" => current_running_time)), output)
        end
        current_running_time = Dates.value(now() - start_time) / 60000.0
    end
    add_to_output(Dict("total_step_attempts" => total_step_attempts), output)
    return output
end

"""
    simulate!(algorithm::RandomWalkMetropolis, input::Dict{String, Any}, output::Dict{String, Vector}) -> (input, output)

Resume an RWM simulation from previous output.

Continues from the last accepted state in `output`, extending the trajectory.
Timestamps are relative to the start of this continuation session.

# Arguments
- `algorithm`: The RWM algorithm with energy, perturbation, and β
- `input`: Dictionary with `"simulation_time_minutes"` and `"target_iterations"`
- `output`: Previous simulation output to continue from

# Returns
Tuple of (input, output) with updated output dictionary.
"""
function simulate!(algorithm::RandomWalkMetropolis, input::Dict{String, Any}, output::Dict{String, Vector})
    if length(output["states"]) == 0
        throw(ArgumentError("output must contain previous simulation states"))
    end

    energy = algorithm.energy
    perturbation = algorithm.perturbation
    β = algorithm.β

    start_time = now()
    x = deepcopy(output["states"][end])
    E = output["Es"][end]

    # Track cumulative time from previous sessions
    previous_time = haskey(output, "timestamps") && !isempty(output["timestamps"]) ? output["timestamps"][end] : 0.0

    total_step_attempts = output["total_step_attempts"][1]

    current_running_time = Dates.value(now() - start_time) / 60000.0
    while current_running_time < input["simulation_time_minutes"] && total_step_attempts < input["target_iterations"]
        total_step_attempts += 1
        x_cand = perturbation(x)
        E_cand, measures = energy(x_cand)

        if rand() < exp(-β*(E_cand - E))
            E = E_cand
            x = x_cand
            add_to_output(merge!(measures, Dict("Es" => E, "states" => x, "αs" => total_step_attempts, "timestamps" => previous_time + current_running_time)), output)
        end
        current_running_time = Dates.value(now() - start_time) / 60000.0
    end
    output["total_step_attempts"] = [total_step_attempts]
    return input, output
end
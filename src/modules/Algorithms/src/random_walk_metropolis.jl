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
output = SimulationOutput()
simulate!(algorithm, initial_state, 60.0, 100000, output)
```
"""
struct RandomWalkMetropolis{E, P}
    energy::E
    perturbation::P
    β::Float64
end

"""
    simulate!(algorithm::RandomWalkMetropolis, x::S, wall_clock_limit_minutes, target_iterations, output::SimulationOutput{S}) -> output

Run a fresh RWM simulation starting from state `x`.

# Arguments
- `algorithm`: The RWM algorithm with energy, perturbation, and β
- `x::S`: Initial state (molecular or point cloud)
- `wall_clock_limit_minutes`: Maximum wall-clock time to run
- `target_iterations`: Maximum number of MCMC steps
- `output`: SimulationOutput to store results

# Output Fields
- `output.E_total`: Accepted energy values
- `output.states`: Accepted configurations
- `output.αs`: Step number at each acceptance
- `output.total_step_attempts`: Total MCMC steps attempted
- `output.measures`: Auto-expanded energy component measures
"""
function simulate!(algorithm::RandomWalkMetropolis, x::S, wall_clock_limit_minutes::Float64, target_iterations::Int, output::SimulationOutput{S}) where S
    start_time = now()
    energy = algorithm.energy
    perturbation = algorithm.perturbation
    β = algorithm.β

    E, measures = energy(x)

    total_step_attempts = 1

    record!(output, E, x, total_step_attempts, measures)

    current_running_time = Dates.value(now() - start_time) / 60000.0
    while current_running_time < wall_clock_limit_minutes && total_step_attempts < target_iterations
        total_step_attempts += 1
        x_cand = perturbation(x)
        E_cand, measures = energy(x_cand)

        if rand() < exp(-β*(E_cand - E))
            E = E_cand
            x = x_cand
            record!(output, E, x, total_step_attempts, measures)
        end
        current_running_time = Dates.value(now() - start_time) / 60000.0
    end
    output.total_step_attempts = total_step_attempts
    return output
end

"""
    simulate!(algorithm::RandomWalkMetropolis, wall_clock_limit_minutes, target_iterations, output::SimulationOutput{S}) -> output

Resume an RWM simulation from previous SimulationOutput.

Continues from the last accepted state, extending the trajectory.

# Arguments
- `algorithm`: The RWM algorithm with energy, perturbation, and β
- `wall_clock_limit_minutes`: Maximum wall-clock time for this continuation
- `target_iterations`: Maximum total MCMC steps (including previous run)
- `output`: Previous SimulationOutput to continue from (must have previous states)
"""
function simulate!(algorithm::RandomWalkMetropolis, wall_clock_limit_minutes::Float64, target_iterations::Int, output::SimulationOutput{S}) where S
    if isempty(output.states)
        throw(ArgumentError("output must contain previous simulation states"))
    end

    energy = algorithm.energy
    perturbation = algorithm.perturbation
    β = algorithm.β

    start_time = now()
    x = deepcopy(output.states[end])
    E = output.E_total[end]

    total_step_attempts = output.total_step_attempts

    current_running_time = Dates.value(now() - start_time) / 60000.0
    while current_running_time < wall_clock_limit_minutes && total_step_attempts < target_iterations
        total_step_attempts += 1
        x_cand = perturbation(x)
        E_cand, measures = energy(x_cand)

        if rand() < exp(-β*(E_cand - E))
            E = E_cand
            x = x_cand
            record!(output, E, x, total_step_attempts, measures)
        end
        current_running_time = Dates.value(now() - start_time) / 60000.0
    end
    output.total_step_attempts = total_step_attempts
    return output
end

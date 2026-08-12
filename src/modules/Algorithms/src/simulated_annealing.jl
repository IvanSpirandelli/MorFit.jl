"""
    SimulatedAnnealing{E, P, C}

Simulated Annealing algorithm for optimization in configuration space.

Uses a cooling schedule to decrease temperature over iterations, with perturbation
step size σ_t scaling as `σ_t_init * √(T/T_init)` (Langevin dynamics justification).

# Fields
- `energy::E`: Energy function `(state) -> (energy, measures::Dict)`
- `perturbation::P`: Perturbation function `(state, σ_t) -> new_state`
- `schedule::C`: Cooling schedule `(step, total_steps) -> T`
- `T_init::Float64`: Initial temperature (for σ_t scaling)
- `σ_t_init::Float64`: Base perturbation magnitude

# Example
```julia
schedule = (step, total) -> quadratic_additive(step, total, 10.0, 0.01)
algorithm = SimulatedAnnealing(
    x -> my_energy(x),
    (x, σ) -> perturb_single_point(x, σ),
    schedule,
    10.0,   # T_init
    0.5     # σ_t_init
)
output = SimulationOutput{PointCloudState}()
simulate!(algorithm, initial_state, 60.0, 10000, output)
```
"""
struct SimulatedAnnealing{E, P, C}
    energy::E
    perturbation::P
    schedule::C
    T_init::Float64
    σ_t_init::Float64
end

function _run_sa_loop!(alg::SimulatedAnnealing, x, E,
                       total_step_attempts::Int, best_state, best_energy::Float64,
                       best_step::Int, all_T::Vector{Float64}, all_E::Vector{Float64},
                       wall_clock_limit_minutes::Float64, target_iterations::Int,
                       output::SimulationOutput{S}) where S
    energy = alg.energy
    perturbation = alg.perturbation
    schedule = alg.schedule
    T_init = alg.T_init
    σ_t_init = alg.σ_t_init

    start_time = now()
    while _elapsed_minutes(start_time) < wall_clock_limit_minutes && total_step_attempts < target_iterations
        total_step_attempts += 1
        T = schedule(total_step_attempts, target_iterations)
        β = 1.0 / T
        σ_t = σ_t_init * sqrt(T / T_init)

        x_cand = perturbation(x, σ_t)
        E_cand, measures = energy(x_cand)

        if rand() < exp(-β * (E_cand - E))
            E = E_cand
            x = x_cand
            measures["T"] = T
            measures["σ_t"] = σ_t
            record!(output, E, x, total_step_attempts, measures)

            if E < best_energy
                best_energy = E
                best_state = deepcopy(x)
                best_step = total_step_attempts
            end
        end

        push!(all_T, T)
        push!(all_E, E)
    end

    output.total_step_attempts = total_step_attempts
    output.metadata["best_state"] = best_state
    output.metadata["best_energy"] = best_energy
    output.metadata["best_step"] = best_step
    output.metadata["all_T"] = all_T
    output.metadata["all_E"] = all_E
    return output
end

"""
    simulate!(alg::SimulatedAnnealing, x::S, wall_clock_limit_minutes, target_iterations, output::SimulationOutput{S}) -> output

Run a fresh Simulated Annealing simulation starting from state `x`.

# Arguments
- `alg`: The SA algorithm with energy, perturbation, schedule, T_init, σ_t_init
- `x::S`: Initial state
- `wall_clock_limit_minutes`: Maximum wall-clock time to run
- `target_iterations`: Maximum number of SA steps
- `output`: SimulationOutput to store results

# Output Fields
- `output.E_total`: Accepted energy values
- `output.states`: Accepted configurations
- `output.αs`: Step number at each acceptance
- `output.total_step_attempts`: Total SA steps attempted
- `output.measures`: Auto-expanded measures including "T" and "σ_t"
- `output.metadata["best_state"]`: Best state found
- `output.metadata["best_energy"]`: Lowest energy found
- `output.metadata["best_step"]`: Step at which best energy was found
- `output.metadata["all_T"]`: Temperature at every step (length = total_step_attempts)
- `output.metadata["all_E"]`: Current energy at every step (after accept/reject)
"""
function simulate!(alg::SimulatedAnnealing, x::S, wall_clock_limit_minutes::Float64,
                   target_iterations::Int, output::SimulationOutput{S}) where S
    E, measures = alg.energy(x)

    T = alg.schedule(1, target_iterations)
    σ_t = alg.σ_t_init * sqrt(T / alg.T_init)
    measures["T"] = T
    measures["σ_t"] = σ_t
    record!(output, E, x, 1, measures)

    all_T = Float64[T]
    all_E = Float64[E]
    sizehint!(all_T, target_iterations)
    sizehint!(all_E, target_iterations)

    _run_sa_loop!(alg, x, E, 1, deepcopy(x), E, 1, all_T, all_E,
                  wall_clock_limit_minutes, target_iterations, output)
end

"""
    simulate!(alg::SimulatedAnnealing, wall_clock_limit_minutes, target_iterations, output::SimulationOutput{S}) -> output

Resume a Simulated Annealing simulation from previous SimulationOutput.

Continues from the last accepted state, extending the step counter from
`output.total_step_attempts`. Per-step tracking vectors (`all_T`, `all_E`)
are restored from metadata and extended.

# Arguments
- `alg`: The SA algorithm
- `wall_clock_limit_minutes`: Maximum wall-clock time for this continuation
- `target_iterations`: Maximum total SA steps (including previous run)
- `output`: Previous SimulationOutput to continue from (must have previous states)
"""
function simulate!(alg::SimulatedAnnealing, wall_clock_limit_minutes::Float64,
                   target_iterations::Int, output::SimulationOutput{S}) where S
    if isempty(output.states)
        throw(ArgumentError("output must contain previous simulation states"))
    end

    x = deepcopy(output.states[end])
    E = output.E_total[end]
    total_step_attempts = output.total_step_attempts

    best_state = get(output.metadata, "best_state", deepcopy(x))
    best_energy = get(output.metadata, "best_energy", E)
    best_step = get(output.metadata, "best_step", total_step_attempts)

    all_T = get(output.metadata, "all_T", Float64[])::Vector{Float64}
    all_E = get(output.metadata, "all_E", Float64[])::Vector{Float64}
    remaining = target_iterations - total_step_attempts
    sizehint!(all_T, length(all_T) + remaining)
    sizehint!(all_E, length(all_E) + remaining)

    _run_sa_loop!(alg, x, E, total_step_attempts, best_state, best_energy, best_step,
                  all_T, all_E, wall_clock_limit_minutes, target_iterations, output)
end

"""
    ConnectedComponentRandomWalkMetropolis{E, P, ICC}

Random Walk Metropolis with connected component tracking for efficient energy updates.

This variant tracks molecular connectivity, allowing incremental energy updates when
only one molecule is perturbed. Useful for systems with many molecules where full
energy recalculation would be expensive.

# Fields
- `energy::E`: Energy function `(components, perturbed_index, state) -> (energy, measures, updated_components)`
- `perturbation::P`: Perturbation function `(state) -> (perturbed_index, new_state)`
- `get_initial_connected_components::ICC`: Function `(state) -> initial_components`
- `β::Float64`: Inverse temperature

# Example
```julia
algorithm = ConnectedComponentRandomWalkMetropolis(
    (cc, i, x) -> cc_energy(cc, i, x),
    x -> get_index_and_perturb_single_randomly_chosen(x, σ_r, σ_t),
    x -> compute_connected_components(x),
    1.0
)
```
"""
struct ConnectedComponentRandomWalkMetropolis{E, P, ICC}
    energy::E
    perturbation::P
    get_initial_connected_components::ICC
    β::Float64
end

function _run_cc_rwm_loop!(alg::ConnectedComponentRandomWalkMetropolis, x, E, ccs,
                           total_step_attempts::Int, wall_clock_limit_minutes::Float64,
                           target_iterations::Int, output::SimulationOutput)
    start_time = now()
    while _elapsed_minutes(start_time) < wall_clock_limit_minutes && total_step_attempts < target_iterations
        total_step_attempts += 1
        i, x_cand = alg.perturbation(x)
        E_cand, measures, updated_ccs = alg.energy(ccs, i, x_cand)

        if rand() < exp(-alg.β * (E_cand - E))
            E = E_cand
            x = x_cand
            ccs = updated_ccs
            record!(output, E, x, total_step_attempts, measures)
        end
    end
    output.total_step_attempts = total_step_attempts
    return output
end

"""
    simulate!(algorithm::ConnectedComponentRandomWalkMetropolis, x, wall_clock_limit_minutes, target_iterations, output::SimulationOutput) -> output

Run a connected component RWM simulation starting from state `x`.

# Arguments
- `algorithm`: The CC-RWM algorithm
- `x`: Initial state as `Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}`
- `wall_clock_limit_minutes`: Maximum wall-clock time to run
- `target_iterations`: Maximum number of MCMC steps
- `output`: SimulationOutput to store results

# Output Fields
Same as `RandomWalkMetropolis`, plus any connected component measures.
"""
function simulate!(algorithm::ConnectedComponentRandomWalkMetropolis, x::MolecularState,
                   wall_clock_limit_minutes::Float64, target_iterations::Int,
                   output::SimulationOutput)
    ccs = algorithm.get_initial_connected_components(x)
    E, measures, _ = algorithm.energy(ccs, 1, x)
    record!(output, E, x, 1, measures)
    _run_cc_rwm_loop!(algorithm, x, E, ccs, 1,
                      wall_clock_limit_minutes, target_iterations, output)
end

simulate!(algorithm::ConnectedComponentRandomWalkMetropolis, x::MolecularState,
          wall_clock_limit_minutes::Float64, output::SimulationOutput) =
    simulate!(algorithm, x, wall_clock_limit_minutes, typemax(Int), output)

"""
    simulate!(algorithm::ConnectedComponentRandomWalkMetropolis, wall_clock_limit_minutes, target_iterations, output::SimulationOutput) -> output

Resume a connected component RWM simulation from previous SimulationOutput.

Continues from the last accepted state, extending the trajectory. Connected
components are recomputed from the resumed state.
"""
function simulate!(algorithm::ConnectedComponentRandomWalkMetropolis,
                   wall_clock_limit_minutes::Float64, target_iterations::Int,
                   output::SimulationOutput)
    if isempty(output.states)
        throw(ArgumentError("output must contain previous simulation states"))
    end
    x = deepcopy(output.states[end])
    E = output.E_total[end]
    ccs = algorithm.get_initial_connected_components(x)
    _run_cc_rwm_loop!(algorithm, x, E, ccs, output.total_step_attempts,
                      wall_clock_limit_minutes, target_iterations, output)
end

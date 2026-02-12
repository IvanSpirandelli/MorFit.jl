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

"""
    simulate!(algorithm::ConnectedComponentRandomWalkMetropolis, x, wall_clock_limit_minutes, output) -> output

Run a connected component RWM simulation starting from state `x`.

# Arguments
- `algorithm`: The CC-RWM algorithm
- `x`: Initial state as `Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}`
- `wall_clock_limit_minutes`: Maximum wall-clock time to run
- `output`: Dictionary to store results

# Output Dictionary Keys
Same as `RandomWalkMetropolis`, plus any connected component measures.
"""
function simulate!(algorithm::ConnectedComponentRandomWalkMetropolis, x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}, wall_clock_limit_minutes::Float64, output::Dict{String, Vector})
    start_time = now()
    energy = algorithm.energy
    perturbation = algorithm.perturbation
    get_initial_connected_components = algorithm.get_initial_connected_components
    β = algorithm.β

    x_cand = deepcopy(x)
    icc = deepcopy(get_initial_connected_components(x))
    E, measures, _ = energy(icc, 1, x)

    total_step_attempts = 1

    add_to_output(merge!(measures, Dict("E_total" => E, "states" => x, "αs" => total_step_attempts)), output)

    current_running_time = Dates.value(now() - start_time) / 60000.0
    while current_running_time < wall_clock_limit_minutes
        total_step_attempts += 1
        i, x_cand = perturbation(x)
        E_cand, measures, ucc = energy(icc, i, x_cand)

        if rand() < exp(-β*(E_cand - E))
            E = E_cand
            x = x_cand
            icc = deepcopy(ucc)
            add_to_output(merge!(measures,Dict("E_total" => E, "states" => x, "αs" => total_step_attempts)), output)
        end
        current_running_time = Dates.value(now() - start_time) / 60000.0
    end
    add_to_output(Dict("total_step_attempts" => total_step_attempts), output)
    return output
end
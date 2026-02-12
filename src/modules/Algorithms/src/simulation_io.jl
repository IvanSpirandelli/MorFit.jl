#=============================================================================
# SimulationOutput: Structured output for MCMC simulations
#
# Replaces manual Dict{String, Vector} initialization with an auto-expanding
# struct that captures all energy function measures without pre-initialization.
=============================================================================#

"""
    SimulationOutput

Structured output container for MCMC simulation trajectories.

Fixed fields store the core simulation data (states, energies, acceptance).
The `measures` dict auto-expands to capture any key returned by the energy
function — no pre-initialization needed.

# Fields
- `states`: Accepted configurations (one per acceptance event)
- `E_total`: Total energy at each accepted state
- `αs`: Step number at each acceptance event
- `total_step_attempts`: Total MCMC steps attempted (scalar)
- `measures`: Auto-expanding dict for energy component measures (E_G, E_O, E_T, Vs, As, etc.)
- `metadata`: Non-trajectory data (rng_state_end, etc.)

# Example
```julia
output = SimulationOutput()
simulate!(algorithm, x_init, 60.0, 10000, output)
output.total_step_attempts  # scalar Int
output.E_total              # Vector{Float64}
output.measures["E_G"]      # Vector{Float64} (auto-created)
```
"""
mutable struct SimulationOutput
    states::Vector{Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}}
    E_total::Vector{Float64}
    αs::Vector{Int}
    total_step_attempts::Int
    measures::Dict{String, Vector{Float64}}
    metadata::Dict{String, Any}
end

"""
    SimulationOutput()

Create an empty SimulationOutput ready to receive simulation data.
"""
function SimulationOutput()
    SimulationOutput(
        Vector{Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}}(),
        Float64[],
        Int[],
        0,
        Dict{String, Vector{Float64}}(),
        Dict{String, Any}(),
    )
end

"""
    SimulationOutput(d::Dict{String, Vector})

Reconstruct a SimulationOutput from a saved Dict (e.g., loaded from JLD2).

Handles legacy key names:
- `"Es"` → `E_total`

Float64 vectors (other than fixed fields) go to `measures`.
Everything else goes to `metadata`.
"""
function SimulationOutput(d::Dict{String, Vector})
    output = SimulationOutput()

    # Fixed fields
    if haskey(d, "states")
        append!(output.states, d["states"])
    end

    e_key = haskey(d, "E_total") ? "E_total" : (haskey(d, "Es") ? "Es" : nothing)
    if !isnothing(e_key)
        append!(output.E_total, d[e_key])
    end

    if haskey(d, "αs")
        append!(output.αs, d["αs"])
    end

    if haskey(d, "total_step_attempts") && !isempty(d["total_step_attempts"])
        output.total_step_attempts = d["total_step_attempts"][1]
    end

    # Remaining keys → measures or metadata
    fixed_keys = Set(["states", "E_total", "Es", "αs", "total_step_attempts"])
    for (k, v) in d
        k in fixed_keys && continue
        if v isa Vector{Float64}
            output.measures[k] = v
        else
            output.metadata[k] = v
        end
    end

    return output
end

"""
    record!(output::SimulationOutput, E, state, step, measures)

Record an accepted MCMC step into the output.

Pushes energy, state, and step number to fixed fields. For each key in
`measures`, auto-creates the vector on first encounter and appends the value.
No pre-initialization required.

# Arguments
- `output`: The simulation output container
- `E::Float64`: Total energy of the accepted state
- `state`: The accepted configuration
- `step::Int`: The MCMC step number at which acceptance occurred
- `measures::Dict`: Energy component measures from the energy function
"""
function record!(output::SimulationOutput, E::Float64,
                 state::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}},
                 step::Int, measures::Dict)
    push!(output.E_total, E)
    push!(output.states, state)
    push!(output.αs, step)
    for (k, v) in measures
        if !haskey(output.measures, k)
            output.measures[k] = Float64[]
        end
        push!(output.measures[k], v)
    end
end

"""
    to_dict(output::SimulationOutput) -> Dict{String, Vector}

Convert SimulationOutput to Dict format for JLD2 serialization.

Produces the same Dict{String, Vector} format that downstream consumers
(distill, visualization, notebooks) expect. Old files remain loadable via
`SimulationOutput(d::Dict{String, Vector})`.
"""
function to_dict(output::SimulationOutput)
    d = Dict{String, Vector}(
        "states" => output.states,
        "E_total" => output.E_total,
        "αs" => output.αs,
        "total_step_attempts" => [output.total_step_attempts],
    )
    for (k, v) in output.measures
        d[k] = v
    end
    for (k, v) in output.metadata
        d[k] = v isa Vector ? v : [v]
    end
    return d
end

"""
    build_input_dict(; molecule_ids, system, sol_params, ol_params, topo_params,
                       num_params, scales, pert_params, perturbation, T_sim,
                       wall_clock_limit_minutes, target_iterations, x_init, seed)

Build the canonical input/config dictionary for saving alongside simulation output.

Centralizes the 5+ identical `saved_config = Dict(...)` blocks scattered across
simulation scripts. The returned Dict{String, Any} is saved to JLD2 as `input`
and used by `continue_simulation` to reconstruct the simulation setup.

All parameters use duck-typing (no Energies module dependency).
"""
function build_input_dict(;
    molecule_ids,
    system,
    sol_params,
    ol_params,
    topo_params,
    num_params,
    scales,
    pert_params,
    perturbation,
    T_sim,
    wall_clock_limit_minutes,
    target_iterations,
    x_init,
    seed=-1,
)
    d = Dict{String, Any}(
        "molecule_ids" => molecule_ids,
        "centers" => system.centers, "radii" => system.radii,
        "rs" => sol_params.rs, "η" => sol_params.η, "bounds" => system.bounds,
        "overlap_jump" => ol_params.jump, "overlap_slope" => ol_params.slope,
        "delaunay_eps" => num_params.delaunay_eps, "exact_delaunay" => num_params.exact_delaunay,
        "persistence_weights" => topo_params.λ,
        "θ_G" => scales.θ_G, "θ_O" => scales.θ_O, "θ_T" => scales.θ_T,
        "σ_r" => pert_params.σ_r, "σ_t" => pert_params.σ_t,
        "perturbation" => perturbation,
        "T_sim" => T_sim,
        "wall_clock_limit_minutes" => wall_clock_limit_minutes,
        "target_iterations" => target_iterations,
        "x_init" => x_init,
    )
    if seed >= 0
        d["seed"] = seed
    end
    return d
end

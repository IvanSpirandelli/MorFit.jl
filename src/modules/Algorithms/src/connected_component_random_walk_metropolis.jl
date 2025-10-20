struct ConnectedComponentRandomWalkMetropolis{E, P, ICC}
    energy::E
    perturbation::P
    get_initial_connected_components::ICC
    β::Float64
end

function simulate!(algorithm::ConnectedComponentRandomWalkMetropolis, x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}, simulation_time_minutes::Float64, output::Dict{String, Vector})
    start_time = now()
    energy = algorithm.energy
    perturbation = algorithm.perturbation
    get_initial_connected_components = algorithm.get_initial_connected_components
    β = algorithm.β

    x_cand = deepcopy(x)
    icc = deepcopy(get_initial_connected_components(x))
    E, measures, _ = energy(icc, 1, x)

    total_step_attempts = 1

    add_to_output(merge!(measures, Dict("Es" => E, "states" => x, "αs" => total_step_attempts, "timestamps" => 0.0)), output)
    
    current_running_time = Dates.value(now() - start_time) / 60000.0
    while current_running_time < simulation_time_minutes
        total_step_attempts += 1
        i, x_cand = perturbation(x)
        E_cand, measures, ucc = energy(icc, i, x_cand)

        if rand() < exp(-β*(E_cand - E))
            E = E_cand
            x = x_cand
            icc = deepcopy(ucc)
            add_to_output(merge!(measures,Dict("Es" => E, "states" => x, "αs" => total_step_attempts, "timestamps" => current_running_time)), output)
        end
        current_running_time = Dates.value(now() - start_time) / 60000.0
    end
    add_to_output(Dict("total_step_attempts" => total_step_attempts), output)
    return output
end
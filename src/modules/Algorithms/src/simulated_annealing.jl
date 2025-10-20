struct SimulatedAnnealing{E, P, TR}
    energy::E
    perturbation::P
    temperature_reduction::TR
end

function simulate!(algorithm::SimulatedAnnealing, x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}, iterations::Int, output::Dict)
    energy = algorithm.energy
    perturbation = algorithm.perturbation
    temperature_reduction = algorithm.temperature_reduction
    
    T = temperature_reduction(0)
    x_cand = deepcopy(x)
    E = energy(x)
    add_to_output(Dict("Es" => E, "states" => x), output)
    accepted_steps = 0
    for i in 1:iterations
        x_cand = perturbation(x)
        E_backup = energy(x_cand)
        if rand() < exp(-(1.0/T)*(E_backup - E))
            E = E_backup
            accepted_steps += 1
            x = deepcopy(x_cand)
            add_to_output(Dict("Es" => E, "states" => x), output)
        end

        T = temperature_reduction(i)
    end

    return x, E, accepted_steps/iterations
end

function simulate!(algorithm::SimulatedAnnealing, x::Vector{Float64}, iterations::Int, output::Dict)
    energy = algorithm.energy
    perturbation = algorithm.perturbation
    temperature_reduction = algorithm.temperature_reduction
    
    T = temperature_reduction(0)
    x_cand = deepcopy(x)
    E = energy(x)
    add_to_output(Dict("Es" => E, "states" => x), output)
    accepted_steps = 0
    for i in 1:iterations
        x_cand = perturbation(x)
        E_backup = energy(x_cand)
        if rand() < exp(-(1.0/T)*(E_backup - E))
            E = E_backup
            accepted_steps += 1
            x = deepcopy(x_cand)
            add_to_output(Dict("Es" => E, "states" => x), output)
        end

        T = temperature_reduction(i)
    end

    return x, E, accepted_steps/iterations
end

function simulate_sample_all!(algorithm::SimulatedAnnealing, x::Vector{Float64}, iterations::Int, output::Dict)
    energy = algorithm.energy
    perturbation = algorithm.perturbation
    temperature_reduction = algorithm.temperature_reduction
    
    T = temperature_reduction(0)
    x_cand = deepcopy(x)
    E = energy(x)
    add_to_output(Dict("Es" => E, "states" => x), output)
    accepted_steps = 0
    for i in 1:iterations
        x_cand = perturbation(x)
        E_backup = energy(x_cand)
        if rand() < exp(-(1.0/T)*(E_backup - E))
            E = E_backup
            accepted_steps += 1
            x = deepcopy(x_cand)
        end
        add_to_output(Dict("Es" => E, "states" => deepcopy(x)), output)

        T = temperature_reduction(i)
    end

    return x, E, accepted_steps/iterations
end

function simulate!(algorithm::SimulatedAnnealing, x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}, simulation_time_minutes::Float64, output::Dict{String, Vector})
    start_time = now()
    energy = algorithm.energy
    perturbation = algorithm.perturbation
    temperature_reduction = algorithm.temperature_reduction
    
    T = temperature_reduction(0.0)

    x_cand = deepcopy(x)

    E, measures = energy(x)

    total_step_attempts = 1
    add_to_output(merge!(measures, Dict("Es" => E, "states" => x, "αs" => total_step_attempts)), output)

    time_passed = Dates.value(now() - start_time) / 60000.0
    while time_passed < simulation_time_minutes
        total_step_attempts += 1
        x_cand = perturbation(x)
        E_cand, measures = energy(x_cand)

        if rand() < exp(-(1.0/T)*(E_cand - E))
            E = E_cand
            x = deepcopy(x_cand)
            add_to_output(merge!(measures,Dict("Es" => E, "states" => x, "αs" => total_step_attempts)), output)
        end
        time_passed = Dates.value(now() - start_time) / 60000.0
        T = temperature_reduction(time_passed)
    end
    add_to_output(Dict("total_step_attempts" => iterations), output)
    return output
end

struct ConnectedComponentSimulatedAnnealing{E, P, TR, ICC}
    energy::E
    perturbation::P
    temperature_reduction::TR
    get_initial_connected_components::ICC
end

function simulate!(algorithm::ConnectedComponentSimulatedAnnealing, x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}, simulation_time_minutes::Float64, output::Dict{String, Vector})
    start_time = now()
    energy = algorithm.energy
    perturbation = algorithm.perturbation
    get_initial_connected_components = algorithm.get_initial_connected_components
    temperature_reduction = algorithm.temperature_reduction
    
    T = temperature_reduction(0.0)

    x_cand = deepcopy(x)
    icc = deepcopy(get_initial_connected_components(x))
    E, measures, _ = energy(icc, 1, x)

    total_step_attempts = 1

    add_to_output(merge!(measures, Dict("Es" => E, "states" => x, "αs" => total_step_attempts)), output)

    time_passed = Dates.value(now() - start_time) / 60000.0
    while time_passed < simulation_time_minutes
        total_step_attempts += 1
        i, x_cand = perturbation(x)
        E_cand, measures, ucc = energy(icc, i, x_cand)

        if rand() < exp(-(1.0/T)*(E_cand - E))
            E = E_cand
            x = deepcopy(x_cand)
            icc = deepcopy(ucc)
            add_to_output(merge!(measures,Dict("Es" => E, "states" => x, "αs" => total_step_attempts)), output)
        end
        time_passed = Dates.value(now() - start_time) / 60000.0
        T = temperature_reduction(time_passed)
    end
    add_to_output(Dict("total_step_attempts" => total_step_attempts), output)
    return output
end

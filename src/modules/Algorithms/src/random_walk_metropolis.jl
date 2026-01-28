struct RandomWalkMetropolis{E, P}
    energy::E
    perturbation::P
    β::Float64
end

function simulate!(algorithm::RandomWalkMetropolis, x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}, simulation_time_minutes::Float64, target_iterations::Int, output::Dict{String, Vector})
    start_time = now()
    energy = algorithm.energy
    perturbation = algorithm.perturbation
    β = algorithm.β

    E, measures = energy(x)

    total_step_attempts = 1

    add_to_output(merge!(measures, Dict("Es" => E, "states" => x, "αs" => total_step_attempts, "timestamps" => 0.0)), output)
    
    current_running_time = Dates.value(now() - start_time) / 60000.0
    while current_running_time < simulation_time_minutes && total_step_attempts <= target_iterations
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

function simulate!(algorithm::RandomWalkMetropolis, x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}, simulation_time_minutes::Float64, output::Dict{String, Vector})
    start_time = now()
    energy = algorithm.energy
    perturbation = algorithm.perturbation
    β = algorithm.β

    E, measures = energy(x)

    total_step_attempts = 1

    add_to_output(merge!(measures, Dict("Es" => E, "states" => x, "αs" => total_step_attempts, "timestamps" => 0.0)), output)
    
    current_running_time = Dates.value(now() - start_time) / 60000.0
    while current_running_time < simulation_time_minutes
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

function simulate!(algorithm::RandomWalkMetropolis, x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}, iterations::Int, output::Dict{String, Vector})
    start_time = now()
    energy = algorithm.energy
    perturbation = algorithm.perturbation
    β = algorithm.β

    E, measures = energy(x)
    add_to_output(merge!(measures, Dict("Es" => E, "states" => x, "αs" => 1, "timestamps" => 0.0)), output)

    for i in 1:iterations
        x_cand = perturbation(x)
        E_cand, measures = energy(x_cand)

        if rand() < exp(-β*(E_cand - E))
            E = E_cand
            x = x_cand
            add_to_output(merge!(measures, Dict("Es" => E, "states" => x, "αs" => i, "timestamps" => Dates.value(now() - start_time) / 60000.0)), output)
        end
    end
    add_to_output(Dict{String, Any}("total_step_attempts" => iterations, "timestamps" => Dates.value(now() - start_time) / 60000.0), output)
    return output
end

function simulate!(algorithm::RandomWalkMetropolis, input::Dict{String, Any}, output::Dict{String, Vector})
    if length(output["states"]) == 0
        @assert "This method expects to have proper simulation output."
    end
    
    println("Resuming Random Walk Metropolis simulation with given input and output dictionaries.")

    energy = algorithm.energy
    perturbation = algorithm.perturbation
    β = algorithm.β


    start_time = now()
    x = deepcopy(output["states"][end])
    E = output["Es"][end]

    total_step_attempts = output["total_step_attempts"][1]
    
    current_running_time = Dates.value(now() - start_time) / 60000.0
    while current_running_time < input["simulation_time_minutes"] && total_step_attempts < input["target_iterations"]
        total_step_attempts += 1
        x_cand = perturbation(x)
        E_cand, measures = energy(x_cand)

        if rand() < exp(-β*(E_cand - E))
            E = E_cand
            x = x_cand
            # The idea is that at entry i of the array it says at which number of steps m it was accepted. Giving i/m acceptance rate
            add_to_output(merge!(measures,Dict("Es" => E, "states" => x, "αs" => total_step_attempts)), output)
        end
        current_running_time = Dates.value(now() - start_time) / 60000.0
    end
    output["total_step_attempts"] = [total_step_attempts]
    return input, output
end
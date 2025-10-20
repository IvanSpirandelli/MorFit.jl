struct RandomWalkMetropolis{E, P}
    energy::E
    perturbation::P
    β::Float64
end

function simulate!(algorithm::RandomWalkMetropolis, x_init::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}, iterations::Int)
    energy = algorithm.energy
    perturbation = algorithm.perturbation
    β = algorithm.β
    x = deepcopy(x_init)

    E = energy(x)
    accepted_steps = 0
    for _ in 1:iterations
        x_cand = perturbation(x)
        E_cand = energy(x_cand)

        if rand() < exp(-β*(E_cand - E))
            E = E_cand
            accepted_steps += 1
            x = x_cand
        end
    end

    return x, E, accepted_steps/iterations
end

function simulate!(algorithm::RandomWalkMetropolis, iterations::Int, output::Dict{String, Vector})
    energy = algorithm.energy
    perturbation = algorithm.perturbation
    β = algorithm.β

    x = deepcopy(output["states"][end])
    E = output["Es"][end]
    
    for _ in 1:iterations
        x_cand = perturbation(x)
        E_cand, measures = energy(x_cand)

        if rand() < exp(-β*(E_cand - E))
            E = E_cand
            x = x_cand
            add_to_output(merge!(measures,Dict("Es" => E, "states" => x)), output)
        end
    end
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

function simulate_no_diagnostics!(algorithm::RandomWalkMetropolis, x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}, simulation_time_minutes::Float64, output::Dict{String, Vector})
    start_time = now()
    energy = algorithm.energy
    perturbation = algorithm.perturbation
    β = algorithm.β

    E = energy(x)

    total_step_attempts = 1

    add_to_output(Dict("Es" => E, "states" => x, "αs" => total_step_attempts, "timestamps" => 0.0), output)
    
    current_running_time = Dates.value(now() - start_time) / 60000.0
    while current_running_time < simulation_time_minutes
        total_step_attempts += 1
        x_cand = perturbation(x)
        E_cand = energy(x_cand)

        if rand() < exp(-β*(E_cand - E))
            # The idea is that at entry i of the array it says at which number of steps m it was accepted. Giving i/m acceptance rate
            E = E_cand
            x = x_cand
            add_to_output(Dict("Es" => E, "states" => x, "αs" => total_step_attempts, "timestamps" => current_running_time), output)
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

function simulate!(algorithm::RandomWalkMetropolis, x::Vector{Float64}, iterations::Int, output::Dict{String, Vector})
    energy = algorithm.energy
    perturbation = algorithm.perturbation
    β = algorithm.β

    E, measures = energy(x)

    accepted_steps = 0
    
    add_to_output(merge!(measures, Dict("Es" => E, "states" => x, "αs" => 1)), output)

    for i in 1:iterations
        x_cand = perturbation(x)
        E_cand, measures = energy(x_cand)

        if rand() < exp(-β*(E_cand - E))
            if !isapprox(E_cand, E)
                accepted_steps += 1
                add_to_output(Dict("αs" => i), output)
            end
            E = E_cand
            x = x_cand
            add_to_output(merge!(measures,Dict("Es" => E, "states" => x)), output)
        end
    end
    add_to_output(Dict("total_step_attempts" => iterations), output)
    return output
end

function simulate_sample_all!(algorithm::RandomWalkMetropolis, x::Vector{Float64}, iterations::Int, output::Dict{String, Vector})
    energy = algorithm.energy
    perturbation = algorithm.perturbation
    β = algorithm.β

    E, measures = energy(x)

    accepted_steps = 0
    
    add_to_output(merge!(measures, Dict("Es" => E, "states" => x, "αs" => 1)), output)

    for i in 1:iterations
        x_cand = perturbation(x)
        E_cand, measures = energy(x_cand)

        if rand() < exp(-β*(E_cand - E))
            if !isapprox(E_cand, E)
                accepted_steps += 1
                add_to_output(Dict("αs" => i), output)
            end
            E = E_cand
            x = x_cand
        end
        add_to_output(merge!(measures,Dict("Es" => E, "states" => deepcopy(x))), output)
    end
    add_to_output(Dict("total_step_attempts" => iterations), output)
    return output
end



function simulate!(input::Dict{String, Any}, output::Dict{String, Vector})
    if length(output["states"]) == 0
        @assert "This method expects to have proper simulation output."
    end

    start_time = now()
    energy = get_energy(input["energy"])
    perturbation = get_perturbation(input["perturbation"])
    β = 1.0 / input["T"]
    x = deepcopy(output["states"])

    total_step_attempts = output["total_step_attempts"]
    
    current_running_time = Dates.value(now() - start_time) / 60000.0
    while current_running_time < simulation_time_minutes
        total_step_attempts += 1
        x_cand = perturbation(x)
        E_cand, measures = energy(x_cand)

        if rand() < exp(-β*(E_cand - E))
            # The idea is that at entry i of the array it says at which number of steps m it was accepted. Giving i/m acceptance rate
            E = E_cand
            x = x_cand
            add_to_output(merge!(measures,Dict("Es" => E, "states" => x, "αs" => total_step_attempts)), output)
        end
        current_running_time = Dates.value(now() - start_time) / 60000.0
    end
    add_to_output(Dict("total_step_attempts" => total_step_attempts), output)
    return input, output
end
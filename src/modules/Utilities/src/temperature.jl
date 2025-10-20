function quadratic_additive(passed_time::Float64, total_time::Float64, T_init::Float64, T_min::Float64)
    T_min + (T_init - T_min)*((total_time - passed_time)/total_time)^2
end

function zig_zag(passed_time::Float64, total_time::Float64, T_init::Float64, T_min::Float64, level::Vector{Float64})
    time_per_level = total_time / length(level)
    current = Int(passed_time ÷ time_per_level) + 1 
    if current == length(level) + 1
        return T_min
    end
    quadratic_additive((current + 1)*time_per_level - passed_time, time_per_level, T_init * level[current], T_min)
end

function get_initial_temperature(input; n_samples=1000, scaling = 0.1)
    energy = get_energy(input)
    n_mol = input["n_mol"]
    bounds = input["bounds"]
    test_Es = [energy(MorphoMol.get_initial_state(n_mol, bounds))[1] for i in 1:n_samples]
    test_Es = [e - minimum(test_Es) for e in test_Es]
    sum(test_Es) / length(test_Es) * scaling
end

# This function takes a sequence of energy evaluations and a a target acceptance rate 
# to compute the temperature needed to achieve the desired targe rate in another simulation.
# This requires other paramters in the simulation to be the same as the original simulation.
function calculate_T0(Es, target_acceptance_rate)
    transitions = []
    for i in 1:length(Es)-1
        if Es[i] > Es[i+1]
            push!(transitions, Es[i])
            push!(transitions, Es[i+1])
        end
    end
    transitions = transitions .- minimum(transitions)
    chi_bar(T) = sum([exp(-transitions[i]/T) for i in 1:2:length(transitions)-1])/sum([exp(-transitions[i]/T) for i in 2:2:length(transitions)])
    χ_0 = target_acceptance_rate
    T_0 = 1.0
    try
        while abs(chi_bar(T_0) - χ_0) > 0.000001
            T_0 = T_0 * (log(chi_bar(T_0)) / log(χ_0 ))
        end
    catch 
        println("No energy decreasing transitions found!")
    end

    if isnan(T_0)
        println("NaN initial energy!")
        return NaN
    else
        return T_0
    end
end
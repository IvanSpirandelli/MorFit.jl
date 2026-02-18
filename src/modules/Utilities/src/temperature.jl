"""
    quadratic_additive(step, iterations, T_init, T_min) -> Float64

Quadratic cooling: `T_min + (T_init - T_min) * ((iterations - step) / iterations)^2`

Ref: http://what-when-how.com/artificial-intelligence/a-comparison-of-cooling-schedules-for-simulated-annealing-artificial-intelligence/
"""
function quadratic_additive(step::Int, iterations::Int, T_init::Float64, T_min::Float64)
    T_min + (T_init - T_min)*((iterations - step)/iterations)^2
end

"""
    zig_zag(zig_cool, step, iterations, T_init, T_min, levels) -> Float64

Zig-zag cooling: divides `iterations` into `length(levels)` equal spans.
Within each span, applies `zig_cool` from `T_init * levels[i]` down to `T_min`.

`zig_cool` is any cooling function with signature `(step, iterations, T_init, T_min)`.

Typical usage: `levels = [1.0, 0.9, 0.8, ..., 0.1]` for 10 decreasing cycles.
"""
function zig_zag(zig_cool, step::Int, iterations::Int, T_init::Float64, T_min::Float64, levels::Vector{Float64})
    iteration_level = [Int(round(i * iterations / length(levels))) for i in 0:length(levels)]
    current = findfirst(x -> x >= step, iteration_level)
    zig_cool(step - iteration_level[current-1], Int(round(iterations / length(levels))), T_init * levels[current-1], T_min)
end

function get_initial_temperature(input; n_samples=1000, scaling = 0.1)
    energy = get_energy(input)
    n_mol = length(input["molecule_ids"])
    bounds = input["bounds"]
    test_Es = [energy(get_initial_state(n_mol, bounds))[1] for i in 1:n_samples]
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

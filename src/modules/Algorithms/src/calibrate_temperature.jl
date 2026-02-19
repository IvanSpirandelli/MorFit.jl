"""
    calibrate_T_init(energy, perturbation, σ_t, x_init, T_candidates;
                     n_steps=1000, target_acceptance=0.8) -> T_init

Calibrate the initial temperature for Simulated Annealing by running short
RWM warmup runs at each candidate temperature with fixed step size `σ_t`.

Returns the interpolated temperature at which the acceptance rate equals
`target_acceptance`. If the target is outside the observed range, returns the
candidate temperature closest to the target rate.

# Arguments
- `energy`: Energy function `(state) -> (energy::Float64, measures::Dict)`
- `perturbation`: Perturbation function `(state, σ_t) -> new_state`
- `σ_t::Float64`: Fixed perturbation magnitude (the desired σ_t_init for SA)
- `x_init::S`: Initial state for the warmup runs (each T starts from a fresh
  copy of this state)
- `T_candidates::Vector{Float64}`: Temperatures to probe (must be sorted ascending)

# Keyword Arguments
- `n_steps::Int=1000`: Number of RWM steps per candidate temperature
- `target_acceptance::Float64=0.8`: Desired acceptance rate

# Returns
- `T_init::Float64`: Interpolated temperature for the target acceptance rate

# Example
```julia
T_init = calibrate_T_init(
    energy, perturbation, 2.5, x_init,
    [0.1, 0.5, 1.0, 2.0, 3.0, 5.0],
    n_steps=1000, target_acceptance=0.8
)
sa = SimulatedAnnealing(energy, perturbation, schedule, T_init, 2.5)
```
"""
function calibrate_T_init(
    energy, perturbation, σ_t::Float64, x_init::S, T_candidates::Vector{Float64};
    n_steps::Int=1000, target_acceptance::Float64=0.8
) where S
    acceptance_rates = Float64[]

    for T in T_candidates
        x = deepcopy(x_init)
        E, _ = energy(x)
        β = 1.0 / T
        n_accepted = 0

        for _ in 1:n_steps
            x_cand = perturbation(x, σ_t)
            E_cand, _ = energy(x_cand)
            if rand() < exp(-β * (E_cand - E))
                E = E_cand
                x = x_cand
                n_accepted += 1
            end
        end

        α = n_accepted / n_steps
        push!(acceptance_rates, α)
    end

    # Log results for diagnostics
    println("Temperature calibration (σ_t=$σ_t, $n_steps steps each):")
    for (T, α) in zip(T_candidates, acceptance_rates)
        marker = abs(α - target_acceptance) < 0.05 ? " <--" : ""
        println("  T=$(rpad(T, 5))  α=$(round(α, digits=4))$marker")
    end

    # Interpolate to find T where α ≈ target_acceptance
    # Acceptance rate is monotonically increasing with T, so scan for bracket
    for i in 1:length(T_candidates)-1
        α_lo, α_hi = acceptance_rates[i], acceptance_rates[i+1]
        T_lo, T_hi = T_candidates[i], T_candidates[i+1]
        if (α_lo <= target_acceptance <= α_hi) || (α_hi <= target_acceptance <= α_lo)
            t = (target_acceptance - α_lo) / (α_hi - α_lo)
            T_init = T_lo + t * (T_hi - T_lo)
            println("  → T_init = $(round(T_init, digits=4)) (interpolated for α=$target_acceptance)")
            return T_init
        end
    end

    # Target outside observed range — pick the closest
    _, best_idx = findmin(abs.(acceptance_rates .- target_acceptance))
    T_init = T_candidates[best_idx]
    println("  → T_init = $(round(T_init, digits=4)) (closest candidate, α=$(round(acceptance_rates[best_idx], digits=4)))")
    return T_init
end

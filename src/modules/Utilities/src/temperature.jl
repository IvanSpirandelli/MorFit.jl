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


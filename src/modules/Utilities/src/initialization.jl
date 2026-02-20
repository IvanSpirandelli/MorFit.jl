function get_initial_state(n_mol::Int, bounds::Float64)
    vcat([[(QuatRotation(exp(Rotations.RotationVecGenerator(randn(3)...))), [rand(Uniform(0.0, bounds)), rand(Uniform(0.0, bounds)), rand(Uniform(0.0, bounds))]) for i in 1:n_mol]...]);
end

function initialize_inhibitor_and_target(bounds::Float64)
   [(QuatRotation(exp(Rotations.RotationVecGenerator(randn(3)...))), [rand(Uniform(0.0, bounds)), rand(Uniform(0.0, bounds)), rand(Uniform(0.0, bounds))]), (QuatRotation(exp(Rotations.RotationVecGenerator(0.0, 0.0, 0.0))), [bounds * 0.5, bounds * 0.5, bounds * 0.5])]
end

function get_initial_point_cloud(n_points::Int, bounds::Float64)
    [rand(Uniform(0.0, bounds), 3) for _ in 1:n_points]
end

function get_initial_point_cloud_sphere(n_points::Int, radius::Float64)
    points = Vector{Float64}[]
    while length(points) < n_points
        p = radius .* (2 .* rand(3) .- 1)
        if sum(p.^2) <= radius^2
            push!(points, p)
        end
    end
    return points
end
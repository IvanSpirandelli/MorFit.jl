"""Random pose: uniformly random position in `[0, bounds]³` with a random orientation."""
random_pose(bounds::Float64) = (random_rotation(1.0), rand(Uniform(0.0, bounds), 3))

function get_initial_state(n_mol::Int, bounds::Float64)
    [random_pose(bounds) for _ in 1:n_mol]
end

function initialize_inhibitor_and_target(bounds::Float64)
    [random_pose(bounds), (one(QuatRotation), fill(bounds / 2, 3))]
end

function get_initial_point_cloud(n_points::Int, bounds::Float64)
    [rand(Uniform(0.0, bounds), 3) for _ in 1:n_points]
end

function get_initial_point_cloud_grid(n_points::Int, bounds::Float64)
    k = round(Int, cbrt(n_points))
    k^3 == n_points || throw(ArgumentError("get_initial_point_cloud_grid requires n_points to be a perfect cube, got $n_points"))
    spacing = bounds / k
    offset = spacing / 2
    points = Vector{Float64}[]
    for i in 0:k-1, j in 0:k-1, l in 0:k-1
        push!(points, [offset + i * spacing, offset + j * spacing, offset + l * spacing])
    end
    return points
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

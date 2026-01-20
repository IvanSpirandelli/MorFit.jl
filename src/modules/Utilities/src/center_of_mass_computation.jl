using LinearAlgebra

const ATOMIC_MASSES = Dict{String, Float64}(
    "H"  => 1.008,
    "C"  => 12.011,
    "N"  => 14.007,
    "O"  => 15.999,
    "P"  => 30.974,
    "S"  => 32.06,
    "SE" => 78.96, 
    "FE" => 55.845,
    "ZN" => 65.38,
    "MG" => 24.305,
    "CA" => 40.078,
    "CL" => 35.45,
)

function get_center_of_mass(coordinates::Vector{Vector{Float64}}, elements::Vector{String})
    num_atoms = length(elements)
    if length(coordinates) != num_atoms
        error("The number of coordinate vectors must match the number of elements.")
    end

    # --- Part A: Calculate Center of Mass (COM) ---
    total_mass = 0.0
    center_of_mass = [0.0, 0.0, 0.0]

    for i in 1:num_atoms
        mass = get(ATOMIC_MASSES, elements[i], 0.0) # Defaults to 0.0 if element not found
        if mass == 0.0
            @warn "Element '$(elements[i])' not found in ATOMIC_MASSES dictionary. It will be ignored."
            continue
        end

        total_mass += mass
        center_of_mass .+= mass .* coordinates[i]
    end

    if total_mass == 0.0
        error("Total mass is zero. Cannot calculate center of mass. Check if elements are in ATOMIC_MASSES.")
    end

    center_of_mass ./= total_mass
    center_of_mass
end

function get_moment_of_inertia(centered_coordinates::Vector{Vector{Float64}}, elements::Vector{String})
    # --- Calculate Moment of Inertia (MOI) Tensor ---
    inertia_tensor = zeros(3, 3)

    for i in 1:num_atoms
        mass = get(ATOMIC_MASSES, elements[i], 0.0)
        if mass == 0.0; continue; end # Skip elements not in the dictionary

        x, y, z = centered_coordinates[i]

        inertia_tensor[1, 1] += mass * (y^2 + z^2)
        inertia_tensor[2, 2] += mass * (x^2 + z^2)
        inertia_tensor[3, 3] += mass * (x^2 + y^2)

        inertia_tensor[1, 2] -= mass * x * y
        inertia_tensor[1, 3] -= mass * x * z
        inertia_tensor[2, 3] -= mass * y * z
    end

    # The tensor is symmetric, so we fill in the other half
    inertia_tensor[2, 1] = inertia_tensor[1, 2]
    inertia_tensor[3, 1] = inertia_tensor[1, 3]
    inertia_tensor[3, 2] = inertia_tensor[2, 3]

    inertia_tensor
end

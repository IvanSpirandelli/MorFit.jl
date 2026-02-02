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

    total_mass = 0.0
    center_of_mass = [0.0, 0.0, 0.0]

    for i in 1:num_atoms
        mass = get(ATOMIC_MASSES, elements[i], 0.0)
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

# Matrix version: coordinates is 3×N matrix
function get_center_of_mass(coordinates::Matrix{Float64}, elements::Vector{String})
    num_atoms = length(elements)
    if size(coordinates, 2) != num_atoms
        error("Number of columns in coordinates must match number of elements.")
    end

    total_mass = 0.0
    center_of_mass = [0.0, 0.0, 0.0]

    for i in 1:num_atoms
        mass = get(ATOMIC_MASSES, elements[i], 0.0)
        if mass == 0.0
            @warn "Element '$(elements[i])' not found in ATOMIC_MASSES dictionary. It will be ignored."
            continue
        end

        total_mass += mass
        center_of_mass .+= mass .* coordinates[:, i]
    end

    if total_mass == 0.0
        error("Total mass is zero. Cannot calculate center of mass. Check if elements are in ATOMIC_MASSES.")
    end

    center_of_mass ./= total_mass
    center_of_mass
end

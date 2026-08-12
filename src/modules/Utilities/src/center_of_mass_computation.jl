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

"""
    get_center_of_mass(coordinates, elements) -> Vector{Float64}

Mass-weighted center of `coordinates` (a 3×N matrix or a vector of `[x, y, z]`
points) with masses looked up by element symbol in `ATOMIC_MASSES`.
Unknown elements are skipped with a warning.
"""
function get_center_of_mass(coordinates::Matrix{Float64}, elements::Vector{String})
    num_atoms = length(elements)
    if size(coordinates, 2) != num_atoms
        throw(ArgumentError("number of columns in coordinates ($(size(coordinates, 2))) must match number of elements ($num_atoms)"))
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
        throw(ArgumentError("total mass is zero — none of the elements were found in ATOMIC_MASSES"))
    end

    center_of_mass ./= total_mass
    center_of_mass
end

get_center_of_mass(coordinates::Vector{Vector{Float64}}, elements::Vector{String}) =
    get_center_of_mass(reduce(hcat, coordinates), elements)

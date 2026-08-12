"""
    get_realization(x, centers; format=:matrix)

Transform template coordinates using rotation-translation pairs to produce realized coordinates.

# Arguments
- `x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}`: Rotation-translation pairs, one per molecule
- `centers::Vector{Matrix{Float64}}`: Template coordinates as 3×N matrices, one per molecule
- `format::Symbol`: Output format (default `:matrix`)

# Formats
- `:matrix` → `Vector{Matrix{Float64}}` - 3×N matrix per molecule (default, most efficient)
- `:flat` → `Vector{Float64}` - all coordinates flattened [x1,y1,z1,x2,y2,z2,...]
- `:points` → `Vector{Vector{Float64}}` - vector of [x,y,z] points
- `:points_per_mol` → `Vector{Vector{Vector{Float64}}}` - points grouped by molecule
- `:point3f` → `Vector{Point3f}` - GeometryBasics Point3f (Float32)
- `:point3f_per_mol` → `Vector{Vector{Point3f}}` - Point3f grouped by molecule

# Example
```julia
x = [(R1, t1), (R2, t2)]  # 2 molecules
templates = [centers1, centers2]  # 3×N matrices

get_realization(x, templates)                    # Vector{Matrix{Float64}}
get_realization(x, templates, format=:flat)      # Vector{Float64}
get_realization(x, templates, format=:points)    # Vector{Vector{Float64}}
```
"""
function get_realization(
    x::Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}},
    centers::Vector{Matrix{Float64}};
    format::Symbol=:matrix
)
    matrices = [R * tc .+ t for ((R, t), tc) in zip(x, centers)]

    if format == :matrix
        return matrices
    elseif format == :flat
        return vec(hcat(matrices...))
    elseif format == :points
        return [Vector{Float64}(col) for mat in matrices for col in eachcol(mat)]
    elseif format == :points_per_mol
        return [[Vector{Float64}(col) for col in eachcol(mat)] for mat in matrices]
    elseif format == :point3f
        return [Point3f(col) for mat in matrices for col in eachcol(mat)]
    elseif format == :point3f_per_mol
        return [[Point3f(col) for col in eachcol(mat)] for mat in matrices]
    else
        throw(ArgumentError("Unknown format: $format. Use :matrix, :flat, :points, :points_per_mol, :point3f, or :point3f_per_mol"))
    end
end

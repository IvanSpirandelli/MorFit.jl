#=============================================================================
# Energy System Types
#
# Core structs for energy calculations, molecular systems, and parameters.
=============================================================================#

#=============================================================================
# Energy Scaling
=============================================================================#

"""
    EnergyScales(; θ_G=1.0, θ_O=1.0, θ_T=1.0)

Scaling factors for the additive energy model:

    E = θ_G · F_sol + θ_O · O + θ_T · Σ(λᵢ · Pᵢ)

# Fields
- `θ_G::Float64`: Geometry scale (solvation free energy F_sol)
- `θ_O::Float64`: Overlap penalty scale
- `θ_T::Float64`: Topology scale (persistence energy)
"""
struct EnergyScales
    θ_G::Float64
    θ_O::Float64
    θ_T::Float64
end

EnergyScales(; θ_G::Float64=1.0, θ_O::Float64=1.0, θ_T::Float64=1.0) = EnergyScales(θ_G, θ_O, θ_T)

#=============================================================================
# Molecular System
=============================================================================#

"""
    MolecularSystem

Molecular system definition containing templates and simulation bounds.

# Fields
- `centers::Vector{Matrix{Float64}}`: 3×N coordinate matrix per molecule (templates)
- `radii::Vector{Vector{Float64}}`: Atomic radii per molecule
- `bounds::Float64`: Simulation box boundary
"""
struct MolecularSystem
    centers::Vector{Matrix{Float64}}
    radii::Vector{Vector{Float64}}
    bounds::Float64
end

"""Number of molecules in the system."""
n_molecules(system::MolecularSystem) = length(system.centers)

#=============================================================================
# Energy Parameters
=============================================================================#

"""
    SolvationParams

Parameters for solvation free energy calculation.
Prefactors are computed internally using the White Bear model.

# Fields
- `rs::Float64`: Probe radius (solvent radius)
- `η::Float64`: Packing fraction
"""
struct SolvationParams
    rs::Float64
    η::Float64
end

"""
    LinearOverlapParams

Parameters for linear overlap penalty calculation.

# Fields
- `jump::Float64`: Penalty jump at contact
- `slope::Float64`: Penalty linear slope with penetration depth
"""
struct LinearOverlapParams
    jump::Float64
    slope::Float64
end

"""
    TopologyParams

Parameters for topology/persistence energy calculation.

# Fields
- `λ::Vector{Float64}`: Persistence weights [λ₀, λ₁, λ₂] for H₀, H₁, H₂
"""
struct TopologyParams
    λ::Vector{Float64}
end

#=============================================================================
# Numerical Parameters
=============================================================================#

"""
    NumericalParams

Numerical/algorithmic parameters for Delaunay triangulation.
These control precision in geometric computations and will be consolidated
when solvation and topology calculations are unified.

# Fields
- `delaunay_eps::Float64`: Tolerance for solvation Delaunay triangulation
- `exact_delaunay::Bool`: Use exact arithmetic for topology Delaunay
"""
struct NumericalParams
    delaunay_eps::Float64
    exact_delaunay::Bool
end

NumericalParams(; delaunay_eps::Float64=1.0, exact_delaunay::Bool=false) = NumericalParams(delaunay_eps, exact_delaunay)

#=============================================================================
# Precomputed Values
=============================================================================#

"""
    Precomputed

Precomputed values for efficient connected component energy calculations.

# Fields
- `bounding_radii::Vector{Float64}`: Bounding sphere radius per molecule
- `single_energies::Vector{Float64}`: Energy of each isolated molecule
- `single_measures::Vector{Dict{String, Any}}`: Measures of each isolated molecule
"""
struct Precomputed
    bounding_radii::Vector{Float64}
    single_energies::Vector{Float64}
    single_measures::Vector{Dict{String, Any}}
end

#=============================================================================
# Perturbation Parameters
=============================================================================#

"""
    PerturbationParams

Parameters for MCMC perturbation moves.

# Fields
- `σ_r::Float64`: Rotation step size (radians)
- `σ_t::Float64`: Translation step size (Ångström)
"""
struct PerturbationParams
    σ_r::Float64
    σ_t::Float64
end

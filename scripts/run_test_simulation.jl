#=
Run a minimal test simulation for 50 iterations.
Saves to Simulation_Results/test.jld2

Usage:
    source ../python_environments/python_tda_313/bin/activate
    julia --project=. scripts/run_test_simulation.jl
=#

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

ENV["PYCALL_JL_RUNTIME_PYTHON"] = Sys.which("python3")

using MorFit
using JLD2
using Rotations

# Configuration
molecule_ids = ["1a30:protein", "1a30:ligand"]
bounds = 100.0
rs = 1.4
η = 0.0073
μ = 0.5
overlap_jump = 1000.0
overlap_slope = 11.0
delaunay_eps = 0.001
temperature = 13.0
σ_r = 0.04
σ_t = 0.9
persistence_weights = [0.0, 0.0, 0.0]
exact_delaunay = false
target_iterations = 50

# Load templates
template_centers = [MorFit.MOLECULE_DATA[mol_id].centers for mol_id in molecule_ids]
template_radii = [MorFit.MOLECULE_DATA[mol_id].radii for mol_id in molecule_ids]

# Compute prefactors and bounding radii
prefactors = MorFit.Energies.get_prefactors(rs, η)
bounding_radii = MorFit.Energies.get_bounding_radii(template_centers, template_radii, rs)

# Bounding overlap check
bol = x -> MorFit.Energies.are_bounding_spheres_overlapping(x, 1, 2, bounding_radii)

# Energy function (cc_fsol - solvation only, no Python required)
ssu_energies, ssu_measures = MorFit.Energies.get_single_subunit_energy_and_measures(
    template_centers, template_radii, rs, prefactors, overlap_jump, overlap_slope, delaunay_eps
)
energy = x -> MorFit.Energies.solvation_free_energy_and_measures_with_bounding_container_check_in_bounds(
    x, template_centers, template_radii, rs, prefactors, overlap_jump, overlap_slope, bounds,
    delaunay_eps, ssu_energies, ssu_measures, bol
)

# Perturbation function
perturbation = x -> MorFit.Utilities.perturb_single_specified(x, σ_r, σ_t; specified_index=2)

# Initial state - place ligand away from protein to avoid overlap
x_init = [
    (QuatRotation(1.0, 0.0, 0.0, 0.0), [0.0, 0.0, 0.0]),   # protein at origin
    (QuatRotation(1.0, 0.0, 0.0, 0.0), [30.0, 0.0, 0.0])   # ligand displaced
]

# Output storage (keys match what energy function returns: Vs, As, Cs, Xs, OLs)
output = Dict{String, Vector}(
    "states" => Vector{Vector{Tuple{QuatRotation{Float64}, Vector{Float64}}}}(),
    "Es" => Float64[],
    "Vs" => Float64[],
    "As" => Float64[],
    "Cs" => Float64[],
    "Xs" => Float64[],
    "OLs" => Float64[],
    "timestamps" => Float16[],
    "αs" => Int[],
    "total_step_attempts" => Int[],
)

# Run simulation
β = 1.0 / temperature
rwm = MorFit.Algorithms.RandomWalkMetropolis(energy, perturbation, β)

println("Running test simulation for $target_iterations iterations...")
MorFit.Algorithms.simulate!(rwm, deepcopy(x_init), 60.0, target_iterations, output)
println("Done! Sampled $(length(output["Es"])) states.")

# Save
input = Dict(
    "molecule_ids" => molecule_ids,
    "centers" => template_centers,
    "radii" => template_radii,
    "rs" => rs,
    "η" => η,
    "μ" => μ,
    "bounds" => bounds,
    "overlap_jump" => overlap_jump,
    "overlap_slope" => overlap_slope,
    "delaunay_eps" => delaunay_eps,
    "T" => temperature,
    "σ_r" => σ_r,
    "σ_t" => σ_t,
    "prefactors" => prefactors,
    "persistence_weights" => persistence_weights,
)

output_dir = joinpath(@__DIR__, "../../Simulation_Results")
mkpath(output_dir)
output_file = joinpath(output_dir, "test.jld2")

@save output_file input output
println("Saved to $output_file")

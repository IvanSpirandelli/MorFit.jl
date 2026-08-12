#=
Run a minimal test simulation for 50 iterations.
Saves to simulation-results/test.jld2

Usage:
    source ../python-environments/python-tda-313/bin/activate
    julia --project=. scripts/run_test_simulation.jl
=#

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

ENV["PYCALL_JL_RUNTIME_PYTHON"] = Sys.which("python3")

using MorFit
using MorFit.Energies
using JLD2
using Rotations

# Configuration
molecule_ids = ["1a30:protein", "1a30:ligand"]
target_iterations = 50

system = MolecularSystem(
    [MorFit.MOLECULE_DATA[mol_id].centers for mol_id in molecule_ids],
    [MorFit.MOLECULE_DATA[mol_id].radii for mol_id in molecule_ids],
    100.0  # bounds
)
sol_params = SolvationParams(1.4, 0.0073)
ol_params = LinearOverlapParams(1000.0, 11.0)
topo_params = TopologyParams([0.0, 0.0, 0.0])
num_params = NumericalParams(delaunay_eps=0.001, exact_delaunay=false)
scales = EnergyScales(θ_G=1.0, θ_O=1.0, θ_T=0.0)
pert_params = PerturbationParams(0.04, 0.9)
temperature = 13.0

# Precompute single-subunit values and bounding overlap check
prefactors = MorFit.Energies.get_wb_prefactors(sol_params.rs, sol_params.η)
bounding_radii = MorFit.Energies.get_bounding_radii(system.centers, system.radii, sol_params.rs)
single_energies, single_measures = MorFit.Energies.get_single_subunit_energy_and_measures(
    system.centers, system.radii, sol_params.rs, prefactors,
    ol_params.jump, ol_params.slope, num_params.delaunay_eps
)
precomputed = Precomputed(bounding_radii, single_energies, single_measures)
bol = x -> MorFit.Energies.are_bounding_spheres_overlapping(x, 1, 2, bounding_radii)

energy = x -> MorFit.Energies.calculate_combined_energy(
    x, system, sol_params, ol_params, topo_params, num_params, scales, precomputed, bol
)

perturbation = x -> MorFit.Utilities.perturb_single_specified(x, pert_params.σ_r, pert_params.σ_t; specified_index=2)

# Initial state - place ligand away from protein to avoid overlap
x_init = [
    (one(QuatRotation), [0.0, 0.0, 0.0]),   # protein at origin
    (one(QuatRotation), [30.0, 0.0, 0.0])   # ligand displaced
]

output = MorFit.Algorithms.SimulationOutput()

β = 1.0 / temperature
rwm = MorFit.Algorithms.RandomWalkMetropolis(energy, perturbation, β)

println("Running test simulation for $target_iterations iterations...")
MorFit.Algorithms.simulate!(rwm, deepcopy(x_init), 60.0, target_iterations, output)
println("Done! Sampled $(length(output.E_total)) states.")

# Save using build_input_dict
saved_config = MorFit.Algorithms.build_input_dict(
    molecule_ids=molecule_ids, system=system, sol_params=sol_params,
    ol_params=ol_params, topo_params=topo_params, num_params=num_params,
    scales=scales, pert_params=pert_params, perturbation="single_specified",
    T_sim=temperature, wall_clock_limit_minutes=60.0,
    target_iterations=target_iterations, x_init=x_init,
)

output_dir = joinpath(@__DIR__, "../../simulation-results")
mkpath(output_dir)
output_file = joinpath(output_dir, "test.jld2")

output_dict = MorFit.Algorithms.to_dict(output)
@save output_file input=saved_config output=output_dict
println("Saved to $output_file")

# MorFit.jl

A Julia package for molecular fitting: sampling rigid-body configurations of
molecules (e.g. protein–ligand pairs) with Markov chain Monte Carlo, driven by
morphometric and topological energy functions.

## Overview

A configuration ("state") assigns each molecule a rigid-body pose — a rotation
and a translation applied to a template of atom positions. MorFit explores this
configuration space with Random Walk Metropolis or Simulated Annealing, scoring
states with an additive energy model:

```
E = θ_G · F_sol + θ_O · O + θ_T · Σᵢ λᵢ · Pᵢ
```

| Term | Meaning |
|------|---------|
| `F_sol` | Solvation free energy via the morphometric approach: a linear combination of the ball union's volume, surface area, and integrated mean and Gaussian curvature, with thermodynamic prefactors from the White Bear mark II functional |
| `O` | Overlap penalty (jump + linear slope in penetration depth) |
| `Pᵢ` | Total persistence of the weighted alpha-shape filtration in homological dimension *i*, weighted by `λᵢ` |
| `θ_G, θ_O, θ_T` | Scaling factors for the three components |

Geometric measures are computed with AlphaMol (P. Koehl et al., via
`AlphaMolWrapper_jll`); persistent homology uses the Python libraries
[diode](https://github.com/mrzv/diode) and [oineus](https://github.com/anigmetov/oineus)
through PyCall.

## Package structure

| Module | Contents |
|--------|----------|
| `MorFit.Algorithms` | `RandomWalkMetropolis`, `SimulatedAnnealing` (with cooling schedules and temperature calibration), `ConnectedComponentRandomWalkMetropolis`, `SimulationOutput` container, simulation I/O |
| `MorFit.Energies` | Parameter structs, combined energy, solvation free energy, overlap, persistence computations |
| `MorFit.Utilities` | State realization, perturbations, initialization, RMSD to experimental reference structures |

The package ships template data for benchmark systems: `MorFit.MOLECULE_DATA`
(origin-centered building blocks, keyed `"pdb_id:component"`) and
`MorFit.EXPERIMENTAL_ASSEMBLIES` (experimentally determined reference
configurations for RMSD evaluation).

## Installation

Requires Julia ≥ 1.7 and a Python environment providing `numpy`, `diode`, and
`oineus` for the topological energy terms.

```julia
using Pkg
Pkg.develop(path="path/to/MorFit.jl")
```

Point PyCall at the Python environment before launching Julia:

```bash
export PYCALL_JL_RUNTIME_PYTHON=/path/to/venv/bin/python3
```

## Quickstart

```julia
using MorFit
using MorFit.Energies
using Rotations

# System: a protein–ligand pair from the shipped template data
molecule_ids = ["1a30:protein", "1a30:ligand"]
system = MolecularSystem(
    [MorFit.MOLECULE_DATA[id].centers for id in molecule_ids],
    [MorFit.MOLECULE_DATA[id].radii for id in molecule_ids],
    100.0,  # simulation box bounds
)

# Energy parameters
sol_params = SolvationParams(1.4, 0.3665)        # probe radius, packing fraction
ol_params  = LinearOverlapParams(1000.0, 11.0)   # overlap jump, slope
topo_params = TopologyParams([1.0, 0.0, 0.0])    # persistence weights λ₀, λ₁, λ₂
num_params = NumericalParams()
scales = EnergyScales(θ_G=1.0, θ_O=1.0, θ_T=0.0)

# Precompute single-molecule energies and the bounding-sphere overlap check
prefactors = Energies.get_wb_prefactors(sol_params.rs, sol_params.η)
bounding_radii = Energies.get_bounding_radii(system.centers, system.radii, sol_params.rs)
single_energies, single_measures = Energies.get_single_subunit_energy_and_measures(
    system.centers, system.radii, sol_params.rs, prefactors,
    ol_params.jump, ol_params.slope, num_params.delaunay_eps)
precomputed = Precomputed(bounding_radii, single_energies, single_measures)
bol_check = x -> Energies.are_bounding_spheres_overlapping(x, 1, 2, bounding_radii)

energy = x -> Energies.calculate_combined_energy(
    x, system, sol_params, ol_params, topo_params, num_params, scales, precomputed, bol_check)
perturbation = x -> MorFit.Utilities.perturb_single_specified(x, 0.04, 0.9; specified_index=2)

# Sample: 10⁴ steps of Random Walk Metropolis at temperature T = 13
x_init = [(one(QuatRotation), [0.0, 0.0, 0.0]), (one(QuatRotation), [30.0, 0.0, 0.0])]
output = MorFit.Algorithms.SimulationOutput()
rwm = MorFit.Algorithms.RandomWalkMetropolis(energy, perturbation, 1.0 / 13.0)
MorFit.Algorithms.simulate!(rwm, x_init, 60.0, 10_000, output)

# Inspect the trajectory
output.E_total                # energies of accepted states
output.states                 # accepted configurations
MorFit.get_min_rmsd(          # RMSD of the minimum-energy state vs. experiment
    Dict("molecule_ids" => molecule_ids, "centers" => system.centers),
    output)
```

A complete runnable version is in [`scripts/run_test_simulation.jl`](scripts/run_test_simulation.jl).

## Tests

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

The energy and integration tests require the Python environment (see
Installation). `scripts/test_persistence_against_gudhi.jl` additionally
validates the persistence computations against [GUDHI](https://gudhi.inria.fr/).

## Related repositories

- **mor-fit-hpc** — runs MorFit simulations on SLURM clusters, including a
  Bayesian-optimization pipeline for parameter tuning
- **mor-fit-analysis** — notebooks for analyzing and visualizing simulation results

## References

- H. Hansen-Goos and R. Roth, *Density functional theory for hard-sphere
  mixtures: the White Bear version mark II* (2006) — solvation prefactors
- P. Koehl et al., *AlphaMol* — geometric measures of unions of balls

## License

MIT — see [LICENSE](LICENSE).

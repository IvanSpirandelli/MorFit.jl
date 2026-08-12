using Rotations
using JLD2

function run_configuration_distance_tests()
    @testset verbose = true "Configuration Distances" begin
        @testset verbose = true "Fixed Target-Inhibitor RMSD" begin
            test_fixed_target_inhibitor_rmsd()
        end
        @testset verbose = true "Simulation RMSD Calculation" begin
            test_simulation_rmsd_calculation()
        end
    end
end

"""
Create an identity state (no rotation, no translation) for n molecules.
"""
function identity_state(n::Int)
    [(one(QuatRotation), [0.0, 0.0, 0.0]) for _ in 1:n]
end

"""
Create a state where molecule `mol_idx` is translated by `translation`.
"""
function translated_state(n::Int, mol_idx::Int, translation::Vector{Float64})
    state = identity_state(n)
    R, t = state[mol_idx]
    state[mol_idx] = (R, t .+ translation)
    return state
end

function test_fixed_target_inhibitor_rmsd()
    # Test: When ligand (molecule 2) is translated by distance d,
    # the RMSD should equal d (since target alignment is identity)

    # Use real templates from MOLECULE_DATA
    molecule_ids = ["1a30:protein", "1a30:ligand"]

    @test haskey(MorFit.MOLECULE_DATA, molecule_ids[1])
    @test haskey(MorFit.MOLECULE_DATA, molecule_ids[2])

    # Get centered templates
    templates = [MorFit.MOLECULE_DATA[mol_id].centers for mol_id in molecule_ids]

    # Reference state: identity (centered templates at origin)
    state_ref = identity_state(2)

    # Test 1: Identical states should have RMSD = 0
    state_sim = identity_state(2)
    rmsd = MorFit.Utilities.get_rmsd_for_fixed_target_inhibitor_pair(
        templates, templates, state_sim, state_ref
    )
    @test rmsd ≈ 0.0 atol=1e-10

    # Test 2: Translate ligand by [1.0, 0.0, 0.0] → RMSD should be 1.0
    state_sim = translated_state(2, 2, [1.0, 0.0, 0.0])
    rmsd = MorFit.Utilities.get_rmsd_for_fixed_target_inhibitor_pair(
        templates, templates, state_sim, state_ref
    )
    @test rmsd ≈ 1.0 atol=1e-10

    # Test 3: Translate ligand by [0.0, 2.0, 0.0] → RMSD should be 2.0
    state_sim = translated_state(2, 2, [0.0, 2.0, 0.0])
    rmsd = MorFit.Utilities.get_rmsd_for_fixed_target_inhibitor_pair(
        templates, templates, state_sim, state_ref
    )
    @test rmsd ≈ 2.0 atol=1e-10

    # Test 4: Translate ligand by [1.0, 1.0, 1.0] → RMSD should be sqrt(3)
    state_sim = translated_state(2, 2, [1.0, 1.0, 1.0])
    rmsd = MorFit.Utilities.get_rmsd_for_fixed_target_inhibitor_pair(
        templates, templates, state_sim, state_ref
    )
    @test rmsd ≈ sqrt(3.0) atol=1e-10

    # Test 5: Translate both molecules together - alignment on target compensates, ligand RMSD = 0
    # (When both molecules are shifted by the same amount, relative position is preserved)
    state_sim = [
        (QuatRotation(1.0, 0.0, 0.0, 0.0), [5.0, 0.0, 0.0]),  # target
        (QuatRotation(1.0, 0.0, 0.0, 0.0), [5.0, 0.0, 0.0])   # ligand
    ]
    rmsd = MorFit.Utilities.get_rmsd_for_fixed_target_inhibitor_pair(
        templates, templates, state_sim, state_ref
    )
    @test rmsd ≈ 0.0 atol=1e-10
end

function test_simulation_rmsd_calculation()
    # Load test simulation file
    test_file = joinpath(@__DIR__, "test_simulation.jld2")
    @test isfile(test_file)

    @load test_file input output

    # Extract simulation data
    molecule_ids = Vector{String}(input["molecule_ids"])
    templates = haskey(input, "centers") ? input["centers"] : input["template_centers"]
    n_states = length(output["states"])

    @test molecule_ids == ["1a30:protein", "1a30:ligand"]
    @test n_states == 31

    # Calculate RMSD for all states
    rmsds = [MorFit.get_rmsd_to_ground_truth(molecule_ids, templates, state)
             for state in output["states"]]

    # All RMSDs should be finite and positive
    @test all(isfinite, rmsds)
    @test all(r -> r > 0, rmsds)

    # Test get_min_rmsd returns the minimum RMSD among all states
    min_rmsd = MorFit.get_min_rmsd(input, output)
    min_energy_idx = argmin(output["E_total"])
    @test min_rmsd ≈ rmsds[min_energy_idx]

    # Test get_min_rmsd_cutoff with a cutoff
    cutoff = output["αs"][15]  # Use iteration of 15th state as cutoff
    rmsd_cutoff = MorFit.get_min_rmsd_cutoff(input, output, cutoff)
    @test isfinite(rmsd_cutoff)
    @test rmsd_cutoff >= minimum(rmsds[1:15])
end

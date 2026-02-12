"""
Integration test for the full simulation pipeline.

Tests:
1. Parameter struct creation
2. Energy function with additive model
3. RWM simulation execution
4. Save/load with JLD2
5. Simulation continuation
"""

function run_integration_tests()
    # Test output directory
    test_dir = joinpath(@__DIR__, "test_output")
    mkpath(test_dir)

    # Load test molecule data
    molecule_ids = ["4ty7:protein", "4ty7:ligand"]
    centers = [MorFit.MOLECULE_DATA[id].centers for id in molecule_ids]
    radii = [MorFit.MOLECULE_DATA[id].radii for id in molecule_ids]

    @testset "Parameter Structs" begin
        # Test struct creation with defaults
        scales = EnergyScales()
        @test scales.θ_G == 1.0
        @test scales.θ_O == 1.0
        @test scales.θ_T == 1.0

        # Test struct creation with custom values
        scales2 = EnergyScales(θ_G=1.0, θ_O=1.1, θ_T=0.5)
        @test scales2.θ_O == 1.1

        # Test other structs
        sol = SolvationParams(1.4, 0.3665)
        @test sol.rs == 1.4
        @test sol.η == 0.3665

        ol = LinearOverlapParams(0.0, 100.0)
        @test ol.jump == 0.0

        topo = TopologyParams([1.0, 0.0, 0.0])
        @test length(topo.λ) == 3

        num = NumericalParams()
        @test num.delaunay_eps == 1.0
        @test num.exact_delaunay == false

        system = MolecularSystem(centers, radii, 250.0)
        @test n_molecules(system) == 2
    end

    @testset "Precomputation" begin
        system = MolecularSystem(centers, radii, 250.0)
        sol = SolvationParams(1.4, 0.3665)
        ol = LinearOverlapParams(0.0, 100.0)
        num = NumericalParams(100.0, false)

        prefactors = Energies.get_wb_prefactors(sol.rs, sol.η)
        @test length(prefactors) == 4

        bounding_radii = Energies.get_bounding_radii(centers, radii, sol.rs)
        @test length(bounding_radii) == 2

        single_energies, single_measures = Energies.get_single_subunit_energy_and_measures(
            centers, radii, sol.rs, prefactors, ol.jump, ol.slope, num.delaunay_eps
        )
        @test length(single_energies) == 2
        @test length(single_measures) == 2
    end

    @testset "Energy Calculation" begin
        system = MolecularSystem(centers, radii, 250.0)
        sol = SolvationParams(1.4, 0.3665)
        ol = LinearOverlapParams(0.0, 100.0)
        topo = TopologyParams([1.0, 0.0, 0.0])
        num = NumericalParams(100.0, false)
        scales = EnergyScales(θ_G=1.0, θ_O=1.1, θ_T=0.5)

        prefactors = Energies.get_wb_prefactors(sol.rs, sol.η)
        bounding_radii = Energies.get_bounding_radii(centers, radii, sol.rs)
        single_energies, single_measures = Energies.get_single_subunit_energy_and_measures(
            centers, radii, sol.rs, prefactors, ol.jump, ol.slope, num.delaunay_eps
        )
        precomputed = Precomputed(bounding_radii, single_energies, single_measures)

        bol_check = x -> Energies.are_bounding_spheres_overlapping(x, 1, 2, bounding_radii)

        # Create initial state (well-separated)
        x = [
            (QuatRotation(1.0, 0.0, 0.0, 0.0), [50.0, 50.0, 50.0]),
            (QuatRotation(1.0, 0.0, 0.0, 0.0), [150.0, 150.0, 150.0])
        ]

        energy, measures = Energies.calculate_combined_energy(
            x, system, sol, ol, topo, num, scales, precomputed, bol_check
        )

        @test isfinite(energy)
        @test haskey(measures, "E_G")
        @test haskey(measures, "E_T")
        @test haskey(measures, "E_O")
    end

    @testset "RWM Simulation" begin
        system = MolecularSystem(centers, radii, 250.0)
        sol = SolvationParams(1.4, 0.3665)
        ol = LinearOverlapParams(0.0, 100.0)
        topo = TopologyParams([1.0, 0.0, 0.0])
        num = NumericalParams(100.0, false)
        scales = EnergyScales(θ_G=1.0, θ_O=1.1, θ_T=0.5)

        prefactors = Energies.get_wb_prefactors(sol.rs, sol.η)
        bounding_radii = Energies.get_bounding_radii(centers, radii, sol.rs)
        single_energies, single_measures = Energies.get_single_subunit_energy_and_measures(
            centers, radii, sol.rs, prefactors, ol.jump, ol.slope, num.delaunay_eps
        )
        precomputed = Precomputed(bounding_radii, single_energies, single_measures)
        bol_check = x -> Energies.are_bounding_spheres_overlapping(x, 1, 2, bounding_radii)

        # Energy function closure
        energy_fn = x -> Energies.calculate_combined_energy(
            x, system, sol, ol, topo, num, scales, precomputed, bol_check
        )

        # Perturbation function
        σ_r, σ_t = 0.2, 1.25
        perturbation_fn = x -> Utilities.perturb_single_specified(x, σ_r, σ_t; specified_index=2)

        # Initial state
        x_init = [
            (QuatRotation(1.0, 0.0, 0.0, 0.0), [50.0, 50.0, 50.0]),
            (QuatRotation(1.0, 0.0, 0.0, 0.0), [150.0, 150.0, 150.0])
        ]

        # Output storage
        output = Algorithms.SimulationOutput()

        # Run simulation
        β = 1.0 / 1.3  # temperature = 1.3
        rwm = Algorithms.RandomWalkMetropolis(energy_fn, perturbation_fn, β)
        Algorithms.simulate!(rwm, deepcopy(x_init), 60.0, 5, output)

        @test output.total_step_attempts == 5
        @test length(output.E_total) > 0
        @test haskey(output.measures, "E_G")

        # Test save/load round-trip via to_dict
        test_file = joinpath(test_dir, "integration_test.jld2")
        output_dict = Algorithms.to_dict(output)
        input_data = Dict(
            "molecule_ids" => molecule_ids,
            "centers" => centers, "radii" => radii,
            "rs" => sol.rs, "η" => sol.η, "bounds" => system.bounds,
            "target_iterations" => 5, "T" => 1.3,
            "wall_clock_limit_minutes" => 60.0,
            "x_init" => x_init,
        )
        @save test_file input=input_data output=output_dict

        @test isfile(test_file)

        # Load and verify dict round-trip
        loaded = JLD2.load(test_file)
        @test loaded["input"]["molecule_ids"] == molecule_ids
        @test loaded["output"]["total_step_attempts"][1] == 5

        # Verify SimulationOutput reconstruction from dict
        reconstructed = Algorithms.SimulationOutput(loaded["output"])
        @test reconstructed.total_step_attempts == 5
        @test reconstructed.E_total == output.E_total
    end

    # Cleanup
    if isdir(test_dir)
        rm(test_dir, recursive=true)
    end
end  # run_integration_tests

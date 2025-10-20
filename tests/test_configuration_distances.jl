function run_configuration_distance_tests()
    @testset verbose = true "Configuration Distances" begin
        @testset verbose = true "Integration" begin
            test_configuration_distance_calls()
        end
    end
end

function test_configuration_distance_calls()
    tc = [0.0 1.0; 0.0 0.0; 0.0 0.0]
    n_mol = 2
    state_a = fill(0.0, 6 * n_mol)
    state_a[4] = 5.0
    state_b = fill(0.0, 6 * n_mol)
    state_b[4] = 5.5

    @test MorphoMol.get_configuration_distance(tc, state_a, state_b, n_mol) == MorphoMol.get_configuration_distance(tc, tc, state_a, state_b) 
    @test MorphoMol.get_configuration_distance(tc, state_a, state_b, n_mol) == 0.5

    n_mol = 3
    state_a = fill(0.0, 6 * n_mol)
    state_a[4] = 5.0
    state_b = fill(0.0, 6 * n_mol)
    state_b[4] = 5.5

    @test MorphoMol.get_configuration_distance(tc, state_a, state_b, n_mol) == MorphoMol.sum_of_permutation(tc, tc, state_a, state_b, [1,2,3], [1,2,3], MorphoMol.get_theta_of_pair) 
    @test MorphoMol.get_configuration_distance(tc, state_a, state_b, n_mol) ≈ 1.0/3.0

    n_mol = 4
    state_a = fill(0.0, 6 * n_mol)
    state_a[4] = 5.0
    state_b = fill(0.0, 6 * n_mol)
    state_b[4] = 5.5

    @test MorphoMol.get_configuration_distance(tc, state_a, state_b, n_mol) ≈ 0.25
end

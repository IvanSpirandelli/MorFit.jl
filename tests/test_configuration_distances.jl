function run_configuration_distance_tests()
    @testset verbose = true "Configuration Distances" begin
        @testset verbose = true "RMSD Align One of Pair" begin
            test_rmsd_align_one_of_pair()
        end
        @testset verbose = true "Global RMSD Target-Inhibitor" begin
            test_global_rmsd_target_inhibitor()
        end
        @testset verbose = true "Fixed Target-Inhibitor RMSD" begin
            test_fixed_target_inhibitor_rmsd()
        end
        @testset verbose = true "Metric to Ground Truth" begin
            test_metric_to_ground_truth()
        end
    end
end

function test_rmsd_align_one_of_pair()
    # Tests get_rmsd_align_one_of_pair with two molecule pairs
    # Should verify RMSD calculation when aligning one molecule of a pair
    @assert false "Requires Implementation."
end

function test_global_rmsd_target_inhibitor()
    # Tests get_global_rmsd_of_target_inhibitor_pair
    # Should verify global RMSD when both target and inhibitor are aligned together
    @assert false "Requires Implementation."
end

function test_fixed_target_inhibitor_rmsd()
    # Tests get_rmsd_for_fixed_target_inhibitor_pair
    # Should verify RMSD when target is fixed and inhibitor follows
    @assert false "Requires Implementation."
end

function test_metric_to_ground_truth()
    # Tests get_metric_to_ground_truth for n_mol=2 and n_mol=3 cases
    # Should verify comparison against experimental assembly data
    @assert false "Requires Implementation."
end

function run_energy_call_tests()
    @testset verbose = true "Energy Calls" begin
        @testset verbose = true "Solvation Free Energy" begin
            test_solvation_free_energy()
        end
        @testset verbose = true "Connected Component Energy (n=2)" begin
            test_cc_energy_two_molecules()
        end
        @testset verbose = true "Connected Component Energy (n>2)" begin
            test_cc_energy_multi_molecules()
        end
        @testset verbose = true "Combined Potential (Fsol + Twasp)" begin
            test_combined_potential()
        end
        @testset verbose = true "Bounding Sphere Overlap" begin
            test_bounding_sphere_overlap()
        end
    end
end

# Helper to create test input with unified template format
function _get_test_input(n_mol::Int)
    mol_type = "6r7m"
    raw_centers = MorFit.TEMPLATES[mol_type]["template_centers"]
    raw_radii = MorFit.TEMPLATES[mol_type]["template_radii"]
    n_atoms = length(raw_centers) ÷ 3
    single_template = reshape(raw_centers, (3, n_atoms))

    template_centers = [single_template for _ in 1:n_mol]
    template_radii = [raw_radii for _ in 1:n_mol]

    rs = 1.4
    η = 0.3665
    prefactors = MorFit.Energies.get_prefactors(rs, η)

    return Dict(
        "template_centers" => template_centers,
        "template_radii" => template_radii,
        "bounds" => 150.0,
        "rs" => rs,
        "prefactors" => prefactors,
        "overlap_jump" => 0.0,
        "overlap_slope" => 1.1,
        "delaunay_eps" => 100.0,
    )
end

function test_solvation_free_energy()
    # Tests MorFit.Energies.solvation_free_energy_and_measures
    # Should verify energy computation for a simple 2-molecule configuration
    @assert false "Requires Implementation."
end

function test_cc_energy_two_molecules()
    # Tests solvation_free_energy_and_measures_with_bounding_container_check
    # Verifies that cc_fsol matches fsol for n_mol=2 when bounding spheres overlap
    @assert false "Requires Implementation."
end

function test_cc_energy_multi_molecules()
    # Tests connected_component_wise_solvation_free_energy_and_measures
    # Verifies correct energy computation for n_mol>2 with connected components
    @assert false "Requires Implementation."
end

function test_combined_potential()
    # Tests calculate_combined_potential (fsol + twasp interpolated)
    # Verifies: 0.5*(e_fsol + e_twasp) ≈ e_combined for μ=0.5
    @assert false "Requires Implementation."
end

function test_bounding_sphere_overlap()
    # Tests are_bounding_spheres_overlapping and get_bounding_radii
    # Should verify overlap detection works correctly
    @assert false "Requires Implementation."
end

function existing_values_equal(d1, d2)
    for k in keys(d1)
        if k in keys(d2)
            if !(d1[k] ≈ d2[k])
                return false
            end
        end
    end
    return true
end

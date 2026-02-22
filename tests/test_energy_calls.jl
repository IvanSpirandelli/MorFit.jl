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
        @testset verbose = true "Bounding Sphere Overlap" begin
            test_bounding_sphere_overlap()
        end
        @testset verbose = true "Overlap vs Edge Penalty" begin
            test_overlap_vs_edge_penalty()
        end
        @testset verbose = true "Periodic Alpha Complex (3x3x3 grid)" begin
            test_periodic_alpha_complex_cubic_grid()
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
    prefactors = MorFit.Energies.get_wb_prefactors(rs, η)

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
    # Tests calculate_combined_energy for n_mol=2
    # Verifies energy computation when bounding spheres overlap
    @assert false "Requires Implementation."
end

function test_cc_energy_multi_molecules()
    # Tests connected_component_wise_solvation_free_energy_and_measures
    # Verifies correct energy computation for n_mol>2 with connected components
    @assert false "Requires Implementation."
end

function test_bounding_sphere_overlap()
    # Tests are_bounding_spheres_overlapping and get_bounding_radii
    # Should verify overlap detection works correctly
    @assert false "Requires Implementation."
end

function test_overlap_vs_edge_penalty()
    # The two calculations do not perfectly match on all seeds. It is time to move away from AlphaMol.
    Random.seed!(999)

    n_points = 25
    radius = 1.0
    points = [rand(3) for _ in 1:n_points]

    # AlphaMol overlap: each point is its own molecule
    flat_coords = vcat(points...)
    molecule_sizes = fill(1, n_points)
    radii = fill(radius, n_points)
    measures = Energies.get_geometric_measures_and_overlap_value(
        flat_coords, molecule_sizes, radii, 1.4, 0.0, 1.0, 1.0
    )
    E_O_alphamol = measures[5]

    # Edge penalty via alpha complex
    edge_penalty = d -> max(0.0, 2 * radius - d)
    _, edge_measures = Energies.point_cloud_persistence_energy_with_edge_penalty(
        points, [0.0, 0.0, 0.0], edge_penalty
    )
    E_O_edges = edge_measures["E_penalty"]

    @test E_O_alphamol ≈ E_O_edges
end

function test_periodic_alpha_complex_cubic_grid()
    # 27 points on a 3×3×3 cubic grid with periodic box [0,3]³.
    # Under periodicity this is a perfect lattice on T³.
    # All finite features have persistence exactly 0.25:
    #   dim 0: 26 = n-1 features, born at 0.0, die at 0.25 (all edges are length 1, α²=0.25)
    #   dim 1: 52 = 2(n-1) features, born at 0.25, die at 0.5
    #   dim 2: 26 = n-1 features, born at 0.5, die at 0.75
    # The infinite features (β₀=1, β₁=3, β₂=3 for T³) are excluded by include_inf_points=False.

    points = [[Float64(i), Float64(j), Float64(k)] for i in 0:2 for j in 0:2 for k in 0:2]
    box_lower = [0.0, 0.0, 0.0]
    box_upper = [3.0, 3.0, 3.0]

    pdgm = Energies.get_periodic_alpha_shape_persistence_diagram(points, box_lower, box_upper, false)

    n = 27
    @test size(pdgm[1]) == (n - 1, 2)
    @test size(pdgm[2]) == (2 * (n - 1), 2)
    @test size(pdgm[3]) == (n - 1, 2)

    P0 = Energies.get_total_persistence(pdgm[1])
    P1 = Energies.get_total_persistence(pdgm[2])
    P2 = Energies.get_total_persistence(pdgm[3])

    @test P0 ≈ 6.5
    @test P1 ≈ 13.0
    @test P2 ≈ 6.5

    # Verify the convenience function returns consistent results
    λ = [1.0, 1.0, 1.0]
    energy, measures = Energies.point_cloud_periodic_persistence_energy(points, λ, box_lower, box_upper, false)
    @test energy ≈ P0 + P1 + P2
    @test measures["P0s"] ≈ P0
    @test measures["P1s"] ≈ P1
    @test measures["P2s"] ≈ P2
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

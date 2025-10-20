function run_morphometric_approach_tests()
    @testset verbose = true "Ball Union Measures" begin
        @testset verbose = true "Integration" begin
            test_ball_union_measures()
        end
        @testset verbose = true "Derivatives" begin
            test_derivatives()
        end
    end
end

function test_ball_union_measures()
    two_spheres = [0.0, 0.0, 0.0, 2.0, 0.0, 0.0] 
    radii = [1.0, 1.0]
    probe_radius = 0.2
    geometric_measures = MorphoMol.Energies.get_geometric_measures(two_spheres, radii, probe_radius, 1.0)
    @test geometric_measures ≈ [14.191621213816292, 33.17521842190822, 23.540276903519608, 12.566370614359174]

    tetrahedron = [
        [1.0, 0.0, -1/sqrt(2)],
        [-1.0, 0.0, -1/sqrt(2)],
        [0.0, 1.0, 1/sqrt(2)],
        [0.0, -1.0, 1/sqrt(2)]
    ] # edge length = 2.0
    atom_coordinates = [(tetrahedron...)...]
    atom_radii = [1.0, 1.0, 1.0, 1.0]
    inflated_radii = [2.4, 2.4, 2.4, 2.4]
    probe_radius = 1.4
    
    geometric_measures_1 = MorphoMol.Energies.get_geometric_measures(atom_coordinates, atom_radii, probe_radius, 1.0)

    geometric_measures_2, _, _, _, _ = MorphoMol.Energies.get_geometric_measures_with_derivatives(
        atom_coordinates, atom_radii, probe_radius, 1.0
    )   
    @test geometric_measures_1 ≈ geometric_measures_2[1:4]

    geometric_measures_3 = MorphoMol.Energies.get_geometric_measures(atom_coordinates, inflated_radii, 0.0, 1.0)
    @test geometric_measures_1 ≈ geometric_measures_3

    hadwiger_measures_and_overlap_1 = MorphoMol.Energies.get_geometric_measures_and_overlap_value(
        atom_coordinates, 1, atom_radii, probe_radius, 0.0, 0.0, 1.0
    )
    @test geometric_measures_1 ≈ hadwiger_measures_and_overlap_1[1:4]
    @test hadwiger_measures_and_overlap_1[5] ≈ 0.0

    hadwiger_measures_and_overlap_2 = MorphoMol.Energies.get_geometric_measures_and_overlap_value(
        atom_coordinates, 1, [1.1, 1.1, 1.1, 1.1], probe_radius, 1.0, 0.0, 1.0
    )
    @test hadwiger_measures_and_overlap_2[5] ≈ 6.0 

    hadwiger_measures_and_overlap_3 = MorphoMol.Energies.get_geometric_measures_and_overlap_value(
        atom_coordinates, 1, [1.1, 1.1, 1.1, 1.1], probe_radius, 0.0, 1.0, 1.0
    )
    @test hadwiger_measures_and_overlap_3[5] ≈ 1.2 

    hadwiger_measures_and_overlap_4, _, _, _, _ = MorphoMol.Energies.get_geometric_measures_and_overlap_value_with_derivatives(
        atom_coordinates, 1, [1.1, 1.1, 1.1, 1.1], probe_radius, 0.0, 1.0, 1.0
    )
    @test hadwiger_measures_and_overlap_3 ≈ hadwiger_measures_and_overlap_4

    hadwiger_measures_and_overlap_5, _, _, _, _ = MorphoMol.Energies.get_geometric_measures_and_overlap_value_with_derivatives(
        atom_coordinates, 2, [1.1, 1.1, 1.1, 1.1], probe_radius, 1.0, 0.0, 1.0
    )
    @test hadwiger_measures_and_overlap_5[5] ≈ 4.0

    hadwiger_measures_and_overlap_6 = MorphoMol.Energies.get_geometric_measures_and_overlap_value(
        atom_coordinates, 2, [1.1, 1.1, 1.1, 1.1], probe_radius, 1.0, 0.0, 1.0
    )
    @test hadwiger_measures_and_overlap_6[5] ≈ 4.0
end

function test_derivatives()
    coordinates = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0]
    ]
    n = length(coordinates)
    atom_coordinates = [(coordinates...)...]
    atom_radii = [1.0, 1.0, 1.0, 1.0]
    probe_radius = 1.4

    geometric_measures, dvol, dsurf, dmean, dgauss = MorphoMol.Energies.get_geometric_measures_with_derivatives(
        atom_coordinates, atom_radii, probe_radius, 1.0
    )   
    dvol_r = reshape(dvol, (3,n))
    @test dvol_r[:,1] ≈ -dvol_r[:,2]

    dsurf_r = reshape(dsurf, (3,n))
    @test dsurf_r[:,1] ≈ -dsurf_r[:,2]

    hadwiger_measures_and_overlap, dvol, dsurf, dmean, dgauss, dlol = MorphoMol.Energies.get_geometric_measures_and_overlap_value_with_derivatives(
        atom_coordinates, 1, [1.0, 1.0, 1.0, 1.0], probe_radius, 0.0, 1.0, 1.0
    )
    dlol_r = reshape(dlol, (3,n))
    @test dlol_r[:,1] ≈ -dlol_r[:,2]
    @test dlol_r[:,1][1] ≈ 1.0

    coordinates = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.8, 0.0],
        [0.0, 1.8, -1.8]
    ]
    n = length(coordinates)
    atom_coordinates = [(coordinates...)...]
    _, _, _, _, _, dlol = MorphoMol.Energies.get_geometric_measures_and_overlap_value_with_derivatives(
        atom_coordinates, 2, [1.0, 1.0, 1.0, 1.0], probe_radius, 0.0, 1.0, 1.0
    )
    dlol_r = reshape(dlol, (3,n))
    @test dlol_r[:,1] ≈ -dlol_r[:,3]
    @test dlol_r[:,2] ≈ zeros(3)
    @test dlol_r[:,4] ≈ zeros(3)
end
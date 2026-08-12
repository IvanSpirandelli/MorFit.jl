using Rotations
using LinearAlgebra

"""
Test suite for get_realization function.

Tests all output formats and validates geometric consistency.
"""

# Deterministic test fixtures
function get_test_fixtures()
    # Simple template: 2 molecules, each with different atom counts
    # Molecule 1: triangle in xy-plane (3 atoms)
    # Molecule 2: line segment (2 atoms)
    template_centers = [
        [0.0 1.0 0.5;    # 3x3 matrix for mol 1
         0.0 0.0 0.866;
         0.0 0.0 0.0],
        [0.0 1.0;        # 3x2 matrix for mol 2
         0.0 0.0;
         0.0 0.0]
    ]

    # State: identity rotation + [1,2,3] translation for mol 1
    #        90° z-rotation + [4,5,6] translation for mol 2
    R1 = QuatRotation(1.0, 0.0, 0.0, 0.0)  # identity
    t1 = [1.0, 2.0, 3.0]

    R2 = QuatRotation(cos(π/4), 0.0, 0.0, sin(π/4))  # 90° around z
    t2 = [4.0, 5.0, 6.0]

    x = [(R1, t1), (R2, t2)]

    return x, template_centers
end

# Pre-computed expected values:
# Mol 1 (identity rotation + [1,2,3]): [0,0,0]→[1,2,3], [1,0,0]→[2,2,3], [0.5,0.866,0]→[1.5,2.866,3]
# Mol 2 (90° z-rot + [4,5,6]): [0,0,0]→[4,5,6], [1,0,0]→[0,1,0]+[4,5,6]=[4,6,6]

function run_realization_tests()
    x, template_centers = get_test_fixtures()

    @testset "format=:matrix (default)" begin
        result = MorFit.Utilities.get_realization(x, template_centers)
        result_explicit = MorFit.Utilities.get_realization(x, template_centers, format=:matrix)

        @test result == result_explicit  # default is :matrix
        @test result isa Vector{Matrix{Float64}}
        @test length(result) == 2
        @test size(result[1]) == (3, 3)  # 3 coords × 3 points
        @test size(result[2]) == (3, 2)  # 3 coords × 2 points

        # Molecule 1 matrix (columns are points)
        @test result[1][:, 1] ≈ [1.0, 2.0, 3.0] atol=1e-10
        @test result[1][:, 2] ≈ [2.0, 2.0, 3.0] atol=1e-10
        @test result[1][:, 3] ≈ [1.5, 2.866, 3.0] atol=1e-3

        # Molecule 2 matrix
        @test result[2][:, 1] ≈ [4.0, 5.0, 6.0] atol=1e-10
        @test result[2][:, 2] ≈ [4.0, 6.0, 6.0] atol=1e-10
    end

    @testset "format=:flat" begin
        result = MorFit.Utilities.get_realization(x, template_centers, format=:flat)

        @test result isa Vector{Float64}
        @test length(result) == 15  # 5 points × 3 coords

        # First molecule points (identity rotation)
        @test result[1] ≈ 1.0 atol=1e-10   # x of point 1
        @test result[2] ≈ 2.0 atol=1e-10   # y of point 1
        @test result[3] ≈ 3.0 atol=1e-10   # z of point 1
        @test result[4] ≈ 2.0 atol=1e-10   # x of point 2
        @test result[5] ≈ 2.0 atol=1e-10   # y of point 2
        @test result[6] ≈ 3.0 atol=1e-10   # z of point 2
        @test result[7] ≈ 1.5 atol=1e-10   # x of point 3
        @test result[8] ≈ 2.866 atol=1e-3  # y of point 3
        @test result[9] ≈ 3.0 atol=1e-10   # z of point 3

        # Second molecule points (90° z-rotation)
        @test result[10] ≈ 4.0 atol=1e-10  # x of point 4
        @test result[11] ≈ 5.0 atol=1e-10  # y of point 4
        @test result[12] ≈ 6.0 atol=1e-10  # z of point 4
        @test result[13] ≈ 4.0 atol=1e-10  # x of point 5
        @test result[14] ≈ 6.0 atol=1e-10  # y of point 5
        @test result[15] ≈ 6.0 atol=1e-10  # z of point 5
    end

    @testset "format=:points" begin
        result = MorFit.Utilities.get_realization(x, template_centers, format=:points)

        @test result isa Vector{Vector{Float64}}
        @test length(result) == 5  # 5 points total

        @test result[1] ≈ [1.0, 2.0, 3.0] atol=1e-10
        @test result[2] ≈ [2.0, 2.0, 3.0] atol=1e-10
        @test result[3] ≈ [1.5, 2.866, 3.0] atol=1e-3
        @test result[4] ≈ [4.0, 5.0, 6.0] atol=1e-10
        @test result[5] ≈ [4.0, 6.0, 6.0] atol=1e-10
    end

    @testset "format=:points_per_mol" begin
        result = MorFit.Utilities.get_realization(x, template_centers, format=:points_per_mol)

        @test result isa Vector{Vector{Vector{Float64}}}
        @test length(result) == 2  # 2 molecules
        @test length(result[1]) == 3  # mol 1 has 3 points
        @test length(result[2]) == 2  # mol 2 has 2 points

        # Molecule 1 points
        @test result[1][1] ≈ [1.0, 2.0, 3.0] atol=1e-10
        @test result[1][2] ≈ [2.0, 2.0, 3.0] atol=1e-10
        @test result[1][3] ≈ [1.5, 2.866, 3.0] atol=1e-3

        # Molecule 2 points
        @test result[2][1] ≈ [4.0, 5.0, 6.0] atol=1e-10
        @test result[2][2] ≈ [4.0, 6.0, 6.0] atol=1e-10
    end

    @testset "format=:point3f" begin
        result = MorFit.Utilities.get_realization(x, template_centers, format=:point3f)

        @test result isa Vector{Point3f}
        @test length(result) == 5

        @test result[1] ≈ Point3f(1.0, 2.0, 3.0) atol=1e-5
        @test result[2] ≈ Point3f(2.0, 2.0, 3.0) atol=1e-5
        @test result[3] ≈ Point3f(1.5, 2.866, 3.0) atol=1e-3
        @test result[4] ≈ Point3f(4.0, 5.0, 6.0) atol=1e-5
        @test result[5] ≈ Point3f(4.0, 6.0, 6.0) atol=1e-5
    end

    @testset "format=:point3f_per_mol" begin
        result = MorFit.Utilities.get_realization(x, template_centers, format=:point3f_per_mol)

        @test result isa Vector{Vector{Point3f}}
        @test length(result) == 2
        @test length(result[1]) == 3
        @test length(result[2]) == 2

        # Molecule 1 points
        @test result[1][1] ≈ Point3f(1.0, 2.0, 3.0) atol=1e-5
        @test result[1][2] ≈ Point3f(2.0, 2.0, 3.0) atol=1e-5
        @test result[1][3] ≈ Point3f(1.5, 2.866, 3.0) atol=1e-3

        # Molecule 2 points
        @test result[2][1] ≈ Point3f(4.0, 5.0, 6.0) atol=1e-5
        @test result[2][2] ≈ Point3f(4.0, 6.0, 6.0) atol=1e-5
    end

    @testset "Invalid format" begin
        @test_throws ArgumentError MorFit.Utilities.get_realization(x, template_centers, format=:invalid)
    end

    @testset "Cross-format consistency" begin
        # All formats should produce equivalent geometric data
        matrices = MorFit.Utilities.get_realization(x, template_centers, format=:matrix)
        flat = MorFit.Utilities.get_realization(x, template_centers, format=:flat)
        points = MorFit.Utilities.get_realization(x, template_centers, format=:points)
        points_per_mol = MorFit.Utilities.get_realization(x, template_centers, format=:points_per_mol)
        point3f = MorFit.Utilities.get_realization(x, template_centers, format=:point3f)

        # flat vs points
        for (i, pt) in enumerate(points)
            @test flat[(i-1)*3+1:(i-1)*3+3] ≈ pt atol=1e-10
        end

        # points vs points_per_mol (flattened)
        flattened_per_mol = vcat(points_per_mol...)
        @test length(flattened_per_mol) == length(points)
        for (i, pt) in enumerate(flattened_per_mol)
            @test pt ≈ points[i] atol=1e-10
        end

        # points vs point3f (allowing Float32 precision loss)
        for (i, p3f) in enumerate(point3f)
            @test Vector{Float64}(p3f) ≈ points[i] atol=1e-5
        end

        # matrices vs points_per_mol
        for (mol_idx, mat) in enumerate(matrices)
            for (pt_idx, col) in enumerate(eachcol(mat))
                @test Vector{Float64}(col) ≈ points_per_mol[mol_idx][pt_idx] atol=1e-10
            end
        end
    end
end

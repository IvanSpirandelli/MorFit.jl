module Tests

export run, run_morphometric_approach_tests, run_realization_tests, run_integration_tests

using Test
using GeometryBasics
using JLD2
using Rotations

using MorFit
using MorFit.Energies
using MorFit.Algorithms
using MorFit.Utilities

include("test_morphometric_approach.jl")
include("test_energy_calls.jl")
include("test_configuration_distances.jl")
include("test_realizations.jl")
include("test_integration.jl")

function run()
    @testset verbose = true "Morphometric Approach Tests" begin
        run_morphometric_approach_tests()
    end
    @testset verbose = true "Energy Call Tests" begin
        run_energy_call_tests()
    end
    @testset verbose = true "Configuration Distances Tests" begin
        run_configuration_distance_tests()
    end
    @testset verbose = true "Realization Tests" begin
        run_realization_tests()
    end
    @testset verbose = true "Integration Tests" begin
        run_integration_tests()
    end
end

end

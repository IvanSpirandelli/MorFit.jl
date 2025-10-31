module MorFit
using JLD2, Rotations
# Modules
include("modules/Utilities/src/Utilities.jl")
include("modules/Algorithms/src/Algorithms.jl")
include("modules/Energies/src/Energies.jl")

# Templates
include("templates/target_and_inhibitor_templates.jl")
include("templates/experimental_assembly.jl")
include("templates/asymmetric_unit_templates.jl")

@load joinpath(@__DIR__, "templates/protein_ligand_data.jld2") protein_ligand_data
const PROTEIN_LIGAND_DATA = protein_ligand_data

# Setup
include("simulation_setup.jl")

# Tests
include("../tests/Tests.jl")

end #module MorFit
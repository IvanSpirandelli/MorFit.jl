module MorFit
using Rotations
# Modules
include("modules/Utilities/src/Utilities.jl")
include("modules/Algorithms/src/Algorithms.jl")
include("modules/Energies/src/Energies.jl")

# Templates
include("templates/target_and_inhibitor_templates.jl")
include("templates/experimental_assembly.jl")
include("templates/asymmetric_unit_templates.jl")

# Setup
include("simulation_setup.jl")

# Tests
include("../tests/Tests.jl")

end #module MorFit
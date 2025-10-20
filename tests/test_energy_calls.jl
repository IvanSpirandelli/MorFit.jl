function run_energy_call_tests()
    @testset verbose = true "Energy Calls" begin
        @testset verbose = true "Integration" begin
            test_energy_calls()
        end
        @testset verbose = true "Previous Errors" begin 
            test_previous_error_calls()
        end
    end
end

function test_previous_error_calls()
    x_flat = [1.8719770104797722, -0.8054259756520944, 2.264161853094224, 64.04252455563663, 47.88893314398397, 31.29679809443243, 3.7196395887896667, 1.3234238947069903, 0.5481448529951869, 74.26445359430221, 43.84510418286682, 80.91055320912551]
    x = MorFit.Utilities.convert_flat_state_to_tuples(x_flat)
    input = get_input(2)
    fsol_e, _, cc_fsol_e, _ = get_energy_calls(input)

    e1, m1 = fsol_e(x)
    @assert m1["OLs"] > 0.0
    e2, m2 = cc_fsol_e(x)
    are_overlapping = MorFit.Energies.are_bounding_spheres_overlapping(x, 1, 2, MorFit.Energies.get_bounding_radii(input["template_centers"], input["template_radii"], input["rs"]))
    @assert are_overlapping
    @assert e1 == e2 
    @assert existing_values_equal(m1, m2)
end

function test_energy_calls()
    n_mol = 2
    input = get_input(n_mol)
    fsol_e, twasp_e, cc_fsol_e, cc_fsol_twasp_e = get_energy_calls(input)
    x = MorFit.Utilities.get_initial_state(n_mol, input["bounds"])

    e1, m1 = fsol_e(x)
    e2, m2 = twasp_e(x)
    e3, m3 = cc_fsol_twasp_e(x)
    @test 0.5*(e1 + e2) ≈ e3
    @test existing_values_equal(m1, m2)
    @test existing_values_equal(m1, m3)
    @test existing_values_equal(m2, m3)

    e4, m4 = cc_fsol_e(x)
    @test e1 ≈ e4
    @test existing_values_equal(m1, m4)


    n_mol = 3
    input = get_input(n_mol)
    fsol_e, twasp_e, cc_fsol_e, cc_fsol_twasp_e = get_energy_calls(input)
    x = MorFit.Utilities.get_initial_state(n_mol, input["bounds"])

    e1, m1 = fsol_e(x)
    e2, m2 = twasp_e(x)

    bol_nmol_l = (x, id1, id2) -> MorFit.Energies.are_bounding_spheres_overlapping(x, id1, id2, MorFit.Energies.get_bounding_radii(input["template_centers"], input["template_radii"], input["rs"]))
    icc = MorFit.Energies.get_initial_connected_component_energies(x, input["template_centers"], input["template_radii"], input["rs"], input["prefactors"], input["overlap_jump"], input["overlap_slope"], input["delaunay_eps"], bol_nmol_l)
    
    e3, m3 = cc_fsol_twasp_e(icc, 1, x)
    @test 0.5*(e1 + e2) ≈ e3
    @test existing_values_equal(m1, m2)
    @test existing_values_equal(m1, m3)
    @test existing_values_equal(m2, m3)
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

function get_input(n_mol::Int)
    bounds = 150.0
    delaunay_eps = 100.0
    overlap_jump = 0.0
    overlap_slope = 1.1
    rs = 1.4
    η = 0.3665
    temperature = 0.0

    comment = "tasp_icosahedron"
    comment = replace(comment, " " => "_")

    T_search_runs = 12
    T_search_time = 15.0

    σ_r = 0.3
    σ_t = 1.25

    μ = 0.5

    simulation_time_minutes = 8 * 60.0
    algorithm = "sa"
    sa_level = "[1.0,0.7,0.5,0.3,0.1,0.0]"
    energy = "_"
    perturbation = "single_random"
    initialization = "random"

    persistence_weights = [1.0, -1.0, -1.0]

    mol_type = "6r7m"
    template_centers = MorFit.TEMPLATES[mol_type]["template_centers"]
    template_radii = MorFit.TEMPLATES[mol_type]["template_radii"]
    
    prefactors = MorFit.Energies.get_prefactors(rs, η)
    x = MorFit.Utilities.get_initial_state(n_mol, bounds)

    input = Dict(
            "algorithm" => algorithm,
            "sa_level" => sa_level,
            "energy" => energy,
            "perturbation" => perturbation,
            "initialization" => initialization,
            "mol_type" => mol_type,
            "template_centers" => template_centers,
            "template_radii" => template_radii,
            "n_mol" => n_mol,
            "x_init" => x,
            "comment" => comment,
            "bounds" => bounds,
            "rs" => rs,
            "η" => η,
            "μ" => μ,
            "prefactors" => prefactors,
            "σ_r" => σ_r,
            "σ_t" => σ_t,
            "T_search_runs" => T_search_runs,
            "T_search_time" => T_search_time,
            "T" => temperature,
            "persistence_weights" => persistence_weights,
            "overlap_jump" => overlap_jump,
            "overlap_slope" => overlap_slope,
            "delaunay_eps" => delaunay_eps,
            "simulation_time_minutes" => simulation_time_minutes,
            "exact_delaunay" => false
        )
    return input
end

function get_energy_calls(input)
    input["energy"] = "fsol"
    fsol_e = MorFit.get_energy(input)
    input["energy"] = "twasp"
    twasp_e = MorFit.get_energy(input)
    input["energy"] = "cc_fsol"
    cc_fsol_e = MorFit.get_energy(input)
    input["energy"] = "cc_fsol_twasp_interpolated"
    cc_fsol_twasp_e = MorFit.get_energy(input)

    return fsol_e, twasp_e, cc_fsol_e, cc_fsol_twasp_e
end
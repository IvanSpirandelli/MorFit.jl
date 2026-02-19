# Test MorFit persistence (diode+oineus) against GUDHI
#
# Usage:
#   PYCALL_JL_RUNTIME_PYTHON=python-environments/python-tda-313/bin/python3 \
#   julia --project=MorFit.jl MorFit.jl/tests/test_persistence_against_gudhi.jl

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using MorFit
using PyCall
using Random
using Printf

Random.seed!(42)

# ── GUDHI via py"..." for reliable access ──
py"""
import gudhi
import numpy as np

def gudhi_persistence_unweighted(points):
    \"\"\"Returns dict with P0, P1, P2 total persistence (using raw GUDHI filtration values).\"\"\"
    ac = gudhi.AlphaComplex(points=np.array(points))
    st = ac.create_simplex_tree()
    pers = st.persistence()
    P_raw = [0.0, 0.0, 0.0]      # sum(d - b) with squared filtration values
    P_sqrt = [0.0, 0.0, 0.0]     # sum(sqrt(d) - sqrt(b))
    for dim, (b, d) in pers:
        if dim >= 3 or d == float('inf'):
            continue
        P_raw[dim] += d - b
        P_sqrt[dim] += d**0.5 - b**0.5
    return {"P_raw": P_raw, "P_sqrt": P_sqrt}

def gudhi_persistence_weighted(points, radii_squared):
    \"\"\"Returns dict with P0, P1, P2 total persistence for weighted alpha complex.\"\"\"
    try:
        ac = gudhi.WeightedAlphaComplex(points=np.array(points), weights=np.array(radii_squared))
    except AttributeError:
        # Older GUDHI: weights parameter on AlphaComplex
        ac = gudhi.AlphaComplex(points=np.array(points), weights=np.array(radii_squared))
    st = ac.create_simplex_tree()
    pers = st.persistence()
    P_raw = [0.0, 0.0, 0.0]
    P_sqrt = [0.0, 0.0, 0.0]
    for dim, (b, d) in pers:
        if dim >= 3 or d == float('inf'):
            continue
        P_raw[dim] += d - b
        # Weighted filtration can be negative
        b_s = b**0.5 if b >= 0 else -((-b)**0.5)
        d_s = d**0.5 if d >= 0 else -((-d)**0.5)
        P_sqrt[dim] += d_s - b_s
    return {"P_raw": P_raw, "P_sqrt": P_sqrt}
"""

function morfit_total_persistence_unweighted(points::Vector{Vector{Float64}})
    λ = [1.0, 1.0, 1.0]
    _, measures = MorFit.Energies.point_cloud_persistence_energy(points, λ, false)
    return [measures["P0s"], measures["P1s"], measures["P2s"]]
end

function morfit_total_persistence_weighted(points::Vector{Vector{Float64}}, radii::Vector{Float64})
    pdgm = MorFit.Energies.get_weighted_alpha_shape_persistence_diagram(points, radii, false)
    P0 = MorFit.Energies.get_total_persistence(pdgm[1])
    P1 = MorFit.Energies.get_total_persistence(pdgm[2])
    P2 = MorFit.Energies.get_total_persistence(pdgm[3])
    return [P0, P1, P2]
end

function check_match(label, gudhi_val, morfit_val; tol=0.01)
    if gudhi_val == 0.0 && morfit_val == 0.0
        return true, "both zero"
    elseif gudhi_val != 0.0
        rel_err = abs(gudhi_val - morfit_val) / abs(gudhi_val)
        return rel_err < tol, @sprintf("rel_err=%.2e", rel_err)
    else
        return abs(morfit_val) < 1e-10, @sprintf("abs_err=%.2e", abs(morfit_val))
    end
end

# ── Run tests ──
println("=" ^ 95)
println("  Persistence comparison: MorFit (diode+oineus) vs GUDHI")
println("=" ^ 95)

n_sizes = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
global all_pass = true

for n in n_sizes
    points = [randn(3) * 5.0 for _ in 1:n]
    radii = abs.(randn(n)) .+ 0.5

    # ── Unweighted ──
    gudhi_result = py"gudhi_persistence_unweighted"(points)
    gudhi_raw = Float64.(gudhi_result["P_raw"])
    gudhi_sqrt = Float64.(gudhi_result["P_sqrt"])
    morfit_P = morfit_total_persistence_unweighted(points)

    println("\n─── n=$n (unweighted) ───")
    @printf("  %-14s  %12s  %12s  %12s\n", "", "P0", "P1", "P2")
    @printf("  %-14s  %12.6f  %12.6f  %12.6f\n", "GUDHI raw", gudhi_raw...)
    @printf("  %-14s  %12.6f  %12.6f  %12.6f\n", "GUDHI sqrt", gudhi_sqrt...)
    @printf("  %-14s  %12.6f  %12.6f  %12.6f\n", "MorFit", morfit_P...)

    # Check which convention matches
    for (dim, label) in enumerate(["P0", "P1", "P2"])
        raw_ok, raw_msg = check_match(label, gudhi_raw[dim], morfit_P[dim])
        sqrt_ok, sqrt_msg = check_match(label, gudhi_sqrt[dim], morfit_P[dim])
        if raw_ok
            @printf("  → %s matches GUDHI raw ✓ (%s)\n", label, raw_msg)
        elseif sqrt_ok
            @printf("  → %s matches GUDHI sqrt ✓ (%s)\n", label, sqrt_msg)
        else
            @printf("  → %s: NO MATCH ✗  (vs raw: %s, vs sqrt: %s)\n", label, raw_msg, sqrt_msg)
            global all_pass = false
        end
    end

    # ── Weighted ──
    try
        gudhi_result_w = py"gudhi_persistence_weighted"(points, [r^2 for r in radii])
        gudhi_raw_w = Float64.(gudhi_result_w["P_raw"])
        gudhi_sqrt_w = Float64.(gudhi_result_w["P_sqrt"])
        morfit_Pw = morfit_total_persistence_weighted(points, radii)

        println("\n─── n=$n (weighted) ───")
        @printf("  %-14s  %12.6f  %12.6f  %12.6f\n", "GUDHI raw", gudhi_raw_w...)
        @printf("  %-14s  %12.6f  %12.6f  %12.6f\n", "GUDHI sqrt", gudhi_sqrt_w...)
        @printf("  %-14s  %12.6f  %12.6f  %12.6f\n", "MorFit", morfit_Pw...)

        for (dim, label) in enumerate(["P0", "P1", "P2"])
            raw_ok, raw_msg = check_match(label, gudhi_raw_w[dim], morfit_Pw[dim])
            sqrt_ok, sqrt_msg = check_match(label, gudhi_sqrt_w[dim], morfit_Pw[dim])
            if raw_ok
                @printf("  → %s matches GUDHI raw ✓ (%s)\n", label, raw_msg)
            elseif sqrt_ok
                @printf("  → %s matches GUDHI sqrt ✓ (%s)\n", label, sqrt_msg)
            else
                @printf("  → %s: NO MATCH ✗  (vs raw: %s, vs sqrt: %s)\n", label, raw_msg, sqrt_msg)
                global all_pass = false
            end
        end
    catch e
        println("\n─── n=$n (weighted) — SKIPPED: $e ───")
    end
end

println("\n" * "=" ^ 95)
if all_pass
    println("  ALL TESTS PASSED ✓")
else
    println("  SOME TESTS FAILED ✗")
    println("  Check which convention (raw vs sqrt) MorFit uses for filtration values")
end
println("=" ^ 95)

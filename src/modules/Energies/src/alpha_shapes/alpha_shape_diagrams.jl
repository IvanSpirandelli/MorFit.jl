#=============================================================================
# Alpha Shape Persistence Diagrams via Python (oineus + diode)
#
# Python functions are defined in _init_python_functions(), called from
# Energies.__init__() so they are available at runtime (not just precompilation).
=============================================================================#

const _python_initialized = Ref(false)

function _init_python_functions()
    _python_initialized[] && return
    py"""
import oineus as oin
import numpy as np
import diode

def _compute_persistence_diagram(simplices):
    fil = oin.Filtration([oin.Simplex(s[0], s[1]) for s in simplices])
    dcmp = oin.Decomposition(fil, True)
    params = oin.ReductionParams()
    dcmp.reduce(params)
    return dcmp.diagram(fil, include_inf_points=False)

def _alpha_shape_diagram(points, exact):
    simplices = diode.fill_alpha_shapes(np.asarray(points), exact=exact)
    return _compute_persistence_diagram(simplices)

def _weighted_alpha_shape_diagram(points, exact):
    simplices = diode.fill_weighted_alpha_shapes(np.asarray(points), exact=exact)
    return _compute_persistence_diagram(simplices)

def _alpha_shape_diagram_and_edges(points, exact):
    simplices = diode.fill_alpha_shapes(np.asarray(points), exact=exact)
    dgm = _compute_persistence_diagram(simplices)
    edges = [tuple(s[0]) for s in simplices if len(s[0]) == 2]
    return dgm, edges

def _periodic_alpha_shape_diagram(points, box_lower, box_upper, exact):
    box_lower = [float(x) for x in box_lower]
    box_upper = [float(x) for x in box_upper]
    simplices = diode.fill_periodic_alpha_shapes(
        np.asarray(points), exact=exact, **{'from': box_lower, 'to': box_upper}
    )
    return _compute_persistence_diagram(simplices)
"""
    _python_initialized[] = true
end

_extract_diagrams(result) = [result[1], result[2], result[3]]

function get_alpha_shape_persistence_diagram(points, exact=false)
    _extract_diagrams(py"_alpha_shape_diagram"(points, exact))
end

function get_weighted_alpha_shape_persistence_diagram(points, radii, exact=false)
    weighted_points = [[p[1], p[2], p[3], r^2] for (p, r) in zip(points, radii)]
    _extract_diagrams(py"_weighted_alpha_shape_diagram"(weighted_points, exact))
end

function get_alpha_shape_persistence_diagram_and_edges(points, exact=false)
    result = py"_alpha_shape_diagram_and_edges"(points, exact)
    dgm = result[1]
    edges = result[2]
    diagrams = [dgm[1], dgm[2], dgm[3]]
    edge_tuples = [(e[1], e[2]) for e in edges]
    return diagrams, edge_tuples
end

function get_periodic_alpha_shape_persistence_diagram(points, box_lower, box_upper, exact=false)
    _extract_diagrams(py"_periodic_alpha_shape_diagram"(points, box_lower, box_upper, exact))
end

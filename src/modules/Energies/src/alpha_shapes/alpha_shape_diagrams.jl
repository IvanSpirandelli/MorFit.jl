function get_alpha_shape_persistence_diagram(points, exact = false)
    py"""
    import oineus as oin
    import numpy as np
    import diode

    def get_alpha_shape_persistence_diagram(points, exact):
        points = np.asarray(points)
        simplices = diode.fill_alpha_shapes(points, exact=exact)
        fil = oin.Filtration([oin.Simplex(s[0], s[1]) for s in simplices])

        dcmp = oin.Decomposition(fil, True)
        params = oin.ReductionParams()
        dcmp.reduce(params)
        dgm = dcmp.diagram(fil, include_inf_points=False)
        return dgm
    """
    asds = py"get_alpha_shape_persistence_diagram"(points, exact)
    [asds[1], asds[2], asds[3]]
end

function get_weighted_alpha_shape_persistence_diagram(points, radii, exact = false)
    py"""
    import oineus as oin
    import numpy as np
    import diode

    def get_weighted_alpha_shape_persistence_diagram(points, exact):
        points = np.asarray(points)
        simplices = diode.fill_weighted_alpha_shapes(points, exact=exact)
        fil = oin.Filtration([oin.Simplex(s[0], s[1]) for s in simplices])

        dcmp = oin.Decomposition(fil, True)
        params = oin.ReductionParams()
        dcmp.reduce(params)
        dgm = dcmp.diagram(fil, include_inf_points=False)
        return dgm
    """
    weighted_points = [[p[1], p[2], p[3], r^2] for (p,r) in zip(points, radii)]
    wasds = py"get_weighted_alpha_shape_persistence_diagram"(weighted_points, exact)
    [wasds[1], wasds[2], wasds[3]]
end

function get_alpha_shape_persistence_diagram_and_edges(points, exact = false)
    py"""
    import oineus as oin
    import numpy as np
    import diode

    def get_alpha_shape_persistence_diagram_and_edges(points, exact):
        points = np.asarray(points)
        simplices = diode.fill_alpha_shapes(points, exact=exact)
        fil = oin.Filtration([oin.Simplex(s[0], s[1]) for s in simplices])

        dcmp = oin.Decomposition(fil, True)
        params = oin.ReductionParams()
        dcmp.reduce(params)
        dgm = dcmp.diagram(fil, include_inf_points=False)

        edges = [tuple(s[0]) for s in simplices if len(s[0]) == 2]
        return dgm, edges
    """
    result = py"get_alpha_shape_persistence_diagram_and_edges"(points, exact)
    dgm = result[1]
    edges = result[2]
    diagrams = [dgm[1], dgm[2], dgm[3]]
    edge_tuples = [(e[1], e[2]) for e in edges]
    return diagrams, edge_tuples
end

function get_weighted_alpha_shape_persistence_diagram_and_edges(points, radii, exact = false)
    py"""
    import oineus as oin
    import numpy as np
    import diode

    def get_weighted_alpha_shape_persistence_diagram_and_edges(points, exact):
        points = np.asarray(points)
        simplices = diode.fill_weighted_alpha_shapes(points, exact=exact)
        fil = oin.Filtration([oin.Simplex(s[0], s[1]) for s in simplices])

        dcmp = oin.Decomposition(fil, True)
        params = oin.ReductionParams()
        dcmp.reduce(params)
        dgm = dcmp.diagram(fil, include_inf_points=False)

        edges = [tuple(s[0]) for s in simplices if len(s[0]) == 2]
        return dgm, edges
    """
    weighted_points = [[p[1], p[2], p[3], r^2] for (p,r) in zip(points, radii)]
    result = py"get_weighted_alpha_shape_persistence_diagram_and_edges"(weighted_points, exact)
    dgm = result[1]
    edges = result[2]
    diagrams = [dgm[1], dgm[2], dgm[3]]
    edge_tuples = [(e[1], e[2]) for e in edges]
    return diagrams, edge_tuples
end

function debug_alpha_shape(points)
    py"""
    import oineus as oin
    import numpy as np
    import diode

    def debug_alpha_shape(points):
        points = np.asarray(points)
        np.save('case2.npy', points)
        simplices = diode.fill_alpha_shapes(points)
        fil = oin.Filtration([oin.Simplex(s[0], s[1]) for s in simplices])
        dcmp = oin.Decomposition(fil, True)
        params = oin.ReductionParams()
        dcmp.reduce(params)
        dgm = dcmp.diagram(fil, include_inf_points=False)
        return simplices, fil, dgm
    """
    py"debug_alpha_shape"(points)
end
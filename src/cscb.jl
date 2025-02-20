using CDDLib

"""
    struct CSCB{T}
   
Represents a centrally symmetric convex body.
Only half of its vertices should be stored in `positive_vertices`,
the other ones are given by `-positive_vertices`.
However, we store all normals.
Normals are normalized such that
the body is given by ``\\{x | h^T x \\leq 1, h \\in normals\\}``.
"""
struct CSCB{T}
    positive_vertices :: Vector{Vector{T}}
    normals ::  Vector{Vector{T}}
end

"""
    calculate_normals(vertices[, symmetric])

Given the vertices of a polytope, compute
the normals ``w`` such that ``w^T v <= 1``.

If `symmetric` is true, it is enough to
include the vertices in one halfspace.
"""
function calculate_normals(vertices; symmetric = true)
    if symmetric
        vertices = [vertices; map(x->-x, vertices)]
    end
    P = polyhedron(vrep(vertices), CDDLib.Library(:exact))
    return [h.a/h.β for h in halfspaces(hrep(P))]
end

"""
    all_vertices(cscb)

Calculate all vertices of the centrally symmetric
convex body cscb.
"""
function all_vertices(cscb :: CSCB)
    return [cscb.positive_vertices; map(x->-x,cscb.positive_vertices)]
end

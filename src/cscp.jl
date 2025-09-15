using CDDLib

"""
    struct CSCP{T}
   
Represents a centrally symmetric convex polytope.
Only half of its vertices should be stored in `positive_vertices`,
the other ones are given by `-positive_vertices`.
We store, however, all normals.
Normals are normalized such that
the body is given by ``\\{x | h^T x \\leq 1, h \\in normals\\}``.
"""
struct CSCP{T}
    positive_vertices::Vector{Vector{T}}
    normals::Vector{Vector{T}}
end

function CSCP(V)
    normals = calculate_normals(V)
    T = typeof(normals)
    return CSCP(T(V), normals)
end

function Base.:(*)(P::AbstractMatrix, cscp::CSCP)
    CSCP([P*v for v in cscp.positive_vertices])
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
        vertices = [vertices; map(x -> -x, vertices)]
    end
    P = polyhedron(vrep(vertices), CDDLib.Library(:exact))
    return [h.a / h.β for h in halfspaces(hrep(P))]
end

"""
    all_vertices(cscp)

Calculate all vertices of the centrally symmetric
convex body cscp.
"""
function all_vertices(cscp::CSCP)
    return [cscp.positive_vertices; map(x -> -x, cscp.positive_vertices)]
end

cube = CSCP(
    [[1, -1, -1], [1, 1, -1], [1, -1, 1], [1, 1, 1]],
    [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]]
)

using Polyhedra
using CDDLib
using DynamicPolynomials
using LinearAlgebra
import SetProg
import LazySets
using SCIP
using JuMP
using Clarabel

"""
    find_shadiness_constant_numerically(cscp, silent)

Uses SCIP to find a numerical approximation
to the shadiness constant of cscp.
If `silent` is true, no information
is printed during the optimization.
"""
function find_shadiness_constant_numerically(cscp, silent=true)
    model = Model(SCIP.Optimizer)
    if silent
        set_silent(model)
    end
    vertices = convert.(Vector{Float64}, cscp.positive_vertices)
    normals = convert.(Vector{Float64}, cscp.normals)
    @variable(model, P[1:3, 1:3])
    @variable(model, β)
    inequalities = [β - w' * P * v for v in vertices for w in normals]
    # β >= w'*p*v
    inequalities = [inequalities; β - 1]
    # We already know that the norm of a projection is at least one
    eq1 = [tr(P) - 2]
    eq2 = reshape(P * P - P, 9)
    equations = [eq1; eq2]
    for c in equations
        @constraint(model, c == 0)
    end
    for c in inequalities
        @constraint(model, c >= 0)
    end
    projection = I
    norm = -1
    try
        @objective(model, Min, β)
        optimize!(model)
        projection = value.(P)
        norm = value.(β)
    catch e
        @error "Problem is numerically problematic."
    end
    return (projection, norm)
end

"""
   determine_matrix_ball_vertices(cscp; rational = true)

Consider the norm with unit ball `cscp`.
This function determines the vertices
of the unit ball of the corresponding operator norm.

If `rational` is true, exact rational arithmetic is used,
otherwise the computation is done in floating point.
"""
function determine_matrix_ball_vertices(cscp::CSCP; rational=true)
    vertices = all_vertices(cscp)
    normals = cscp.normals
    d = length(vertices[1])
    @polyvar p[1:d, 1:d]
    ineq = [w' * p * v for v in vertices for w in normals]
    half_spaces = [Polyhedra.HalfSpace(coefficients(ie, reshape(p, 9)), 1) for ie in ineq]

    P = (
        if rational
            polyhedron(hrep(half_spaces), CDDLib.Library(:exact))
        else
            polyhedron(hrep(half_spaces), CDDLib.Library(:float))
        end
    )

    return [reshape(A,3,3) for A in points(P)]
end

"""
    determine_squared_spectralnorm_bound(cscp)

Consider the norm ``||.||_C`` whose unit ball is `cscp`.
Then
``r ||.||_C \\leq ||.||_2 \\leq R ||.||_C``.
This implies
``||A||_2 = \\sup_x ||Ax||_2 / ||x||_2 \\leq R/r ||A||_C``.
as well as 
``||A||_C = \\sup_x ||Ax||_C / ||x||_C \\leq R/r ||A||_2``.

This method deterimens R^2/r^2.
The ratio is returned squared in order to return exact results for rational inputs.
"""
function determine_squared_spectralnorm_bound(cscp::CSCP)
    vertices = cscp.positive_vertices
    normals = cscp.normals
    # squared radius of circumscribed centered sphere
    R_squared = maximum(sum(a^2 for a in x) for x in vertices)
    # Let b be a point where the inscribed sphere touches the cscp. Then
    # we have for the corresponding normal w
    # w^T*b = 1, b = a w
    # => a ||w||_2^2 = 1
    # => ||b||_2 = a ||w||_2 = 1/||w||_2
    # => ||b||_2^2 = 1/||w||_2^2

    # squared radius of inscribed centered circle
    r_squared = minimum(Rational{BigInt}(1) / sum(t^2 for t in w) for w in normals)
    return R_squared / r_squared
end

"""
    approximate_john_position(cscp; prec=10^1)

Given a CSCP M=``cscp``
this functions calculates a rational approximation (up to precision `prec`) of a matrix ``Q``,
such that ``QM`` is approximately in John position, i.e.
the maximal volume inscribed ellipsoid is the Euclidean unit ball.
"""
function approximate_john_position(cscp; prec=10^1)
    P = polyhedron(vrep(all_vertices(cscp)), CDDLib.Library(:float))
    lazypolytope = LazySets.HPolytope(P)
    john_ellipsoid = LazySets.underapproximate(
        lazypolytope, LazySets.Ellipsoid; backend=Clarabel.Optimizer,interior_point=[0.0; 0.0; 0.0]
    )
    L = eigvecs(LazySets.shape_matrix(john_ellipsoid))
    R = diagm(map(sqrt, eigvals(LazySets.shape_matrix(john_ellipsoid))))
    Q = map(x -> fast_round(x, prec), inv(L * R))
    return Q
end

"""
    approximate_projection(A;atol, prec)

Given a floating point approximation of
a projection `A`, calculate a
rational approximation (using precision `prec`)
to `A` which is a exact projection.

Only works in dimension 3.

`atol` is used as a bound for the minimal
singular value to determine the kernel of `A`.
"""
function approximate_projection(A; atol=1e-5, prec=10^4)
    d = size(A, 1)
    if d != 3
        error("Only 3x3 matrices A are supported by approximate_projection.")
    end
    T = typeof(A[1, 1])
    v = nullspace(A; atol=atol)[:, 1]

    # If A e_j = 0, then v is a multiple of e_j and
    # v_j would be in absoulte value the maximal entry of v 
    j = argmax(map(abs, v))
    image_vectors = [A[:, i] for i in 1:d if i != j]

    w = cross(image_vectors[1], image_vectors[2])
    v = fast_round.(v, prec)
    w = fast_round.(w, prec)
    return I - v * w' / (w' * v)
end

"""
    operator_norm(cscp, A)

Determines the operator norm of the linear map `A`
relative to the norm whose unit ball
is `cscp`.
"""
function operator_norm(cscp::CSCP, A)
    vertices = cscp.positive_vertices
    return maximum([h' * A * v for h in cscp.normals for v in vertices])
end

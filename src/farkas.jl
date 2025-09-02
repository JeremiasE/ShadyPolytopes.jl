using Polyhedra
using CDDLib
using HiGHS
using CSV
using DataFrames

"""
    find_single_farkas_certificate(cscb, w, α)

Determine a rational certificate showing that there is no
projection with
image ``w^⟂``, whose norm relative to cscb is
larger then α
The certificate consists of a non-negative rational vector ``r`` such
that ``A^⊤ r = w`` for an appropriately defined matrix ``A``.
We can always get vectors where only d entries
are non-zero, so we return these ``d`` indices
and the ``d`` corresponding entries.
"""
function find_single_farkas_certificate(cscb::CSCB, w, α)
    A = matrix_from_projection_normal(cscb, w, α)
    d = size(A, 1)
    model = Model(HiGHS.Optimizer)
    set_silent(model)
    @variable(model, y[1:d] >= 0)
    @constraint(model, A' * y == w)
    optimize!(model)
    z = value.(y)
    if !(is_solved_and_feasible(model))
        error("Problem is infeasible for w=", w, ", α=", α)
    end
    # Caratheodory's theorem ensures the existence
    # of a non-negative solution to A^⊤ y = w
    # where only n entries of y are non-zero.
    # It turns out that in our application HiHGS actually finds these solutions.
    # We only get an approximate solution from HiGHS, but
    # knowing the indices where the solution is non-zero,
    # we can try to get a rational solution using the resulting d×d system of linear equations
    bv = (z .!= 0)
    r = A'[:, bv] \ w
    return (1:d)[bv], r
end

"""
    matrix_from_projection_normal(cscb, w, α)

Calculate the matrix A such that
The existence of a positive solution of ``A^T y = w``
certififies the non-existence of a projection
described in `find_single_farkas_certificate`
"""
function matrix_from_projection_normal(cscb::CSCB, w, α)
    V = cscb.positive_vertices
    H = cscb.normals
    return Matrix{Rational{BigInt}}(
        vcat([h' * x * w' - w' * x * h' - α * w' for x in V for h in H]...)
    )
end

"""
    find_nonnegative_solution(A, w)

Finds a non-negative solution
to ``A^⊤ y = w`` using rational arithmetic CDDLib.
This sadly is very slow and therefore not used in the methods above.
"""
function find_nonnegative_solution(A, w)
    # We want A^⊤ y = w, y ≥ 0
    n = size(A, 1)
    eq = [HyperPlane(Vector{Rational{Int}}(r), w[i]) for (i, r) in enumerate(eachcol(A))]
    # HalfSpace([1, 1], 1) == x+y ≤ 1

    ineq = [
        HalfSpace(Vector{Rational{Int}}([i == j ? -1//1 : 0 for i in 1:n]), 0) for j in 1:n
    ]
    P = polyhedron(hrep(ineq), CDDLib.Library(:exact))
    vrep(P)
    return points(P)
end

"""
    add_one_at_position_k(v,k)

Adds a ``1`` into v at position k. 
"""
function add_one_at_position_k(v, k)
    w = copy(v)
    insert!(w, k, 1)
    return w
end

"""
    generate_farkas_certificates(cscb[, k, α; io])

Generate rational certificates as in `find_single_farkas_certificate`
by sampling the surface of the unit-cube. Output to `io`.
"""
function generate_farkas_certificates(cscb, k=1000, α=84//83, pos=1; io=stdout)
    for i in (-k):k
        for j in (-k):k
            w = add_one_at_position_k([i//k, j//k], pos)
            indices, r = find_single_farkas_certificate(cscb, w, α)
            println(io, join(w, ";"), ";", join(indices, ";"), ";", join(r, ";"))
            flush(io)
        end
    end
end

"""
     check_single_farkas_certificate(cscb, w, y, α)

Verify a single certificate generated
by `find_single_farkas_cerfificate`, the output of this
method is translated into `y`.
"""
function check_single_farkas_certificate(cscb, w, y, α)
    if !all(y .>= 0)
        error(y, " is not non-negative.")
    end
    A = matrix_from_projection_normal(cscb, w, α)
    if !(A' * y == w)
        error(y, " is not a solution.")
    end
end

"""
    check_farkas_certificate(cscb, file, α)

Verify the output of `generate_farkas_certificate`.
"""

function check_farkas_cerfificate_file(cscb, file, α=84//83; silent=true)
    df = CSV.read(file, DataFrame; header=false)
    V = cscb.positive_vertices
    H = cscb.normals
    n = length(V) * length(H)
    for (i, row) in enumerate(eachrow(df))
        if !silent
            println(i)
        end
        w = parse.(Rational{BigInt}, collect(row[1:3]))
        indices = collect(row[4:6])
        yentries = parse.(Rational{BigInt}, collect(row[7:9]))
        y = zeros(Rational{BigInt}, n)
        y[indices] = yentries
        check_single_farkas_certificate(cscb, w, y, α)
    end
end

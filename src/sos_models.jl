"""
    shadiness_via_projection_matrix(cscb, add_beta_bound, feasibility)

Produces a semialgebraic set which allows to determine
the shadiness constant of `cscb` via straightforward
representation of the projection as a matrix.

If `add_beta_bound` is true, we add the (redundant)
condition that the norm of a projection is at least 1.

If `feasibility` is true, the produced
set is used to show that the shadiness
constant is at least `bound`.

Otherwise the produced set can be used to calculate a lower bound for
the shadiness constant.
"""
function shadiness_via_projection_matrix(cscb;
                                         add_beta_bound=true, 
                                         feasibility=true, bound=101//100)
    @polyvar p[1:3,1:3]
    @polyvar β
    vertices = cscb.positive_vertices
    normals = cscb.normals
    
    if feasibility
        vars = reshape(p,9)
    else
        vars = [reshape(p,9);β]
    end

    ineq = [(1+β)-w'*p*v for v in vertices for w in normals]
    eq1 = [tr(p)-2]
    eq2 = reshape(p*p-p,9)
    eq = [eq1;eq2]

    if add_beta_bound
        ineq = [β;ineq]
    end

    L = algebraic_set(eq)
    K = basicsemialgebraicset(L,ineq)
    return L, K, vars
end


"""
    shadiness_via_vw(cscb;feasibility, v_index, w_index, v_entry=1,
                          bound, quadratic)

Produces a semialgebraic set which allows to determine
the shadiness constant of `cscb` via the
representation of a rank d-1 projection by a normal to the image
and a vector spanning the kernel.

If `quadratic` is true, we represent the bounds
on the entries of the vector by a quadratic ineqality,
otherwise we use two linear ones.

If `feasibility` is true, the produced
set is used to show that the shadiness
constant is at least `bound`.

Otherwise the produced set can be used to calculate a lower bound for
the shadiness constant.
"""
function shadiness_via_vw(cscb;
                          feasibility = true,
                          v_index=3, w_index=3, v_entry=1,
                          bound=101//100, quadratic=true)
    @polyvar w[1:2]
    @polyvar v[1:2]
    @polyvar β
    vertices = cscb.positive_vertices
    normals = cscb.normals

    if feasibility
        vars = [w;v]
    else
        vars = [w;v;β]
    end

    vv, ww = v_w_vectors(v_index,w_index, v_entry, v, w)

    if feasibility
        ineq1=[(1+β)*ww'*vv+ww'*x*h'*vv-h'*x*ww'*vv for x in vertices for h in normals]
    else
        ineq1=[bound*ww'*vv+ww'*x*h'*vv-h'*x*ww'*vv for x in vertices for h in normals]
    end

    ineq2 = [ww'*vv]

    if quadratic
        ineq3 = [1-k^2 for k in [v;w]]
    else
        ineq3 = [[1-k for k in [v;w]]; [k+1 for k in [v;w]]]
    end

    ineq = [β; ineq3; ineq2; ineq1]

    L = FullSpace()
    K = basic_semialgebraic_set(L,ineq)
    return L, K, vars
end


"""
    v_w_vectors(v_index, w_index, v_entry, v_variables, w_variables)

Inserts `v_entry` into the vector `v_variables` at position `v_index`
and inserts 1 into the vector `w_variables` at position `w_index`.
"""
function v_w_vectors(v_index, w_index, v_entry, v_variables, w_variables)
    v = Vector{typeof(1+v_variables[1])}([v_variables...])
    w = Vector{typeof(1+w_variables[1])}([w_variables...])
    insert!(v,v_index,v_entry)
    insert!(w,w_index,1)
    return v,w
end

"""
    putinar_model(solver, vars, K, obj; feasibility, maxdegree)

Returns a JuMP Model which can be used to show that the function `obj`
using the variables `vars` is non-negative over the basic
semialgebraic set `K` if `feasibility` is true and otherwise find a
proveable lower bound for `obj` over `K`.  This is done using
Putinar's Positivstellensatz using polynomials and
a SOS decomposition using polynomials of maximal degree `maxdegree`.
"""

function putinar_model(solver,vars,K,obj;feasibility=true, maxdegree=4)
    model = SOSModel(solver)
    if feasibility
        @constraint(model, c, obj >= 0, domain = K, maxdegree = maxdegree)
    else
        @variable(model, τ)
        @objective(model, Max, τ)
        @constraint(model, c, obj >= τ, domain = K, maxdegree = maxdegree)
    end
    return model
end

"""
    solve_via_vw(cscb, solver; feasibility=true, bound=101//100, maxdegree=10)

If `feasibility` is true, calculate a SOS decomposition
that certifies that the shadiness constant of `cscb` is at least `bound`.
Otherwise determine a lower bound for the shadiness constant.

Uses `shadiness_via_vw`, more on the other parameters there.
"""


function solve_via_vw(cscb, solver; feasibility=true, bound=101//100, maxdegree=10)
    L, K, vars = shadiness_via_vw(cscb; feasibility = feasibility, bound=bound)
    if feasibility
        obj = -1+0*vars[1]
    else
        obj = vars[end]
    end
    model = putinar_model(solver, vars, K, obj; feasibility=feasibility, maxdegree=maxdegree)
    println(model)
    optimize!(model)
    println(solution_summary(model))
    return L, K, vars, model
end

"""
    solve_via_projection_matrix(cscb, solver; feasibility=true, bound=101//100, maxdegree=10)

If `feasibility` is true, calculate a SOS decomposition
that certifies that the shadiness constant of `cscb` is at least `bound`.
Otherwise determine a lower bound for the shadiness constant.

Uses `shadiness_via_projection_matrix`, more on the other parameters there.
"""
function solve_via_projection_matrix(cscb, solver; feasibility=true, bound=101//100, maxdegree=10)
    L, K, vars = shadiness_via_projection_matrix(cscb; feasibility = feasibility, bound=bound)
    if feasibility
        obj = -1+0*vars[1]
    else
        obj = vars[end]
    end
    model = putinar_model(solver, vars, K, obj; feasibility=feasibility, maxdegree=maxdegree)
    println(model)
    optimize!(model)
    println(solution_summary(model))
    return L, K, vars, model
end


"""
    check_putinar_certificate(model,K)

Calculates the weighted sos decomposition calculated in `model`, this
should approximately equal the objective function used in the model
minus a lower bound.
    """
function check_putinar_certificate(model,K)
    return polynomial(sos_decomposition(model[:c],K))
end







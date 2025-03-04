using DynamicPolynomials
using ShadyPolytopes
using SumOfSquares
using SCS
using Clarabel

function small_projection_sos_set(feasibility = true; bound=93//100)
    vertices = [[1,1,0],[1,-1,0],[0,0,1]]
    normals = calculate_normals(vertices, symmetric = true)
    cscb = CSCB{Rational{BigInt}}(vertices,normals)
    return shadiness_via_projection_matrix(cscb;
                                           use_beta_bound=false,
                                           use_norm_bound=true,
                                           use_equations=true,
                                           feasibility=feasibility,
                                           bound=bound)
end

function tiny_sos_set(feasibility = true, equation = false; bound=174//100)
    @polyvar x y     
    vars = [x,y]
    # We are looking to
    # min(y) over our set
    # To show that min(y)>bound, we add y ≤ bound
    # and show that this is impossible
    
    if feasibility
        # 4 ≥ x^2+y^2, y ≤ bound
        ineq = [4-x^2-y^2,4-(x-2)^2-y^2,bound-y]
    else
        # 4 ≥ x^2+y^2, 4 ≥ (x-2)^2+y^2
        ineq = [4-x^2-y^2,4-(x-2)^2-y^2]
    end
    if equation
        L = algebraic_set([x-y])
    else
        L = FullSpace()
    end
    K = basic_semialgebraic_set(L,ineq)
    return L, K, vars
end

function optimize_small_projection_sos_model(feasibility = true; bound=95//100, maxdegree=3, solver = Clarabel.Optimizer)
    L, K, vars, squared_variable_bound = small_projection_sos_set(feasibility, bound=bound)
    if feasibility
        obj = -1+0*vars[1]
    else
        obj = vars[end]
    end
    model = putinar_model(solver,vars,K,obj;feasibility=feasibility, maxdegree=maxdegree)
    optimize!(model)
    return L, K, vars, squared_variable_bound, obj, model
end

function optimize_tiny_sos_model(feasibility = true, equation = false; bound=174//100, maxdegree=2, solver=SCS.Optimizer)
    L, K, vars = tiny_sos_set(feasibility, equation, bound=bound)
    if feasibility
        obj = -1+0*vars[1]
    else
        obj = vars[end]
    end
    model = putinar_model(solver,vars,K,obj;feasibility=feasibility, maxdegree=maxdegree)
    optimize!(model)
    return L, K, vars, obj, model
end


function certificate_embryo_feasibility()
    L, K, vars, model = test_embryo_feasibility_model()
    sos, offset_ineq_sos, given_ineq_sos, given_ineq = round_sos_decomposition(model, K, -1, vars)
    print_certificate(sos, offset_ineq_sos, given_ineq_sos, given_ineq, vars)
end

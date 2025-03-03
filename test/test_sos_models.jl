using DynamicPolynomials
using ShadyPolytopes
using SumOfSquares
using SCS

function tiny_sos_set(feasibility = true, equation = false; bound=174//100)
    @polyvar x y     
    vars = [x,y]
    if feasibility
        # 4 ≥ x^2+y^2, y ≥ bound
        ineq = [4-x^2-y^2,4-(x-2)^2-y^2,y-bound]
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

function small_sos_set(feasibility = true)
    @polyvar x y
    vars = [x,y]
    ineq = [4-x^2-y^2]
    eq = [4-(x-2)^2-y^2]
    L = algebraic_set(eq)
    K = basic_semialgebraic_set(L,ineq)
    return L, K, vars
end    

function optimize_tiny_sos_model(feasibility = true, equation = false; bound=174//100)
    L, K, vars = tiny_sos_set(feasibility, equation, bound=bound)
    if feasibility
        obj = -1+0*vars[1]
    else
        obj = vars[end]
    end
    model = putinar_model(SCS.Optimizer,vars,K,obj;feasibility=feasibility, maxdegree=2)
    optimize!(model)
    return L, K, vars, model
end


function certificate_embryo_feasibility()
    L, K, vars, model = test_embryo_feasibility_model()
    sos, offset_ineq_sos, given_ineq_sos, given_ineq = round_sos_decomposition(model, K, -1, vars)
    print_certificate(sos, offset_ineq_sos, given_ineq_sos, given_ineq, vars)
end

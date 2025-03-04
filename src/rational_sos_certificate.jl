import MultivariatePolynomials
import SemialgebraicSets

"""
    gram_to_sos(RG, monos)

Takes a rational gram matrix and turns it into a rational SOS
decomposition.
"""
function gram_to_sos(RG, monos)
    BK = bunchkaufman(RG)
    d = diag(BK.D)
    L = BK.U'*BK.P
    return RationalSOSDecomposition(L*monos, d)
end

"""
    round_gram_matrix(gram, obj, degree,vars;prec)

Takes an approximate gram matrix and produces a rational gram_matrix
which represents `obj` as a rational SOS decomposition using monomials
in the variables `vars` up to degree `degree`. The rounding to
rational values is done up to precision `prec`.
"""
function round_gram_matrix(gram, obj, degree,vars; prec=10^7)
    G = value_matrix(gram)
    RG = fast_round.(G,prec)
    n = size(RG,1)
    monos = collect(gram.basis)
    new_monos = [1;monomials(vars,1:degree)]
    rounded_obj = monos'*RG*monos
    error = rounded_obj - obj
    m = length(new_monos)
    projected_G = zeros(Rational{BigInt},m,m)
    mono_matrix = new_monos*new_monos'
    for i in 1:m
        for j in 1:m
            current_mono = new_monos[i]*new_monos[j]
            coeff = coefficients(error,[current_mono])[1]
            nc = count(mono_matrix .== current_mono)
            ipos = (monos .== new_monos[i])
            jpos = (monos .== new_monos[j])
            if (count(ipos)==1 && count(jpos) ==1)
                projected_G[i,j] = RG[ipos,jpos][1,1] - coeff // nc
            else
                projected_G[i,j] = -coeff // nc
            end
        end
    end
    return projected_G, new_monos
end
 



"""
    reduce_one_minus_monomial(t, vars)

Returns ``r_1,...,r_k``
such that ``1-t = \\sum_i r_i*(1-vars_i)``.
"""
function reduce_one_minus_monomial(mono, vars)
    if degree(mono) == 0
        return zeros(typeof(vars[1]+big(1)//1), length(vars))
    end
    
    contained_vars = effective_variables(mono)
    s = contained_vars[1]
    if degree(mono) == 1
        qv = zeros(typeof(vars[1]+big(1)//1), length(vars))
        qv[vars .== s] .= 1
        return qv
    end
    # 1-xy = (1-x) + x(1-y)
    rv = reduce_one_minus_monomial(monomials(mono ÷ s)[1], vars)
    qv = rv * s
    qv[vars .== s] = qv[vars.==s] .+ 1
    return qv
end

"""
    reduce_n_minus_monomial(t, vars, n)

Returns ``N``, ``r_1,…,r_k``
such that ``N-t = \\sum_i r_i*(n-vars_i)``.

If ``(1-x^α) = \\sum_i (1-x_i)*p_i(x_1,…,x_k)``,
then ``(n^{|α|}-x^α = \\sum_i (n-x_i) n^α p_i(x_1/n,…,x_k/n)``.
"""
function reduce_n_minus_monomial(mono, vars, n)
    α = degree(mono)
    if α == 0
        return 1, zeros(typeof(vars[1]+big(1)//1), length(vars))
    end
    N = n^α
    qv = reduce_one_minus_monomial(mono, vars)
    rv = [n^(α-1)*q(vars => vars*1//n) for q in qv]
    return N, rv
end

"""
    offset_sum_of_monomials(mono_vector,vars, n)

If `mono_vector` = (m_1,...,m_k) and `vars`=x_1,...,x_l
and the result is N and p_1,...,p_l,
then ``N - m_1^2+m_2^2+...+m_k^2`` as ``p_1 (factor-x_1^2)+...+ p_l (factor-x_l^2)``.

"""
function offset_sum_of_monomials(mono_vector,vars; n=1//1)
    sumN = 0
    sumpv = zeros(typeof(vars[1]+big(1)//1), length(vars))

    for m in mono_vector
        N, pv = reduce_n_minus_monomial(m, vars,n)
        sumN  += N
        sumpv += pv
    end   
    #println(typeof(p))
    #println(typeof(polynomial.(Rational{BigInt}(1)*collect(monomials(p[1])))))
    #println(typeof(coefficients(p[1])))
     
    return sumN, [RationalSOSDecomposition(Rational{BigInt}(1)*collect(monomials(monosum)),
                                     coefficients(monosum))
            for monosum in sumpv]
end

"""
    rationally_reduce_sos_decomp(model, K, obj, vars; feasibility, prec)

Calculate rational approxmiations of all parts of the computed
weighted SOS decomposition of `obj` besides the main SOS part. 
"""
function rationally_reduce_sos_decomp(model, K, obj, vars; feasibility=true, prec = 10^7)
    ineqs = SemialgebraicSets.inequalities(K)
    ineqs_indices = eachindex(ineqs)
    eqs = SemialgebraicSets.equalities(K)
    eqs_indices = eachindex(eqs)
    
    if feasibility
        new_obj = obj
    else
        # obj = τ + sos + sum(ineq_part)+sum(eq_part)
        obj = new_obj - model[:τ]
    end
    
    lang_mul = lagrangian_multipliers(model[:c])
    
    ineqs_part = sum(SOSDecomposition(lang_mul[i])*ineqs[i] for i in ineqs_indices)
    sos_part = sos_decomposition(model[:c])
    eqs_part = new_obj-sos_part-ineqs_part

    new_L = algebraic_set(equalities(K))
    I = ideal(new_L)
    SemialgebraicSets.compute_gröbner_basis!(I)
    new_eqs = equalities(new_L)
    eqs_poly_coeffs = div(eqs_part,I.p)
    new_eqs_indices = eachindex(new_eqs)
    
    #here we have to add the transformation from the gröbner basis to the original generators
    A,B = lift_polynomials(eqs, new_eqs, vars)
    # A'*eqs = new_eqs
    # new_eqs' * new_eqs_poly_coeffs = eqs'*A*new_eqs_poly_coeffs
    
    rounded_new_eqs_poly_coeffs = [round_poly(p,prec) for p in eqs_poly_coeffs]
    rounded_eqs_poly_coeffs = [sum(A[k,i]*rounded_new_eqs_poly_coeffs[i]
                                   for i in eachindex(rounded_new_eqs_poly_coeffs))
                               for k in eachindex(eqs)]
    
    # BUG: This does not work - A*rounded_new_eqs_poly_coeffs                     

    rounded_eqs_part = sum(rounded_eqs_poly_coeffs[i]*eqs[i] for i in eqs_indices)
    
    rounded_ineqs_sos_decomp = [round_sos(SOSDecomposition(lang_mul[i]), prec) for i in ineqs_indices]
    rounded_ineqs_part = sum(rounded_ineqs_sos_decomp[i]*ineqs[i] for i in ineqs_indices)

    return (new_obj-rounded_ineqs_part-rounded_eqs_part,
            rounded_ineqs_sos_decomp,
            rounded_eqs_poly_coeffs,
            rounded_ineqs_part,
            rounded_eqs_part,
            eqs_part,
            new_L
            )
end

"""
    round_sos_decomposition(model, K, obj, vars, squared_variable_bound, offset; prec)

Turn the calculated approximate weighted SOS decomposition into
an exact rational weighted SOS decomposition.
"""
function round_sos_decomposition(model, K, obj, vars, squared_variable_bound, offset=1//10^4; prec=10^2)
    r = rationally_reduce_sos_decomp(model, K, obj, vars; prec=prec, feasibility=true)
    remobj, ineqs_sos, eqs_poly_coeffs, ineqs_part, eqs_part = r
    
    remobj_degree = maximum(map(degree, terms(remobj)))
    approx_gram = gram_matrix(model[:c])
    gram_degree = maximum(map(degree, collect(approx_gram.basis)))
    comb_degree = max(ceil(Int64, remobj_degree / 2), gram_degree)   
    RG, new_monos = round_gram_matrix(approx_gram, remobj, comb_degree,vars; prec=prec)

    # nm' RG nm = remobj
    rounded_sos = gram_to_sos(RG+offset*I, new_monos)
    rounded_sos_part = polynomial(rounded_sos)
    
    #return new_monos,soss,ineqs,offset,vars
    N, offset_sos = offset_sum_of_monomials(new_monos[2:end],vars; n=squared_variable_bound)
    for i in eachindex(vars)
        offset_sos[i].wv *= offset
    end
    #k = length(new_monos)
    
    #offsetpart = k-newmonos'*newmonos
    offset_part = (sum((squared_variable_bound-t^2)*polynomial(p) for (p,t) in zip(offset_sos,vars)))
    left_hand_side = rounded_sos_part + offset_part + ineqs_part + eqs_part
    println(left_hand_side)
    cert = RationalPutniarCertificate(vars,
                                      squared_variable_bound,
                                      equalities(K),
                                      inequalities(K),
                                      rounded_sos,
                                      offset_sos,
                                      eqs_poly_coeffs,
                                      ineqs_sos,
                                      coefficients(left_hand_side)[1])
    return cert
    # newobj = -1 - ineq_part = sos_part + offset_part 
    # return  sos_part + offset_part + ineq_part
end

"""
    RationalPutniarCertificate

Collects all the information to produce a rational putinar certificate
that shows that a set is empty. 
"""

mutable struct RationalPutniarCertificate
    vars # the variable
    squared_variable_bound
    eqs # the equalities
    ineqs # the inequalities
    sos # the SOS part
    offset_sos # weighted sos for the variable bounds
    eqs_polys
    ineqs_sos # SOS
    left_hand_side # a negative rational
end

"""
    check_rational_putinar_certificate(rpc)

Check that the certificate is correct.
"""

function check_rational_putinar_certificate(rpc)
    if rpc.left_hand_side < 0
        println("Right hand side negative 😄")
    else
        println("Right hand side not negativ 🥲")
    end

    sos_part = polynomial(rpc.sos)
    offset_part = sum(polynomial(p)*(rpc.squared_variable_bound-t^2)
                      for (t,p) in zip(rpc.vars,rpc.offset_sos))
    eqs_part = sum(p*q
                   for (q,p) in zip(rpc.eqs, rpc.eqs_polys))
    ineqs_part = sum(polynomial(p)*q
                     for (q,p) in zip(rpc.ineqs, rpc.ineqs_sos))
    
    right_hand_side = sos_part + offset_part + eqs_part + ineqs_part
    if (right_hand_side- rpc.left_hand_side)==0
        println("Left hand side = right hand side 😄")
    else
        println("Left hand side ≠ right hand side 🥲")
    end
end

function print_certificate(sos,  offset_ineq_sos, given_ineq_sos, given_ineq, vars; io=stdout)
    println(io,"(")
    println(io,sos)
    for i in eachindex(vars)
        println(io,"+(",(1-vars[i]^2),")*(",offset_ineq_sos[i],")")
    end
    for i in eachindex(given_ineq)
        println(io,"+(",given_ineq[i],")*(",given_ineq_sos[i],")")
    end
    println(io,")")
end

function print_certificate_julia(sos, offset_ineq_sos, given_ineq_sos, given_ineq, vars; io=stdout)
    println(io,"using DynamicPolynomials")
    println(io,"@polyvar v[1:2] w[1:2]")
    println(io,"const PolyType = Polynomial{DynamicPolynomials.Commutative{DynamicPolynomials.CreationOrder}, Graded{LexOrder}, Rational{BigInt}}")
    print(io,"sos_result=")
    show_big(io, sos)
    println(io,"")
    
    println(io,"offset_ineq = Vector{PolyType}()")
    println(io,"offset_ineq_sos =  Vector{PolyType}()")   
    
    for i in eachindex(vars)
        println(io, "println(",i,")")
        println(io,"push!(offset_ineq,",(1-vars[i]^2),")")
        print(io,"push!(offset_ineq_sos,")
        show_big(io,offset_ineq_sos[i])
        println(io,")")
    end
    println(io,"offset_ineq_result=sum(offset_ineq[i]*offset_ineq_sos[i] for i in eachindex(offset_ineq))")

    println(io,"given_ineq = Vector{PolyType}()")
    println(io,"given_ineq_sos =  Vector{PolyType}()")
   
    for i in eachindex(given_ineq)
        println(io, "println(",i,")")
        print(io, "push!(given_ineq,")
        show_big(io, given_ineq[i])
        println(io, ")")
        print(io,"push!(given_ineq_sos,")
        show_big(io, given_ineq_sos[i])
        println(io, ")")
    end
    println(io,"given_ineq_result=sum(given_ineq[i]*given_ineq_sos[i] for i in eachindex(given_ineq))")
    println(io,"println(sos_result+offset_ineq_result+given_ineq_result)")
end


function show_big(io::IO, a :: Integer)
    print(io, "big(",a,")")
end

function show_big(io::IO, a :: Rational)
    show_big(io, numerator(a))
    print(io, "//")
    show_big(io, denominator(a))
end

function show_big(io::IO, a :: MultivariatePolynomials.Polynomial)
    monos = monomials(a)
    coeffs = coefficients(a)
    n = length(monos)
    for i in 1:n
        show_big(io,coeffs[i])
        print(io,"*",monos[i])
        if i<n
            print(io,"+")
        end
    end
end

function show_big(io::IO, d :: RationalSOSDecomposition)
    for (i, p, w) in zip(eachindex(d.pv),d.pv,d.wv)
        show_big(io, w)
        print(io, "*")
        print(io, "(")
        show_big(io, polynomial(p))
        print(io, ")^2")
        if i != length(d.pv)
            print(io, " + ")
        end
    end
end

function show_big(io::IO, d :: SOSDecomposition)
    for (i, p) in enumerate(d.ps)
        print(io, "(")
        show_big(io, p)
        print(io, ")^2")
        if i != length(d.ps)
            print(io, " + ")
        end
    end
end

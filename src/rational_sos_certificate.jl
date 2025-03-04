import MultivariatePolynomials
import SemialgebraicSets
import DynamicPolynomials

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

Assumes that the first |`vars`| inequalities in K are of the form
(squared_variable_bound-x_i^2). 
"""
function round_sos_decomposition(model, K, obj, vars, squared_variable_bound,
                                 offset=big(1)//10^4;
                                 prec=big(10^2))
    r = rationally_reduce_sos_decomp(model, K, obj, vars; prec=prec, feasibility=true)
    remobj, ineqs_sos, eqs_poly_coeffs, ineqs_part, eqs_part = r
    println("Finished rounding equality and inequalitiy weights.")
    
    remobj_degree = maximum(map(degree, terms(remobj)))
    approx_gram = gram_matrix(model[:c])
    gram_degree = maximum(map(degree, collect(approx_gram.basis)))
    comb_degree = max(ceil(Int64, remobj_degree / 2), gram_degree)   
    RG, new_monos = round_gram_matrix(approx_gram, remobj, comb_degree,vars; prec=prec)
    println("Finished projecting Gram matrix.")
    
    # nm' RG nm = remobj = rounded_sos - offset *new_monos'*new_monos
    rounded_sos = gram_to_sos(RG+offset*I, new_monos)
    rounded_sos_part = polynomial(rounded_sos)
    if remobj ≠ rounded_sos_part- offset *new_monos'*new_monos
        println(remobj-rounded_sos_part)
        error("Offset is to small, rounded Gram matrix + offset is not pd.")
    end
    println("Finished converting Gram matrix to SOS.")
    
    #return new_monos,soss,ineqs,offset,vars
    N, offset_sos = offset_sum_of_monomials(new_monos[2:end],vars; n=squared_variable_bound)
    for i in eachindex(vars)
        offset_sos[i].wv *= offset
    end
    #k = length(new_monos)
    
    #offsetpart = k-newmonos'*newmonos
    offset_part = (sum((squared_variable_bound-t^2)*polynomial(p) for (p,t) in zip(offset_sos,vars)))
    println("Finished calculating offset weights.")
    left_hand_side = rounded_sos_part + offset_part + ineqs_part + eqs_part

    for i in eachindex(vars)
        append!(offset_sos[i].pv, ineqs_sos[i].ps)
        append!(offset_sos[i].wv, ones(Rational{BigInt},length(ineqs_sos[i].ps)))
    end
    
    ineqs = inequalities(K)
    remaining_ineqs = ineqs[length(vars)+1:end]
    remaining_ineqs_sos = ineqs_sos[length(vars)+1:end]
    
    cert = RationalPutniarCertificate(vars,
                                      squared_variable_bound,
                                      equalities(K),
                                      remaining_ineqs,
                                      rounded_sos,
                                      offset_sos,
                                      eqs_poly_coeffs,
                                      remaining_ineqs_sos,
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
    # the variables
    vars 
    # each variable x must satify x^2 ≤ squared_variable_bound
    squared_variable_bound 
    # the equalities defining our set
    eqs
    # the inequalities defining our set
    ineqs
    # the main rational SOS part
    sos
    # rational SOS for the variable bounds
    offset_sos
    # polynomial coefficients for the equalities
    eqs_polys
    # the SOS for the inequalities
    ineqs_sos 
    # a negative rational number
    left_hand_side 
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
    offset_part = sum((rpc.squared_variable_bound-t^2)*polynomial(p)
                      for (t,p) in zip(rpc.vars,rpc.offset_sos))
    eqs_part = sum(q*p
                   for (q,p) in zip(rpc.eqs, rpc.eqs_polys))
    ineqs_part = sum(q*polynomial(p)
                     for (q,p) in zip(rpc.ineqs, rpc.ineqs_sos))
    
    right_hand_side = sos_part + offset_part + eqs_part + ineqs_part
    if (right_hand_side- rpc.left_hand_side)==0
        println("Left hand side = right hand side 😄")
    else
        println("Left hand side ≠ right hand side 🥲")
    end
    return (;sos_part, offset_part, eqs_part,ineqs_part,right_hand_side,rpc.left_hand_side)
end

function print_certificate(rpc; io=stdout)
    println(io,"(")
    println(io,rpc.sos)
    for i in eachindex(rpc.vars)
        println(io,"+(",(rpc.squared_variable_bound-rpc.vars[i]^2),")*(",rpc.offset_sos[i],")")
    end
    for i in eachindex(rpc.eqs_polys)
        println(io,"+(",rpc.eqs[i],")*(",rpc.eqs_polys[i],")")
    end
    for i in eachindex(rpc.ineqs)
        println(io,"+(",rpc.ineqs[i],")*(",rpc.ineqs_sos[i],")")
    end
    println(io,")")
end

function print_certificate_julia(rpc; io=stdout)
    println(io,"using DynamicPolynomials")
    println(io,"@polyvar p[1:3,1:3]")
    println(io,"const PolyType = Polynomial{DynamicPolynomials.Commutative{DynamicPolynomials.CreationOrder}, Graded{LexOrder}, Rational{BigInt}}")
       
    println(io,"println(\"Filling offset vectors\")")
    println(io,"offset_ineqs = Vector{PolyType}()")
    println(io,"offset_ineqs_sos =  Vector{PolyType}()")

    for i in eachindex(rpc.vars)
        println(io,"push!(offset_ineqs,",(rpc.squared_variable_bound-rpc.vars[i]^2),")")
        print(io,"push!(offset_ineqs_sos,")
        show_big(io,rpc.offset_sos[i])
        println(io,")")
    end
    println(io,"offset_ineq_result=sum(offset_ineqs[i]*offset_ineqs_sos[i] for i in eachindex(offset_ineqs))")

    
    println(io,"println(\"Filling equation vectors\")")
    println(io,"eqs = Vector{PolyType}()")
    println(io,"eqs_polys =  Vector{PolyType}()")
    for i in eachindex(rpc.eqs)
        print(io, "push!(eqs,")
        show_big(io, rpc.eqs[i])
        println(io, ")")
        print(io,"push!(eqs_polys,")
        show_big(io, rpc.eqs_polys[i])
        println(io, ")")
    end

    println(io,"println(\"Filling inequalitiy vectors\")")
    println(io,"ineqs = Vector{PolyType}()")
    println(io,"ineqs_sos =  Vector{PolyType}()")
    for i in eachindex(rpc.ineqs)
        print(io, "push!(ineqs,")
        show_big(io, rpc.ineqs[i])
        println(io, ")")
        print(io,"push!(ineqs_sos,")
        show_big(io, rpc.ineqs_sos[i])
        println(io, ")")
    end
    
    println(io,"println(\"Calculating sos_part\")")
    print(io,"sos_part=")
    show_big(io, rpc.sos)
    println(io,"")
    
    println(io,"println(\"Calculating offset_part\")")
    println(io,"offset_part=sum(offset_ineqs[i]*offset_ineqs_sos[i] for i in eachindex(offset_ineqs))")

    println(io,"println(\"Calculating eqs_part\")")
    println(io,"eqs_part=sum(eqs[i]*eqs_polys[i] for i in eachindex(eqs))")
    
    println(io,"println(\"Calculating ineqs_part\")")
    println(io,"ineqs_part=sum(ineqs[i]*ineqs_sos[i] for i in eachindex(ineqs))")
    
    println(io,"println(sos_part+offset_part+eqs_part+ineqs_part)")
end


function show_big(io::IO, a :: Integer)
    print(io, "big(",a,")")
end

function show_big(io::IO, a :: Rational)
    show_big(io, numerator(a))
    print(io, "//")
    show_big(io, denominator(a))
end

function show_big(io::IO, a :: DynamicPolynomials.Polynomial)
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

function show_big(io::IO, a :: DynamicPolynomials.Term)
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

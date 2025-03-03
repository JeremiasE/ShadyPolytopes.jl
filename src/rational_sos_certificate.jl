import MultivariatePolynomials

"""
    gram_to_sos(RG, monos)

Takes a rational gram matrix and turns it into a rational weighted SOS
decomposition.
"""
function gram_to_sos(RG, monos)
    BK = bunchkaufman(RG)
    d = diag(BK.D)
    L = BK.U'*BK.P
    return WeightedSOSDecmposition(L*monos, d)
end

"""
    round_gram_matrix(gram, obj, degree,vars;prec)

Takes an approximate gram matrix and produces a rational gram_matrix
which represents `obj` as a weighted SOS decomposition using monomials
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
    reduce_one_minus_monomial(monomialm, vars)

If we set ``r := reduce_square_monomial(s,vars)``
then ``1-s = \\sum_i r_i*(1-vars_i)``.
"""
function reduce_one_minus_monomial(monomialm, vars)
    contained_vars = effective_variables(monomialm)
    t = contained_vars[1]
    if degree(monomialm) == 1
        m = zeros(typeof(zero(Rational{BigInt}(1)*vars[1])), length(vars))
        m[vars .== t] .= 1
        return m
    end
    # 1-xy = (1-x) + x(1-y)
    r = reduce_one_minus_monomial(monomials(monomialm ÷ t)[1], vars)
    q = r * t
    q[vars .== t] = q[vars.==t] .+ 1
    return q
end


"""
    offset_sum_of_monomials(mono_vector,vars)

If `mono_vector` = (m_1,...,m_k) and `vars`=x_1,...,x_l
and the result is p_1,...,p_l,
then ``k - m_1^2+m_2^2+...+m_k^2`` as ``p_1 (1-x_1^2)+...+ p_l (1-x_l^2)``.
"""
function offset_sum_of_monomials(mono_vector,vars)
    n = length(mono_vector)
    p = sum(reduce_one_minus_monomial(m, vars) for m in mono_vector)
    println(typeof(p))
    println(typeof(polynomial.(Rational{BigInt}(1)*collect(monomials(p[1])))))
    println(typeof(coefficients(p[1])))
     
    return [WeightedSOSDecomposition(Rational{BigInt}(1)*collect(monomials(monosum)),
                                     coefficients(monosum))
            for monosum in p]
end

"""
    
"""
function rationally_reduce_sos_decomp(model, K, obj; feasibility=true, prec = 10^7)
    ie = SemialgebraicSets.inequalities(K)
    indices = eachindex(ie)
    M = lagrangian_multipliers(model[:c])
    sos = [fast_round_sos(SOSDecomposition(M[i]), prec)*ie[i] for i in indices]
    sos_decomp = [fast_round_sos(SOSDecomposition(M[i]), prec) for i in indices]
    if feasibility
        return obj-polynomial(sos), sos_decomp
    end
    # obj = τ + sos + sum(...)
    τ = fast_round(value(model[:τ]), prec)
    return obj-τ-polynomial(sos), sos_decomp
end

function round_sos_decomposition(model, K, obj, vars, offset=1//10^4; prec=10^2)
    given_ineq = SemialgebraicSets.inequalities(K)
    newobj, given_ineq_sos = rationally_reduce_sos_decomp(model, K, obj; prec=prec, feasibility=true)
    given_ineq_part = sum(polynomial(given_ineq_sos[i])*given_ineq[i] for i in eachindex(given_ineq))
    
    newobj_degree = maximum(map(degree, terms(newobj)))
    approx_gram = gram_matrix(model[:c])
    gram_degree = maximum(map(degree, collect(approx_gram.basis)))
    comb_degree = max(ceil(Int64, newobj_degree / 2), gram_degree)   
    RG, new_monos = round_gram_matrix(approx_gram, newobj, comb_degree,vars; prec=prec)
    # nm' RG nm = newobj
    pv, w = gram_to_sos(RG+offset*I, new_monos)
    sos = WeightedSOSDecomposition(pv,w)
    sos_part = polynomial(sos)
    
    #return new_monos,soss,ineqs,offset,vars
    offset_ineq_sos = offset_sum_of_monomials(new_monos[2:end],vars)
    println(typeof(offset_ineq_sos))
    for i in eachindex(vars)
        offset_ineq_sos[i].wv *= offset
    end
    k = length(new_monos)
    #offsetpart = k-newmonos'*newmonos
    offset_part = (sum((1-t^2)*polynomial(p) for (p,t) in zip(offset_ineq_sos,vars)))
    println(sos_part + offset_part + given_ineq_part)
    return sos, offset_ineq_sos, given_ineq_sos, given_ineq
    # newobj = -1 - ineq_part = sos_part + offset_part 
    # return  sos_part + offset_part + ineq_part
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

function show_big(io::IO, d :: WeightedSOSDecomposition)
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

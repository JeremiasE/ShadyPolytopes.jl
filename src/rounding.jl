using DynamicPolynomials
using SumOfSquares

"""
    threshold(x,tol)

Return 0 if the absolute value of `x` is below tol,
otherwise return `x`.
"""
function threshold(x, tol)
    return abs(x) > tol ? x : 0
end

"""
    max_coeff(poly)

Find the maximal coefficient in the polynomial `poly`.
This is a fast way to check if a polynomial is
approximately 0, so it can quantify the
quality of an approximation of a polynomial.
"""
function max_coeff(poly)
    return maximum(coefficients(poly))
end

"""
    fast_round(x[, n])

Approximate a number by a rational
with numerator `n`.
"""
function fast_round(x, n=10^7)
    return trunc(BigInt, x * n)//n
end

"""
    round_poly(p[, n])

Approximate a floating point polynomial `p` by a rational polynomial whose coefficients
    have numerator `n`
"""
function round_poly(p, n=10^7)
    return mapreduce(*, +, map(x -> fast_round(x, n),
                     coefficients(p)),
                     monomials(p))
end

"""
    round_sos(p[, n])

Approximate the polynomials in an SOS decomposition `sos`
by rationals with numerator `n`.
"""
function round_sos(sos, n=10^7)
    return SOSDecomposition([round_poly(p, n) for p in sos])
end

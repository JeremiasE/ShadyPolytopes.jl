import MultivariatePolynomials as MP
using DynamicPolynomials

"""
    WeightedSOSDecomposition(pv, wv)

Takes a vector of polynomials `pv` and a vector of (positive) weights `wv`
to represent the SOS
``\\sum_{i} wv_i \\cdot pv_i^2``.

"""

mutable struct WeightedSOSDecomposition{T}
    pv :: Vector{<:MP.AbstractPolynomialLike{T}}
    wv :: Vector{T}
end


function MP.polynomial(d::WeightedSOSDecomposition)
    return sum(d.wv[i]*d.pv[i]^2 for i in eachindex(d.pv))
end

function Base.show(io::IO, d::WeightedSOSDecomposition)
    for (i, p, w) in zip(eachindex(d.pv),d.pv,d.wv)
        print(io, w)
        print(io, "*")
        print(io, "(")
        print(io, p)
        print(io, ")^2")
        if i != length(d.pv)
            print(io, " + ")
        end
    end
end

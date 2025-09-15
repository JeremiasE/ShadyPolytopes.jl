using Clarabel
using ShadyPolytopes

@time begin
    L, K,
    vars,
    squared_variable_bound,
    obj,
    model = solve_via_projection_matrix(
        optimal_icosahedron, Clarabel.Optimizer; maxdegree = 5, bound = 101//100
    )

    rpc = ShadyPolytopes.round_sos_decomposition(
        model, K, obj, vars, squared_variable_bound; prec = big(10^5)
    )

    datapath = abspath(joinpath(@__DIR__, "..", "data"))
    filename = joinpath(datapath, "sos_certificate.jl")

    open(filename, "w") do f
        ShadyPolytopes.print_certificate_julia(rpc; io = f)
    end
end

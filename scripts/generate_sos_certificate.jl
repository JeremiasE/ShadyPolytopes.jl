using Clarabel
using ShadyPolytopes

L,K,vars,squared_variable_bound, obj,model = solve_via_projection_matrix(
    optimal_icosahedron,
    Clarabel.Optimizer;
    maxdegree = 5,
    bound = 101//100)

rpc = ShadyPolytopes.round_sos_decomposition(model,K,obj,vars,squared_variable_bound,1//10^6; prec=10^4)

datapath = abspath(joinpath(@__DIR__,"..","data"))
filename = joinpath(datapath,"sos-certificate.jl")

open(filename,"w") do f
    print_certificate_julia(rpc; io = f)
end

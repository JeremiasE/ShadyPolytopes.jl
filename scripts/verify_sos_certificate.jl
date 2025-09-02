@time begin
    datapath = abspath(joinpath(@__DIR__, "..", "data"))
    filename = joinpath(datapath, "sos_certificate.jl")
    
    include(filename)
    
    println("[Calculating offset_part]")
    offset_part = sum(offset_ineqs[i]*offset_ineqs_sos[i] for i in eachindex(offset_ineqs))
    
    println("[Calculating eqs_part]")
    eqs_part = sum(eqs[i]*eqs_polys[i] for i in eachindex(eqs))
    
    println("[Calculating ineqs_part]")
    ineqs_part = sum(ineqs[i]*ineqs_sos[i] for i in eachindex(ineqs))
    
    println("Result: ", sos_part+offset_part+eqs_part+ineqs_part)
end

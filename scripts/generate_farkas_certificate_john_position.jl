using ShadyPolytopes

datapath = abspath(joinpath(@__DIR__,"..","data"))

Threads.@threads for i=1:3
    filename = joinpath(datapath,"farkas-certificates-opt-ico-jp-$(i).csv")
    open(filename,"w") do f
        generate_farkas_certificates(optimal_icosahedron_john_position, 1400, 84//83, i; io=f)
    end
end

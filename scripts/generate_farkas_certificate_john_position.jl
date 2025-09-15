using ShadyPolytopes
using CodecZlib

println("Using ", Threads.nthreads(), " threads")

@time begin
    datapath = abspath(joinpath(@__DIR__, "..", "data"))
    k = 1400
    α = 84//83
    
    Threads.@threads for i in 1:3
        filename = joinpath(
            datapath,
            "farkas-certificates-opt-ico-jp-$(k)-$(numerator(α))_$(denominator(α))-$(i).csv.gz",
        )
        open(GzipCompressorStream, filename, "w") do f
            generate_farkas_certificates(optimal_icosahedron_john_position, k, α, i; io=f)
        end
    end
end

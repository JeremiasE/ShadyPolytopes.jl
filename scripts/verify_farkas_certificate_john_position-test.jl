using ShadyPolytopes
using CodecZlib

println("Using ", Threads.nthreads(), " threads")

@time begin
    datapath = abspath(joinpath(@__DIR__, "..", "data"))
    k = 11
    α = 84//83
    @Threads.threads for i in 1:3
        filename = joinpath(
            datapath,
            "farkas-certificates-opt-ico-jp-$(k)-$(numerator(α))_$(denominator(α))-$(i).csv.gz",
        )
        check_farkas_certificate_file(optimal_icosahedron_john_position, filename, α; silent=false)
    end
end

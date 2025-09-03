using ShadyPolytopes
using CodecZlib

@time begin
    datapath = abspath(joinpath(@__DIR__, "..", "data"))
    k = 1400
    α = 84//83
    
    for i in 1:3
        filename = joinpath(
            datapath,
            "farkas-certificates-opt-ico-jp-$(k)-$(numerator(α))_$(denominator(α))-$(i).csv.gz",
        )
        open(GzipDecompressorStream, filename, "r") do f
            check_farkas_certificate_file(optimal_icosahedron_john_position, f, α; silent=false)
        end
    end
end

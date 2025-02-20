using ShadyPolytopes
using Test

@testset "ShadyPolytopes.jl" begin
    @testset "Single Farkas Certificate" begin
        w =  [5//10,1//10,3//10]
        α = 84//83
        indices, yentries = find_single_farkas_certificate(optimal_icosahedron, w, α)
        n = length(optimal_icosahedron.positive_vertices) * length(optimal_icosahedron.normals)
        y = zeros(Rational{BigInt},n)
        y[indices]=yentries
        @test_nowarn check_single_farkas_certificate(optimal_icosahedron, w, y, α)
        @test_throws ErrorException find_single_farkas_certificate(optimal_icosahedron, w, 105//100)
    end
    @testset "Set of Farkas Certificates" begin
        α = 84//83
        buf = IOBuffer()
        generate_farkas_certificates(optimal_icosahedron, 3, α; io=buf)
        @test_nowarn check_farkas_cerfificate_file(optimal_icosahedron, take!(buf), α) 
    end
    
end

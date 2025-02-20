using ShadyPolytopes
using Test

@testset "ShadyPolytopes.jl" begin
    @testset "Farkas Certificate" begin
        w =  [5//10,1//10,3//10]
        α = 84//83
        indices, yentries = find_single_farkas_certificate(optimal_icosahedron, w, α)
        n = length(optimal_icosahedron.positive_vertices) * length(optimal_icosahedron.normals)
        y = zeros(Rational{BigInt},n)
        y[indices]=yentries
        @test_nowarn check_single_farkas_certificate(optimal_icosahedron, w, y, α)
    end
end

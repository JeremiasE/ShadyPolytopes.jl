using ShadyPolytopes
using Test

include("test_sos_models.jl")

@testset "ShadyPolytopes.jl" begin
    @testset "Single Farkas Certificate" begin
        w = [5//10, 1//10, 3//10]
        α = 84//83
        indices, yentries = find_single_farkas_certificate(optimal_icosahedron, w, α)
        n =
            length(optimal_icosahedron.positive_vertices) *
            length(optimal_icosahedron.normals)
        y = zeros(Rational{BigInt}, n)
        y[indices] = yentries
        @test_nowarn check_single_farkas_certificate(optimal_icosahedron, w, y, α)
        @test_throws ErrorException find_single_farkas_certificate(
            optimal_icosahedron, w, 105//100
        )
    end
    @testset "Set of Farkas Certificates" begin
        α = 84//83
        buf = IOBuffer()
        generate_farkas_certificates(optimal_icosahedron, 3, α; io=buf)
        @test_nowarn check_farkas_certificate_file(optimal_icosahedron, take!(buf), α)
    end
    @testset "Simple SOS Models" begin
        L, K, vars, obj, model = optimize_tiny_sos_model(true, false; bound=-174//100)
        @test is_solved_and_feasible(model)

        L, K, vars, obj, model = optimize_tiny_sos_model(true, false; bound=-173//100)
        @test !is_solved_and_feasible(model)

        L, K, vars, obj, model = optimize_tiny_sos_model(false, false)
        @test is_solved_and_feasible(model)
        @test value(model[:τ]) ≤ -sqrt(3) + 0.001

        L, K, vars, obj, model = optimize_tiny_sos_model(true, true; bound=-1//100)
        @test is_solved_and_feasible(model)

        L, K, vars, obj, model = optimize_tiny_sos_model(false, true)
        @test is_solved_and_feasible(model)
        @test value(model[:τ]) ≥ -0.00001
    end
    @testset "Lifting polynomials" begin
        L, K, dp_vars = shadiness_via_projection_matrix(optimal_icosahedron)
        eqs = copy(SemialgebraicSets.equalities(L))
        I = SemialgebraicSets.ideal(L)
        SemialgebraicSets.compute_gröbner_basis!(I)
        new_eqs = SemialgebraicSets.equalities(L)
        A, B = lift_polynomials(eqs, new_eqs, dp_vars)
        @test A' * eqs == new_eqs
    end
end

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
    @testset "Analysing norm balls" begin
        P, norm = find_shadiness_constant_numerically(optimal_icosahedron)
        @test norm ≈ 1.0127 atol=1e-4
        @test P^2 ≈ P atol=1e-4
        on = operator_norm(optimal_icosahedron,P)
        @test 1.01 < on < 1.02
        PQ = approximate_projection(P)
        @test P ≈ PQ atol=1e-04
        @test PQ^2 ≈ PQ atol=1e-04
        
        @test Set(calculate_normals(cube.positive_vertices))==Set(cube.normals)

        extremal_matrices = determine_matrix_ball_vertices(cube, rational=true)
        @test length(extremal_matrices) == 216
        @test determine_squared_spectralnorm_bound(cube) == 3
        A = convert(Array{Rational{BigInt}},[1 4 3; 5 3 5; 1 2 7])
        B = approximate_john_position(A*cube; prec=10^4)
        @test determine_squared_spectralnorm_bound(B*A*cube) ≈ 3 atol=1e-2
    end
    @testset "Farkas certificates" begin
        w = [1,2,3]
        α = 84//83
        A = matrix_from_projection_normal(optimal_icosahedron, w, α)
        n = size(A,1)
        indices, y_entries = find_single_farkas_certificate(
            optimal_icosahedron, w, α
        )
        y = zeros(Rational{BigInt}, n)
        y[indices] = y_entries
        @test all(y.>=0)
        @test A'*y == w
    end
    
end

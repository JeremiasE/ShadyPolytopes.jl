module ShadyPolytopes

include("cscp.jl")
include("rounding.jl")
include("analyse_ball.jl")
include("perturbed_icosahedron.jl")
include("farkas.jl")
include("sos_models.jl")
include("rational_sos_decomposition.jl")
include("rational_sos_certificate.jl")
include("lift_polynomials.jl")

export fast_round
export CSCP, calculate_normals, all_vertices, cube
export optimal_icosahedron, optimal_icosahedron_john_position, nice_icosahedron
export find_shadiness_constant_numerically, determine_squared_spectralnorm_bound
export determine_matrix_ball_vertices
export approximate_john_position, approximate_projection, operator_norm
export find_single_farkas_certificate, matrix_from_projection_normal
export generate_farkas_certificates,
    check_single_farkas_certificate, check_farkas_certificate_file
export putinar_model, solve_via_projection_matrix, solve_via_vw
export shadiness_via_projection_matrix
export shadiness_via_vw
export lift_polynomials
export max_coeff
export RationalSOSDecomposition

end

module ShadyPolytopes

include("cscb.jl")
include("rounding.jl")
include("analyse_ball.jl")
include("perturbed_icosahedron.jl")
include("farkas.jl")
include("sos_models.jl")

export fast_round
export CSCB, calculate_normals, all_vertices
export optimal_icosahedron, optimal_icosahedron_john_position, nice_icosahedron
export determine_matrix_ball_vertices, determine_squared_spectralnorm_bound
export approximate_john_position, approximate_projection, operator_norm
export find_single_farkas_certificate, matrix_from_projection_normal
export generate_farkas_certificates, check_single_farkas_certificate, check_farkas_cerfificate_file
export putinar_model, solve_via_projection_matrix, solve_via_vw
export shadiness_via_projection_matrix
export shadiness_via_vw

end

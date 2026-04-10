module DecoratedRandomPlanarMaps

using LinearAlgebra
using SparseArrays
using Random
using Statistics
using JSON3
using YAML

include("core/Base.jl")
include("core/HalfEdge.jl")
include("core/Graph.jl")
include("core/Samplers.jl")
include("core/Timing.jl")

include("models/Uniform.jl")
include("models/Schnyder.jl")
include("models/FK.jl")

include("render/Common.jl")
include("render/STL.jl")
include("render/Web.jl")
include("render/SVG.jl")

include("layout/Tutte.jl")
include("layout/SFDP.jl")
include("layout/Problem.jl")

include("Pipeline.jl")

export AbstractRandomPlanarMap, num_vertices, num_faces, layout_edges
export HalfEdge, PlanarMap, assign_vertex_ids!, extract_faces,
    WeightedGraph, build_weighted_graph, collapse_undirected_edges, sanitize_edge_array,
    sample_uniform_dyck_path, sample_uniform_primitive_dyck_path,
    TimingRecorder, StepTiming, track!, write_json, pretty_lines
export UniformMap, UniformQuadrangulationMap, UniformQuadrangulation, generate_uniform_map, build_uniform_quadrangulation_map
export SchnyderMap, SchnyderWoodMap, generate_schnyder_map
export FKMap, FKDecoratedMap, HCMap, SpanningTreeMap, build_fk_map, generate_fk_map, build_spanning_tree_map, build_fk_map_from_word
export build_hc_map, build_hc_map_from_word, sample_fk_word, sample_hc_word,
    sample_fk_word_exact_rejection, sample_hc_word_exact_rejection,
    estimate_fk_exact_rejection_acceptance, recommended_fk_approx_sampler_params
export LayoutProblem, prepare_layout_problem
export DirichletSolveInfo, solve_dirichlet_laplacian, harmonic_measure_boundary_positions
export compute_tutte_layout, compute_sfdp_layout
export EDGE_COLOR_HINTS, sanitize_triangles, fan_triangulate_faces, surface_triangles,
    grouped_edges, pad_positions_3d, should_include_exploration, exploration_segment_points,
    fk_q_from_p, model_parameter_pairs, model_parameter_string, model_display_name, default_render_title
export export_stl_binary, export_web_binaries, export_svg_preview
export render_makie_2d, render_makie_3d, render_2d, render_3d
export load_config, build_map_from_config, run_pipeline

const VERSION = v"0.3.3"

render_2d(args...; kwargs...) = render_makie_2d(args...; kwargs...)
render_3d(args...; kwargs...) = render_makie_3d(args...; kwargs...)

function render_makie_2d(args...; kwargs...)
    error("Makie rendering lives in the optional GLMakie extension. Install GLMakie and GeometryBasics, load them, and then call `render_makie_2d` again.")
end

function render_makie_3d(args...; kwargs...)
    error("Makie rendering lives in the optional GLMakie extension. Install GLMakie and GeometryBasics, load them, and then call `render_makie_3d` again.")
end

end

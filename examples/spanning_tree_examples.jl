using DecoratedRandomPlanarMaps

function timed_step(f, label)
    t0 = time()
    value = f()
    println(rpad(string(label) * ":", 22), round(time() - t0; digits=3), " s")
    return value
end

# Comment out any block you do not want.
# For spanning-tree maps the default 2D boundary mode is `face`.

faces = 2000
seed = 17
outroot = joinpath(@__DIR__, "out", "spanning_tree")
mkpath(outroot)

meta_base = Dict(
    "model" => "spanning_tree",
    "variant" => "spanning_tree",
    "faces" => faces,
    "seed" => seed,
    "p" => 0.0,
    "q" => 0.0,
)

map_data = timed_step(() -> build_spanning_tree_map(; faces=faces, seed=seed), "generate_map")

# ------------------------------------------------------------------
# 2D Tutte layout with a triangular fixed boundary (equilateral outer face)
# ------------------------------------------------------------------
problem2 = timed_step(() -> prepare_layout_problem(map_data; dimension=2, options=Dict("hc_boundary_mode" => "face")), "prepare_2d")
pos2, layout2_meta = timed_step(() -> compute_tutte_layout(
        problem2.num_vertices,
        problem2.edges,
        problem2.boundary_vertices;
        boundary_positions=problem2.boundary_positions,
        solver="auto",
        return_metadata=true,
    ), "layout_2d")
meta2 = merge(meta_base, problem2.metadata, layout2_meta)

timed_step(() -> export_svg_preview(
        problem2.render_map_data,
        pos2,
        joinpath(outroot, "spanning_tree_2d.svg");
        edge_groups=problem2.edge_groups,
        faces=problem2.faces,
        triangles=problem2.surface_triangles,
        metadata=meta2,
        title=default_render_title(problem2.render_map_data; metadata=meta2, dimension=2, engine="Tutte"),
    ), "export_svg")

timed_step(() -> export_web_binaries(
        problem2.render_map_data,
        pos2,
        joinpath(outroot, "spanning_tree_2d_web");
        edge_groups=problem2.edge_groups,
        faces=problem2.faces,
        triangles=problem2.surface_triangles,
        metadata=meta2,
    ), "export_web")

# ------------------------------------------------------------------
# 3D SFDP layout
# ------------------------------------------------------------------
problem3 = timed_step(() -> prepare_layout_problem(map_data; dimension=3), "prepare_3d")
pos3 = timed_step(() -> compute_sfdp_layout(problem3.num_vertices, problem3.edges; scale=1.0, seed=seed), "layout_3d")
meta3 = merge(meta_base, problem3.metadata)

timed_step(() -> export_web_binaries(
        problem3.render_map_data,
        pos3 .* 10.0,
        joinpath(outroot, "spanning_tree_3d_web");
        edge_groups=problem3.edge_groups,
        faces=problem3.faces,
        triangles=problem3.surface_triangles,
        metadata=meta3,
    ), "export_web_3d")

timed_step(() -> export_stl_binary(map_data, pos3 .* 100.0, joinpath(outroot, "spanning_tree_3d.stl")), "export_stl")

# ------------------------------------------------------------------
# Optional Makie
# ------------------------------------------------------------------
# using GLMakie
# using GeometryBasics
# render_makie_2d(
#     problem2.render_map_data,
#     pos2;
#     edge_groups=problem2.edge_groups,
#     faces=problem2.faces,
#     triangles=problem2.surface_triangles,
#     metadata=meta2,
# )
# render_makie_3d(
#     problem3.render_map_data,
#     pos3;
#     edge_groups=problem3.edge_groups,
#     faces=problem3.faces,
#     triangles=problem3.surface_triangles,
#     metadata=meta3,
# )

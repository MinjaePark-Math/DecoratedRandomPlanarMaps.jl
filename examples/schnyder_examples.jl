using DecoratedRandomPlanarMaps

function timed_step(f, label)
    t0 = time()
    value = f()
    println(rpad(label * ":", 22), round(time() - t0; digits=3), " s")
    return value
end

# Comment out any block you do not want.

faces = 10000
seed = 17
outroot = joinpath(@__DIR__, "out", "schnyder")
mkpath(outroot)

meta_base = Dict("model" => "schnyder", "faces" => faces, "seed" => seed)
map_data = timed_step("generate_map") do
    generate_schnyder_map(; size_faces=faces, seed=seed)
end

# ------------------------------------------------------------------
# 2D Tutte layout
# ------------------------------------------------------------------
problem2 = timed_step("prepare_2d") do
    prepare_layout_problem(map_data; dimension=2)
end
pos2, layout2_meta = timed_step("layout_2d") do
    compute_tutte_layout(
        problem2.num_vertices,
        problem2.edges,
        problem2.boundary_vertices;
        boundary_positions=problem2.boundary_positions,
        solver="auto",
        return_metadata=true,
    )
end
meta2 = merge(meta_base, problem2.metadata, layout2_meta)

timed_step("export_svg") do
    export_svg_preview(
        problem2.render_map_data,
        pos2,
        joinpath(outroot, "schnyder_2d.svg");
        edge_groups=problem2.edge_groups,
        faces=problem2.faces,
        triangles=problem2.surface_triangles,
        metadata=meta2,
        title=default_render_title(problem2.render_map_data; metadata=meta2, dimension=2, engine="Tutte"),
    )
end

timed_step("export_web") do
    export_web_binaries(
        problem2.render_map_data,
        pos2,
        joinpath(outroot, "schnyder_2d_web");
        edge_groups=problem2.edge_groups,
        faces=problem2.faces,
        triangles=problem2.surface_triangles,
        metadata=meta2,
    )
end

# ------------------------------------------------------------------
# 3D SFDP layout
# ------------------------------------------------------------------
problem3 = timed_step("prepare_3d") do
    prepare_layout_problem(map_data; dimension=3)
end
pos3 = timed_step("layout_3d") do
    compute_sfdp_layout(problem3.num_vertices, problem3.edges; scale=1.0, seed=seed)
end
meta3 = merge(meta_base, problem3.metadata)

timed_step("export_web_3d") do
    export_web_binaries(
        problem3.render_map_data,
        pos3 .* 10.0,
        joinpath(outroot, "schnyder_3d_web");
        edge_groups=problem3.edge_groups,
        faces=problem3.faces,
        triangles=problem3.surface_triangles,
        metadata=meta3,
    )
end

timed_step("export_stl") do
    export_stl_binary(map_data, pos3 .* 100.0, joinpath(outroot, "schnyder_3d.stl"))
end

# ------------------------------------------------------------------
# Optional Makie
# ------------------------------------------------------------------
using GLMakie
using GeometryBasics
render_makie_2d(
    problem2.render_map_data,
    pos2;
    edge_groups=problem2.edge_groups,
    faces=problem2.faces,
    triangles=problem2.surface_triangles,
    metadata=meta2,
)
# render_makie_3d(
#     problem3.render_map_data,
#     pos3;
#     edge_groups=problem3.edge_groups,
#     faces=problem3.faces,
#     triangles=problem3.surface_triangles,
#     metadata=meta3,
# )

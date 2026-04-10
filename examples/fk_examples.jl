using DecoratedRandomPlanarMaps

function timed_step(f, label)
    t0 = time()
    value = f()
    println(rpad(label * ":", 22), round(time() - t0; digits=3), " s")
    return value
end

# Comment out any block you do not want.
# Change `hc_boundary_mode` to "c_gasket" if you want the c-gasket instead.

faces = 10000
seed = 11
p = 0.25
outroot = joinpath(@__DIR__, "out", "fk")
mkpath(outroot)

sampling_method = :auto
pool_size = nothing           # e.g. 128
mh_steps = nothing            # e.g. 512

sfdp_K = nothing              # e.g. 0.35
sfdp_repulsiveforce = nothing # e.g. 1.2
sfdp_iterations = nothing     # e.g. 300

map_data = timed_step("generate_map") do
    build_fk_map(
        ;
        faces=faces,
        p=p,
        seed=seed,
        sampling_method=sampling_method,
        pool_size=pool_size,
        mh_steps=mh_steps,
    )
end

# ------------------------------------------------------------------
# 2D Tutte layout on the maximal h-gasket
# ------------------------------------------------------------------
problem2 = timed_step("prepare_2d") do
    prepare_layout_problem(map_data; dimension=2, options=Dict("hc_boundary_mode" => "h_gasket"))
end

pos2, meta2 = timed_step("layout_2d") do
    compute_tutte_layout(
        problem2.num_vertices,
        problem2.edges,
        problem2.boundary_vertices;
        boundary_positions=problem2.boundary_positions,
        solver="auto",
        return_metadata=true,
    )
end

meta2_full = merge(problem2.metadata, meta2, Dict("faces" => faces, "seed" => seed, "p" => p, "q" => fk_q_from_p(p)))

timed_step("export_svg") do
    export_svg_preview(
        problem2.render_map_data,
        pos2,
        joinpath(outroot, "fk_h_gasket_2d.svg");
        edge_groups=problem2.edge_groups,
        faces=problem2.faces,
        triangles=problem2.surface_triangles,
        metadata=meta2_full,
        title="fk · 2D h-gasket",
    )
end

timed_step("export_web") do
    export_web_binaries(
        problem2.render_map_data,
        pos2,
        joinpath(outroot, "fk_h_gasket_2d_web");
        edge_groups=problem2.edge_groups,
        faces=problem2.faces,
        triangles=problem2.surface_triangles,
        metadata=meta2_full,
    )
end

# ------------------------------------------------------------------
# 3D SFDP layout on the full decorated triangulation
# ------------------------------------------------------------------
problem3 = timed_step("prepare_3d") do
    prepare_layout_problem(map_data; dimension=3)
end

pos3 = timed_step("layout_3d") do
    compute_sfdp_layout(
        problem3.num_vertices,
        problem3.edges;
        scale=1.0,
        seed=seed,
        K=sfdp_K,
        repulsiveforce=sfdp_repulsiveforce,
        iterations=sfdp_iterations,
    )
end

meta3 = merge(problem3.metadata, Dict("faces" => faces, "seed" => seed, "p" => p, "q" => fk_q_from_p(p)))
sfdp_K !== nothing && (meta3["K"] = sfdp_K)
sfdp_repulsiveforce !== nothing && (meta3["repulsiveforce"] = sfdp_repulsiveforce)
sfdp_iterations !== nothing && (meta3["iterations"] = sfdp_iterations)

timed_step("export_web_3d") do
    export_web_binaries(
        problem3.render_map_data,
        pos3 .* 10.0,
        joinpath(outroot, "fk_3d_web");
        edge_groups=problem3.edge_groups,
        faces=problem3.faces,
        triangles=problem3.surface_triangles,
        metadata=meta3,
    )
end

timed_step("export_stl") do
    export_stl_binary(map_data, pos3 .* 100.0, joinpath(outroot, "fk_3d.stl"))
end

# ------------------------------------------------------------------
# Optional Makie
# ------------------------------------------------------------------
using GLMakie
using GeometryBasics
# render_makie_2d(
#     problem2.render_map_data,
#     pos2;
#     edge_groups=problem2.edge_groups,
#     faces=problem2.faces,
#     triangles=problem2.surface_triangles,
#     metadata=meta2_full,
#     title="fk · 2D Makie",
# )

render_makie_3d(
    problem3.render_map_data,
    pos3;
    edge_groups=problem3.edge_groups,
    faces=problem3.faces,
    triangles=problem3.surface_triangles,
    metadata=meta3,
    title="fk · 3D Makie",
)
using DecoratedRandomPlanarMaps

function timed_step(f, label)
    t0 = time()
    value = f()
    println(rpad(label * ":", 24), round(time() - t0; digits=3), " s")
    return value
end

function projected_render_positions(problem, pos)
    if problem.render_vertex_indices === nothing
        return pos
    end
    return pos[Int.(problem.render_vertex_indices).+1, :]
end

function strand_edge_groups(problem, label)
    groups = problem.edge_groups
    if !haskey(groups, "generic")
        groups = grouped_edges(problem.render_map_data)
    end
    haskey(groups, "generic") || error("$(label) is missing the gray strand backbone edge group")
    println(label, " strand segments: ", size(groups["generic"], 1))
    return groups
end

# Comment out any block you do not want.
# The half-plane example is 2D; the closed-system examples are 3D.

half_plane_order = 20000
boundary_half_length = nothing
closed_order = 20000
meander_order = 10000
seed = 712

outroot = joinpath(@__DIR__, "out", "meandric")
mkpath(outroot)

# ------------------------------------------------------------------
# Half-plane meandric system, 2D Tutte
# ------------------------------------------------------------------

half_plane = timed_step("generate_half_plane") do
    build_half_plane_meandric_system(
        ;
        order=half_plane_order,
        boundary_half_length=boundary_half_length,
        seed=seed,
    )
end

half_plane_problem = timed_step("prepare_half_plane") do
    prepare_layout_problem(half_plane; dimension=2)
end
half_plane_groups = strand_edge_groups(half_plane_problem, "half-plane")

half_plane_pos, half_plane_layout_meta = timed_step("layout_half_plane") do
    compute_tutte_layout(
        half_plane_problem.num_vertices,
        half_plane_problem.edges,
        half_plane_problem.boundary_vertices;
        boundary_positions=half_plane_problem.boundary_positions,
        solver="auto",
        seed=seed,
        return_metadata=true,
    )
end

half_plane_meta = merge(
    Dict(
        "model" => "half_plane_meandric",
        "pairs" => half_plane_order,
        "seed" => seed,
        "boundary_half_length" => Int(half_plane.boundary_half_length),
    ),
    half_plane_problem.metadata,
    half_plane_layout_meta,
)

timed_step("svg_half_plane") do
    export_svg_preview(
        half_plane_problem.render_map_data,
        half_plane_pos,
        joinpath(outroot, "half_plane_tutte.svg");
        edge_groups=half_plane_groups,
        faces=half_plane_problem.faces,
        triangles=half_plane_problem.surface_triangles,
        metadata=half_plane_meta,
        title=default_render_title(half_plane_problem.render_map_data; metadata=half_plane_meta, dimension=2, engine="Tutte") * " · gray strand backbone",
    )
end

timed_step("web_half_plane") do
    export_web_binaries(
        half_plane_problem.render_map_data,
        half_plane_pos,
        joinpath(outroot, "half_plane_tutte_web");
        edge_groups=half_plane_groups,
        faces=half_plane_problem.faces,
        triangles=half_plane_problem.surface_triangles,
        metadata=half_plane_meta,
    )
end

# ------------------------------------------------------------------
# Uniform meandric system, 3D SFDP via glued trees
# ------------------------------------------------------------------

closed_system = timed_step("generate_closed") do
    build_uniform_meandric_system(; order=closed_order, seed=seed)
end
println("closed-system loops: ", length(closed_system.component_sizes))

closed_problem = timed_step("prepare_closed") do
    prepare_layout_problem(closed_system; dimension=3)
end
closed_groups = strand_edge_groups(closed_problem, "uniform meandric")

closed_aux_pos = timed_step("layout_closed") do
    compute_sfdp_layout(closed_problem.num_vertices, closed_problem.edges; scale=1.0, seed=seed)
end
closed_render_pos = projected_render_positions(closed_problem, closed_aux_pos)

closed_meta = merge(
    Dict("model" => "uniform_meandric", "pairs" => closed_order, "seed" => seed),
    closed_problem.metadata,
)

timed_step("web_closed") do
    export_web_binaries(
        closed_problem.render_map_data,
        closed_render_pos .* 10.0,
        joinpath(outroot, "uniform_meandric_sfdp_web");
        edge_groups=closed_groups,
        faces=closed_problem.faces,
        triangles=closed_problem.surface_triangles,
        metadata=closed_meta,
    )
end

# ------------------------------------------------------------------
# Uniform meander, 3D SFDP after tempered MCMC
# ------------------------------------------------------------------

meander = timed_step("generate_meander") do
    build_uniform_meander(
        ;
        order=meander_order,
        seed=seed,
        sweeps_per_temperature=40,
        search_sweeps=160,
        mixing_sweeps=40,
        restarts=120,
    )
end
println("meander loops: ", length(meander.component_sizes))

meander_problem = timed_step("prepare_meander") do
    prepare_layout_problem(meander; dimension=3)
end
meander_groups = strand_edge_groups(meander_problem, "uniform meander")

meander_aux_pos = timed_step("layout_meander") do
    compute_sfdp_layout(meander_problem.num_vertices, meander_problem.edges; scale=1.0, seed=seed)
end
meander_render_pos = projected_render_positions(meander_problem, meander_aux_pos)

meander_meta = merge(
    Dict("model" => "uniform_meander", "pairs" => meander_order, "seed" => seed),
    meander_problem.metadata,
)

timed_step("web_meander") do
    export_web_binaries(
        meander_problem.render_map_data,
        meander_render_pos .* 10.0,
        joinpath(outroot, "uniform_meander_sfdp_web");
        edge_groups=meander_groups,
        faces=meander_problem.faces,
        triangles=meander_problem.surface_triangles,
        metadata=meander_meta,
    )
end

# ------------------------------------------------------------------
# Optional Makie
# ------------------------------------------------------------------
# using GLMakie
# using GeometryBasics
# render_makie_2d(
#     half_plane_problem.render_map_data,
#     half_plane_pos;
#     edge_groups=half_plane_groups,
#     faces=half_plane_problem.faces,
#     triangles=half_plane_problem.surface_triangles,
#     metadata=half_plane_meta,
#     title="half_plane_meandric · 2D · gray strand backbone",
# )
# render_makie_3d(
#     closed_problem.render_map_data,
#     closed_render_pos;
#     edge_groups=closed_groups,
#     faces=closed_problem.faces,
#     triangles=closed_problem.surface_triangles,
#     metadata=closed_meta,
#     title="uniform_meandric · 3D · gray strand backbone",
# )
# render_makie_3d(
#     meander_problem.render_map_data,
#     meander_render_pos;
#     edge_groups=meander_groups,
#     faces=meander_problem.faces,
#     triangles=meander_problem.surface_triangles,
#     metadata=meander_meta,
#     title="uniform_meander · 3D · gray strand backbone",
# )

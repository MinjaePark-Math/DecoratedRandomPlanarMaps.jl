using Test
using DecoratedRandomPlanarMaps

function _edge_set(edge_array)
    arr = Int32.(edge_array)
    size(arr, 1) == 0 && return Set{Tuple{Int,Int}}()
    out = Set{Tuple{Int,Int}}()
    for i in 1:size(arr, 1)
        u = Int(arr[i, 1])
        v = Int(arr[i, 2])
        push!(out, (min(u, v), max(u, v)))
    end
    return out
end

function _max_triangle_edge_id_incidence(green_edges, diag_edges)
    counts = Dict{Int,Int}()
    for i in 1:size(green_edges, 1)
        for e in green_edges[i, :]
            e < 0 && continue
            counts[Int(e)] = get(counts, Int(e), 0) + 1
        end
        d = Int(diag_edges[i])
        d < 0 && continue
        counts[d] = get(counts, d, 0) + 1
    end
    isempty(counts) && return 0
    return maximum(values(counts))
end

@testset "Prepare Schnyder 2D problem" begin
    m = generate_schnyder_map(; size_faces=6, seed=1)
    problem = prepare_layout_problem(m; dimension=2)

    @test problem.boundary_vertices !== nothing
    @test length(problem.boundary_vertices) == 3
    @test problem.boundary_positions !== nothing
    @test size(problem.boundary_positions) == (3, 2)
    @test problem.faces !== nothing
end

@testset "Prepare uniform 2D problem" begin
    m = generate_uniform_map(; faces=20, seed=3)
    problem = prepare_layout_problem(m; dimension=2)

    @test problem.boundary_vertices !== nothing
    @test length(problem.boundary_vertices) == 4
    @test problem.boundary_positions !== nothing
    @test size(problem.boundary_positions) == (4, 2)
    @test problem.faces !== nothing
    @test problem.surface_triangles !== nothing
    @test problem.surface_triangle_edge_ids !== nothing
    @test size(problem.surface_triangles) == size(problem.surface_triangle_edge_ids)
end

@testset "Pinched uniform quads collapse before triangulation" begin
    collapsed = DecoratedRandomPlanarMaps._pinched_quad_collapsed_edge(Int32[1, 2, 1, 3])
    @test collapsed == Int32[1, 3]

    m = generate_uniform_map(; faces=20, seed=2)
    problem = prepare_layout_problem(m; dimension=2)

    @test get(problem.metadata, "circle_packing_collapsed_degenerate_quads", 0) > 0
    @test (0, 2) in _edge_set(problem.edges)
    @test all(length(unique(problem.surface_triangles[i, :])) == 3 for i in 1:size(problem.surface_triangles, 1))
end

@testset "Sphere circle packing preparation exists for sphere-topology models" begin
    schnyder = prepare_layout_problem(generate_schnyder_map(; size_faces=8, seed=2); dimension=3)
    @test schnyder.packing_boundary_vertices !== nothing
    @test length(schnyder.packing_boundary_vertices) == 3
    @test schnyder.packing_faces !== nothing

    uniform = prepare_layout_problem(generate_uniform_map(; faces=20, seed=4); dimension=3)
    @test uniform.packing_boundary_vertices === nothing
    @test uniform.packing_triangles === nothing
    @test uniform.packing_triangle_edge_ids === nothing
    @test uniform.faces !== nothing

    st_map = build_spanning_tree_map(; faces=40, seed=3)
    spanning_tree = prepare_layout_problem(st_map; dimension=3, options=Dict("engine" => "circle_packing"))
    @test spanning_tree.packing_boundary_vertices !== nothing
    @test length(spanning_tree.packing_boundary_vertices) >= 3
    @test length(unique(spanning_tree.packing_boundary_vertices)) == length(spanning_tree.packing_boundary_vertices)
    @test spanning_tree.packing_triangles !== nothing
    @test spanning_tree.packing_triangle_edge_ids !== nothing
    @test size(spanning_tree.packing_triangles) == size(spanning_tree.packing_triangle_edge_ids)
    @test 0 <= spanning_tree.metadata["circle_packing_removed_outer_vertex"] < num_vertices(st_map)
    @test spanning_tree.num_vertices <= num_vertices(st_map) - 1

    pos3, _, pack_meta = compute_circle_packing_layout(
        spanning_tree.num_vertices,
        spanning_tree.edges,
        spanning_tree.packing_boundary_vertices;
        triangles=spanning_tree.packing_triangles,
        triangle_edge_ids=spanning_tree.packing_triangle_edge_ids,
        project_to_sphere=true,
        return_metadata=true,
    )
    layout_meta = merge(copy(spanning_tree.metadata), pack_meta)
    render_pos, render_groups, render_geometry = DecoratedRandomPlanarMaps._augment_sphere_circle_render(
        pos3,
        spanning_tree.edge_groups,
        pack_meta["sphere_circle_geometry"],
        layout_meta,
    )
    outer = Int32(layout_meta["circle_packing_render_outer_vertex"])
    boundary = Int32.(layout_meta["circle_packing_outer_boundary_vertices"])
    star_keys = Set((min(outer, v), max(outer, v)) for v in boundary)

    finite_rows = count(i -> all(isfinite, render_pos[i, :]), 1:size(render_pos, 1))
    @test finite_rows == spanning_tree.num_vertices + 1
    @test render_geometry !== nothing
    @test size(render_geometry["normals"], 1) == size(render_pos, 1)
    generic_keys = Set(
        (min(render_groups["generic"][i, 1], render_groups["generic"][i, 2]),
         max(render_groups["generic"][i, 1], render_groups["generic"][i, 2])) for i in 1:size(render_groups["generic"], 1)
    )
    @test star_keys ⊆ generic_keys
    for (name, arr) in render_groups
        name == "generic" && continue
        keys = Set((min(arr[i, 1], arr[i, 2]), max(arr[i, 1], arr[i, 2])) for i in 1:size(arr, 1))
        @test isempty(intersect(star_keys, keys))
    end

    fk_map = build_hc_map(; faces=300, p=0.5, seed=7)
    fk_problem = prepare_layout_problem(fk_map; dimension=3, options=Dict("engine" => "circle_packing"))
    @test fk_problem.packing_boundary_vertices !== nothing
    @test length(fk_problem.packing_boundary_vertices) >= 3
    @test length(unique(fk_problem.packing_boundary_vertices)) == length(fk_problem.packing_boundary_vertices)
    @test fk_problem.packing_triangles !== nothing
    @test fk_problem.packing_triangle_edge_ids !== nothing
    @test size(fk_problem.packing_triangles) == size(fk_problem.packing_triangle_edge_ids)
    @test 0 <= fk_problem.metadata["circle_packing_removed_outer_vertex"] < num_vertices(fk_map)
end

@testset "FK triangle incidences stay manifold at the edge-id level" begin
    m = build_hc_map_from_word("hhcHHcCF")
    @test _max_triangle_edge_id_incidence(m.triangle_green_edges, m.triangle_diag_edge) <= 2

    prod_idx, start_step, end_step = DecoratedRandomPlanarMaps.maximal_gasket_span(m, "h")
    keep_steps = [s for s in Int(start_step):Int(end_step) if s != Int(m.order_face[Int(prod_idx) + 1])]
    gasket_green = Int32.(m.triangle_green_edges[keep_steps .+ 1, :])
    gasket_diag = Int32.(m.triangle_diag_edge[keep_steps .+ 1])
    @test _max_triangle_edge_id_incidence(gasket_green, gasket_diag) <= 2

    random_map = build_hc_map(; faces=80, p=0.25, seed=4)
    @test _max_triangle_edge_id_incidence(random_map.triangle_green_edges, random_map.triangle_diag_edge) <= 2
end

@testset "FK layout problems preserve edge multiplicity for Tutte and SFDP" begin
    m = build_hc_map_from_word("hhcHHcCF")

    problem2 = prepare_layout_problem(m; dimension=2, options=Dict("hc_boundary_mode" => "h_gasket"))
    @test problem2.metadata["layout_num_edges"] > problem2.metadata["layout_num_collapsed_edges"]
    pos2 = compute_tutte_layout(
        problem2.num_vertices,
        problem2.edges,
        problem2.boundary_vertices;
        boundary_positions=problem2.boundary_positions,
        solver="direct",
    )
    @test size(pos2) == (problem2.num_vertices, 2)
    @test all(isfinite, pos2)

    problem3 = prepare_layout_problem(m; dimension=3)
    @test problem3.metadata["layout_num_edges"] > problem3.metadata["layout_num_collapsed_edges"]
    pos3 = compute_sfdp_layout(problem3.num_vertices, problem3.edges; seed=5, iterations=30, scale=1.0)
    @test size(pos3) == (problem3.num_vertices, 3)
    @test all(isfinite, pos3)
end

@testset "Unsupported layout combinations are rejected explicitly" begin
    st = build_spanning_tree_map(; faces=20, seed=2)
    uniform = generate_uniform_map(; faces=20, seed=2)
    fk = build_hc_map(; faces=40, p=0.5, seed=2)

    err_st = @test_throws ArgumentError prepare_layout_problem(st; dimension=2)
    @test occursin("spanning_tree maps do not currently support 2D embeddings", sprint(showerror, err_st.value))

    err_uniform_2d = @test_throws ArgumentError prepare_layout_problem(uniform; dimension=2, options=Dict("engine" => "circle_packing"))
    @test occursin("uniform maps do not currently support circle_packing", sprint(showerror, err_uniform_2d.value))

    err_uniform_3d = @test_throws ArgumentError prepare_layout_problem(uniform; dimension=3, options=Dict("engine" => "circle_packing"))
    @test occursin("uniform maps do not currently support circle_packing", sprint(showerror, err_uniform_3d.value))

    err_fk_face = @test_throws ArgumentError prepare_layout_problem(fk; dimension=2, options=Dict("hc_boundary_mode" => "face"))
    @test occursin("FK 2D layouts only support", sprint(showerror, err_fk_face.value))
end

@testset "H gasket uses red boundary loop" begin
    m = build_hc_map_from_word("hhcHHcCF")
    problem = prepare_layout_problem(m; dimension=2, options=Dict("hc_boundary_mode" => "h_gasket"))

    original_vertices = Int.(problem.metadata["gasket_original_vertex_ids"])
    boundary = Int.(problem.boundary_vertices)
    boundary_old = original_vertices[boundary .+ 1]
    root = Int(problem.metadata["gasket_root_prod_index"])

    @test all(m.vertex_color[boundary_old .+ 1] .== 0)
    @test boundary_old[1] == Int(m.opposite_vertex_at_production[root + 1])
    @test !(Int(m.predecessor_vertex[root + 1]) in Set(original_vertices))
    @test Int(m.vertex_of_prod[root + 1]) in Set(original_vertices)
    @test problem.metadata["gasket_num_triangles"] == 3
    @test problem.metadata["gasket_num_collapsed_triangles"] == 3
    @test problem.metadata["gasket_num_edges"] == problem.metadata["gasket_num_collapsed_edges"]
    @test problem.metadata["gasket_removed_parallel_strip_vertices"] == 1
    @test problem.metadata["gasket_collapsed_residual_parallel_edges"] == 1

    red_edges = _edge_set(problem.edge_groups["red"])
    orange_edges = _edge_set(problem.edge_groups["orange"])
    for i in 1:(length(boundary) - 1)
        @test (min(boundary[i], boundary[i + 1]), max(boundary[i], boundary[i + 1])) in red_edges
    end
    @test (min(boundary[end], boundary[1]), max(boundary[end], boundary[1])) in orange_edges

    pos = compute_tutte_layout(
        problem.num_vertices,
        problem.edges,
        problem.boundary_vertices;
        boundary_positions=problem.boundary_positions,
        solver="direct",
    )
    @test size(pos) == (problem.num_vertices, 2)
end

@testset "C gasket uses blue boundary loop" begin
    m = build_hc_map_from_word("cchCChHF")
    problem = prepare_layout_problem(m; dimension=2, options=Dict("hc_boundary_mode" => "c_gasket"))

    original_vertices = Int.(problem.metadata["gasket_original_vertex_ids"])
    boundary = Int.(problem.boundary_vertices)
    boundary_old = original_vertices[boundary .+ 1]
    root = Int(problem.metadata["gasket_root_prod_index"])

    @test all(m.vertex_color[boundary_old .+ 1] .== 1)
    @test boundary_old[1] == Int(m.opposite_vertex_at_production[root + 1])
    @test !(Int(m.predecessor_vertex[root + 1]) in Set(original_vertices))
    @test Int(m.vertex_of_prod[root + 1]) in Set(original_vertices)

    blue_edges = _edge_set(problem.edge_groups["blue"])
    purple_edges = _edge_set(problem.edge_groups["purple"])
    for i in 1:(length(boundary) - 1)
        @test (min(boundary[i], boundary[i + 1]), max(boundary[i], boundary[i + 1])) in blue_edges
    end
    @test (min(boundary[end], boundary[1]), max(boundary[end], boundary[1])) in purple_edges
end

@testset "Degenerate FK gasket holes collapse in 2D" begin
    base = prepare_layout_problem(
        build_hc_map_from_word("cchCCF");
        dimension=2,
        options=Dict("hc_boundary_mode" => "c_gasket"),
    )
    reduced = prepare_layout_problem(
        build_hc_map_from_word("cchCCcFF");
        dimension=2,
        options=Dict("hc_boundary_mode" => "c_gasket"),
    )

    @test reduced.metadata["gasket_cleanup_mode"] == "fk_word_small_holes"
    @test haskey(reduced.metadata, "gasket_clean_fk_word")
    @test reduced.metadata["gasket_removed_size1_holes"] + reduced.metadata["gasket_removed_size2_holes"] > 0
    @test reduced.metadata["gasket_num_triangles"] == 3

    @test reduced.num_vertices == base.num_vertices
    @test reduced.boundary_vertices == base.boundary_vertices
    @test reduced.surface_triangles == base.surface_triangles
    @test reduced.edges == base.edges
    @test reduced.edge_groups == base.edge_groups
    @test all(reduced.edges[:, 1] .!= reduced.edges[:, 2])
end

@testset "HC gasket export" begin
    m = build_hc_map(; faces=80, p=0.25, seed=4)
    problem = prepare_layout_problem(m; dimension=2, options=Dict("hc_boundary_mode" => "h_gasket"))
    pos = compute_tutte_layout(
        problem.num_vertices,
        problem.edges,
        problem.boundary_vertices;
        boundary_positions=problem.boundary_positions,
        seed=7,
        solver="direct",
    )

    mktempdir() do out_dir
        export_web_binaries(
            problem.render_map_data,
            pos,
            out_dir;
            edge_groups=problem.edge_groups,
            faces=problem.faces,
            triangles=problem.surface_triangles,
            metadata=problem.metadata,
        )

        @test isfile(joinpath(out_dir, "data_verts.bin"))
        @test isfile(joinpath(out_dir, "web_meta.json"))
        @test isfile(joinpath(out_dir, "index.html"))
    end
end

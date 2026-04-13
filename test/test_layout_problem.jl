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

    problem2 = prepare_layout_problem(m; dimension=2, options=Dict("hc_boundary_mode" => "face"))
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
    @test problem.metadata["gasket_num_triangles"] == 5
    @test problem.metadata["gasket_num_collapsed_triangles"] == 4
    @test problem.metadata["gasket_num_edges"] > problem.metadata["gasket_num_collapsed_edges"]

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

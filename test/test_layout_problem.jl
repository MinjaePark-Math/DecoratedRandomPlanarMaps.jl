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

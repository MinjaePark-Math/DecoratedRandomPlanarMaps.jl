using Test
using Random
using DecoratedRandomPlanarMaps

@testset "Mated-CRT parameter resolution" begin
    params = DecoratedRandomPlanarMaps._resolve_mated_crt_parameters(; gamma=sqrt(2.0))
    @test params["kappa"] ≈ 2.0 atol=1.0e-12
    @test params["kappa_prime"] ≈ 8.0 atol=1.0e-12
    @test params["gamma_prime"] ≈ 2sqrt(2.0) atol=1.0e-12
    @test params["correlation"] ≈ 0.0 atol=1.0e-12

    endpoint = DecoratedRandomPlanarMaps._resolve_mated_crt_parameters(; gamma=2.0)
    @test endpoint["gamma"] ≈ 2.0 atol=1.0e-12
    @test endpoint["kappa"] ≈ 4.0 atol=1.0e-12
    @test endpoint["kappa_prime"] ≈ 4.0 atol=1.0e-12
    @test endpoint["gamma_prime"] ≈ 2.0 atol=1.0e-12
    @test endpoint["correlation"] ≈ 1.0 atol=1.0e-12

    endpoint_from_prime = DecoratedRandomPlanarMaps._resolve_mated_crt_parameters(; gamma_prime=2.0, kappa_prime=4.0)
    @test endpoint_from_prime["gamma"] ≈ 2.0 atol=1.0e-12

    cfg_map = build_map_from_config(Dict(
        "type" => "mated_crt",
        "vertices" => 12,
        "topology" => "disk",
        "kappa_prime" => 8.0,
        "seed" => 1,
        "refinement" => 3,
        "burnin_sweeps" => 8,
        "gibbs_sweeps" => 16,
    ))
    @test cfg_map isa MatedCRTMap
    @test cfg_map.topology == :disk
    @test cfg_map.gamma ≈ sqrt(2.0) atol=1.0e-12
    @test cfg_map.kappa_prime ≈ 8.0 atol=1.0e-12

    approx_cfg_map = build_map_from_config(Dict(
        "type" => "mated_crt",
        "vertices" => 18,
        "topology" => "sphere",
        "gamma" => 0.6,
        "seed" => 3,
        "sphere_sampler" => "approx",
        "sphere_exact_tail_cutoff" => 8,
    ))
    @test approx_cfg_map isa MatedCRTMap
    @test approx_cfg_map.topology == :sphere
    @test occursin("approx_hybrid", approx_cfg_map.sampler)
    @test approx_cfg_map.sphere_layout_map !== nothing
end

@testset "Mated-CRT render groups use red/blue without FK exploration" begin
    m = build_mated_crt_map(; vertices=12, topology="disk", gamma=sqrt(2.0), seed=1)
    groups = DecoratedRandomPlanarMaps.grouped_edges(m)

    @test Set(["generic", "red", "blue"]) ⊆ Set(keys(groups))
    @test !haskey(groups, "upper")
    @test !haskey(groups, "lower")
    @test DecoratedRandomPlanarMaps.edge_group_color("red") == "#c0392b"
    @test DecoratedRandomPlanarMaps.edge_group_color("blue") == "#2980b9"
    @test !DecoratedRandomPlanarMaps.should_include_exploration(
        m;
        edge_groups=groups,
        metadata=Dict("model" => "mated_crt"),
    )
end

@testset "Mated-CRT disk topology supports harmonic-boundary Tutte layout" begin
    m = build_mated_crt_map(
        ;
        vertices=12,
        topology="disk",
        gamma=sqrt(2.0),
        seed=1,
        refinement=3,
        burnin_sweeps=8,
        gibbs_sweeps=16,
    )
    problem = prepare_layout_problem(m; dimension=2)

    @test problem.num_vertices == num_vertices(m)
    @test problem.boundary_vertices == m.boundary_vertices
    @test problem.boundary_positions === nothing
    @test problem.surface_triangles !== nothing
    @test problem.surface_triangle_edge_ids !== nothing
    @test size(problem.surface_triangles) == size(problem.surface_triangle_edge_ids)
    @test size(problem.surface_triangles, 2) == 3
    @test length(Set((min(m.edge_u[i], m.edge_v[i]), max(m.edge_u[i], m.edge_v[i])) for i in eachindex(m.edge_u))) == length(m.edge_u)

    pos, meta = compute_tutte_layout(
        problem.num_vertices,
        problem.edges,
        problem.boundary_vertices;
        boundary_positions=problem.boundary_positions,
        solver="direct",
        return_metadata=true,
    )
    @test size(pos) == (problem.num_vertices, 2)
    @test all(isfinite, pos)
    @test meta["boundary_positions_mode"] == "harmonic_measure"
end

@testset "Mated-CRT disk outer face is recovered on larger boundary-heavy samples" begin
    m = build_mated_crt_map(
        ;
        vertices=120,
        topology="disk",
        gamma=sqrt(2.0),
        seed=1,
        refinement=4,
        burnin_sweeps=64,
        gibbs_sweeps=128,
    )
    @test Int(m.outer_face_index) >= 0
    @test DecoratedRandomPlanarMaps._same_cyclic_cycle(
        m.faces[Int(m.outer_face_index) + 1],
        m.boundary_vertices,
    )
    @test all(length(face) == 3 for (i, face) in enumerate(m.faces) if i != Int(m.outer_face_index) + 1)
    @test length(m.edge_u) == 3 * num_vertices(m) - 3 - length(m.boundary_vertices)
    @test num_faces(m) - 1 == 2 * num_vertices(m) - 2 - length(m.boundary_vertices)

    problem = prepare_layout_problem(m; dimension=2)
    @test problem.surface_triangles !== nothing
    @test size(problem.surface_triangles, 2) == 3
    @test size(problem.surface_triangles) == size(problem.surface_triangle_edge_ids)
end

@testset "Mated-CRT disk quotient recovers collapsed triangles from Brownian seed 71" begin
    m = build_mated_crt_map(
        ;
        vertices=200,
        topology="disk",
        gamma=sqrt(2.0),
        seed=71,
        refinement=4,
        burnin_sweeps=64,
        gibbs_sweeps=128,
    )
    @test Int(m.outer_face_index) >= 0
    @test all(length(face) == 3 for (i, face) in enumerate(m.faces) if i != Int(m.outer_face_index) + 1)
    @test length(m.edge_u) == 3 * num_vertices(m) - 3 - length(m.boundary_vertices)
    @test num_faces(m) - 1 == 2 * num_vertices(m) - 2 - length(m.boundary_vertices)
end

@testset "Mated-CRT sphere topology supports sphere circle packing" begin
    m = build_mated_crt_map(; vertices=12, topology="sphere", gamma=sqrt(2.0), seed=2)
    problem = prepare_layout_problem(m; dimension=3, options=Dict("engine" => "circle_packing"))
    raw_edges = sanitize_edge_array(hcat(m.edge_u, m.edge_v))
    collapsed_edges = first(collapse_undirected_edges(raw_edges; drop_loops=true))
    exact_word = DecoratedRandomPlanarMaps.sample_uniform_excursion_word(
        DecoratedRandomPlanarMaps.QuadrantBridgeSampler(num_vertices(m) - 2),
        MersenneTwister(2),
    )

    @test isempty(m.boundary_vertices)
    @test Int(m.outer_face_index) == -1
    @test all(length(face) == 3 for face in m.faces)
    @test m.sphere_layout_map !== nothing
    @test DecoratedRandomPlanarMaps.variant(m.sphere_layout_map) == "spanning_tree"
    @test length(m.sphere_layout_map.word) == 2 * (num_vertices(m) - 2)
    @test m.sphere_layout_map.word == exact_word
    @test occursin("exact_spanning_tree_word", m.sampler)
    @test size(raw_edges, 1) > size(collapsed_edges, 1)

    @test problem.num_vertices <= num_vertices(m)
    @test problem.packing_boundary_vertices !== nothing
    @test length(problem.packing_boundary_vertices) >= 3
    @test length(unique(problem.packing_boundary_vertices)) == length(problem.packing_boundary_vertices)
    @test problem.packing_triangles !== nothing
    @test problem.packing_triangle_edge_ids !== nothing
    @test size(problem.packing_triangles) == size(problem.packing_triangle_edge_ids)
    @test problem.metadata["sphere_construction"] == "exact_cone_walk_spanning_tree"
    @test problem.metadata["layout_graph_mode"] == "spanning_tree_raw"
    if haskey(problem.metadata, "circle_packing_removed_outer_vertex")
        @test 0 <= problem.metadata["circle_packing_removed_outer_vertex"] < num_vertices(m)
        @test problem.num_vertices <= num_vertices(m) - 1
    else
        @test haskey(problem.metadata, "circle_packing_outer_triangle_index")
        @test 0 <= problem.metadata["circle_packing_outer_triangle_index"] < num_faces(m)
    end
    @test length(m.edge_u) == 3 * num_vertices(m) - 6
    @test num_faces(m) == 2 * num_vertices(m) - 4

    pos3, radii, meta = compute_circle_packing_layout(
        problem.num_vertices,
        problem.edges,
        problem.packing_boundary_vertices;
        triangles=problem.packing_triangles,
        triangle_edge_ids=problem.packing_triangle_edge_ids,
        project_to_sphere=true,
        maxiter=100,
        return_metadata=true,
    )

    @test size(pos3) == (problem.num_vertices, 3)
    @test length(radii) == problem.num_vertices
    @test all(isfinite, pos3)
    @test all(radii .> 0.0)
    @test meta["packing_topology"] == "sphere"
    @test meta["sphere_projection"] == "inverse_stereographic"
end

@testset "Mated-CRT sphere SFDP uses the raw spanning-tree multigraph" begin
    m = build_mated_crt_map(
        ;
        vertices=80,
        topology="sphere",
        gamma=sqrt(2.0),
        seed=2,
        refinement=4,
        burnin_sweeps=64,
        gibbs_sweeps=128,
    )
    problem = prepare_layout_problem(m; dimension=3, options=Dict("engine" => "sfdp"))
    raw_edges = sanitize_edge_array(hcat(m.edge_u, m.edge_v))
    collapsed_edges = first(DecoratedRandomPlanarMaps.collapse_undirected_edges(raw_edges; drop_loops=true))

    @test problem.num_vertices == num_vertices(m)
    @test problem.render_vertex_indices === nothing
    @test problem.edges == raw_edges
    @test size(problem.edges, 1) > size(collapsed_edges, 1)
    @test problem.metadata["layout_graph_mode"] == "spanning_tree_raw"
    @test problem.metadata["sphere_construction"] == "exact_cone_walk_spanning_tree"
    @test m.sphere_layout_map !== nothing
    @test length(m.sphere_layout_map.production_steps) == num_vertices(m) - 2

    pos3 = compute_sfdp_layout(problem.num_vertices, problem.edges; seed=5, iterations=30, scale=1.0)
    @test size(pos3) == (problem.num_vertices, 3)
    @test all(isfinite, pos3)
end

@testset "Mated-CRT sphere supports approximate hybrid cone-walk sampling" begin
    m = build_mated_crt_map(
        ;
        vertices=36,
        topology="sphere",
        gamma=0.6,
        seed=7,
        sphere_sampler=:approx,
        sphere_exact_tail_cutoff=10,
    )
    problem = prepare_layout_problem(m; dimension=3, options=Dict("engine" => "sfdp"))

    @test m.sphere_layout_map !== nothing
    @test length(m.sphere_layout_map.word) == 2 * (num_vertices(m) - 2)
    @test occursin("approx_hybrid", m.sampler)
    @test problem.metadata["sphere_construction"] == "approx_hybrid_cone_walk_spanning_tree"
    @test problem.metadata["layout_graph_mode"] == "spanning_tree_raw"

    L, R = DecoratedRandomPlanarMaps._mated_crt_exact_word_coordinates(m.sphere_layout_map.word)
    @test L[end] == 0.0
    @test R[end] == 0.0
    @test minimum(L) >= 0.0
    @test minimum(R) >= 0.0
end

@testset "Mated-CRT sphere supports the gamma equals two endpoint" begin
    m = build_mated_crt_map(
        ;
        vertices=12,
        topology="sphere",
        gamma=2.0,
        seed=9,
        sphere_sampler=:approx,
    )
    problem = prepare_layout_problem(m; dimension=3, options=Dict("engine" => "sfdp"))

    @test m.gamma ≈ 2.0 atol=1.0e-12
    @test m.brownian_correlation ≈ 1.0 atol=1.0e-12
    @test occursin("approx_perfect_corr_diagonal_excursion", m.sampler)
    @test problem.metadata["sphere_construction"] == "approx_hybrid_cone_walk_spanning_tree"
    @test problem.metadata["layout_graph_mode"] == "spanning_tree_raw"

    word = m.sphere_layout_map.word
    @test length(word) == 2 * (num_vertices(m) - 2)
    for i in 1:2:length(word)
        @test word[i:min(i + 1, end)] in ("hc", "ch", "HC", "CH")
    end
end

@testset "Mated-CRT disk rejects the gamma equals two endpoint" begin
    @test_throws ArgumentError build_mated_crt_map(
        ;
        vertices=12,
        topology="disk",
        gamma=2.0,
        seed=1,
    )
end

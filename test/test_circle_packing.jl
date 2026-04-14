using Test
using LinearAlgebra
using DecoratedRandomPlanarMaps

function _reduce_fk_fragment_test(word::AbstractString)
    chars = collect(replace(String(word), r"\s+" => ""))
    prod_symbols = Char[]
    prev_stack = Int[]
    next_stack = Int[]
    prev_same = Int[]
    top_any = -1
    top_by_type = Dict('h' => -1, 'c' => -1)
    unmatched_orders = Char[]

    for ch in chars
        if ch == 'h' || ch == 'c'
            idx = length(prod_symbols)
            push!(prod_symbols, ch)
            push!(prev_stack, top_any)
            push!(next_stack, -1)
            if top_any != -1
                next_stack[top_any + 1] = idx
            end
            top_any = idx
            push!(prev_same, top_by_type[ch])
            top_by_type[ch] = idx
            continue
        end

        target_type, idx = if ch == 'H'
            ('h', top_by_type['h'])
        elseif ch == 'C'
            ('c', top_by_type['c'])
        elseif ch == 'F'
            top_any == -1 ? ('\0', -1) : (prod_symbols[top_any + 1], top_any)
        else
            throw(ArgumentError("unexpected FK symbol $(repr(ch))"))
        end

        if idx == -1
            push!(unmatched_orders, ch)
            continue
        end

        ps = prev_stack[idx + 1]
        ns = next_stack[idx + 1]
        if ps != -1
            next_stack[ps + 1] = ns
        end
        if ns != -1
            prev_stack[ns + 1] = ps
        else
            top_any = ps
        end
        top_by_type[target_type] = prev_same[idx + 1]
    end

    unmatched_burgers = Char[]
    idx = top_any
    while idx != -1 && prev_stack[idx + 1] != -1
        idx = prev_stack[idx + 1]
    end
    while idx != -1
        push!(unmatched_burgers, prod_symbols[idx + 1])
        idx = next_stack[idx + 1]
    end

    return String(vcat(unmatched_orders, unmatched_burgers))
end

function _count_small_fresh_subgaskets_test(word::AbstractString)
    chars = collect(replace(String(word), r"\s+" => ""))
    prod_symbols = Char[]
    prod_steps = Int[]
    prev_stack = Int[]
    next_stack = Int[]
    prev_same = Int[]
    top_any = -1
    top_by_type = Dict('h' => -1, 'c' => -1)
    spans = Tuple{Int,Int}[]

    for (idx1, ch) in enumerate(chars)
        pos = idx1 - 1
        if ch == 'h' || ch == 'c'
            idx = length(prod_symbols)
            push!(prod_symbols, ch)
            push!(prod_steps, pos)
            push!(prev_stack, top_any)
            push!(next_stack, -1)
            if top_any != -1
                next_stack[top_any + 1] = idx
            end
            top_any = idx
            push!(prev_same, top_by_type[ch])
            top_by_type[ch] = idx
        elseif ch == 'H' || ch == 'C'
            target_type = ch == 'H' ? 'h' : 'c'
            idx = top_by_type[target_type]
            idx == -1 && continue
            ps = prev_stack[idx + 1]
            ns = next_stack[idx + 1]
            if ps != -1
                next_stack[ps + 1] = ns
            end
            if ns != -1
                prev_stack[ns + 1] = ps
            else
                top_any = ps
            end
            top_by_type[target_type] = prev_same[idx + 1]
        elseif ch == 'F'
            top_any == -1 && continue
            idx = top_any
            push!(spans, (prod_steps[idx + 1], pos))
            target_type = prod_symbols[idx + 1]
            ps = prev_stack[idx + 1]
            ns = next_stack[idx + 1]
            if ps != -1
                next_stack[ps + 1] = ns
            end
            if ns != -1
                prev_stack[ns + 1] = ps
            else
                top_any = ps
            end
            top_by_type[target_type] = prev_same[idx + 1]
        end
    end

    isempty(spans) && return 0
    root_start = minimum(first.(spans))
    root_end = maximum(last.(spans))
    count = 0
    for (start_step, end_step) in spans
        start_step == root_start && end_step == root_end && continue
        if length(_reduce_fk_fragment_test(String(chars[start_step + 1:end_step + 1]))) <= 1
            count += 1
        end
    end
    return count
end

@testset "Circle packing on a wheel graph" begin
    edges = Int32[
        0 1;
        1 2;
        2 3;
        3 0;
        0 4;
        1 4;
        2 4;
        3 4;
    ]
    boundary = Int32[0, 1, 2, 3]
    triangles = Int32[
        0 1 4;
        1 2 4;
        2 3 4;
        3 0 4;
    ]

    pos, radii, meta = compute_circle_packing_layout(
        5,
        edges,
        boundary;
        triangles=triangles,
        maxiter=200,
        tol=1.0e-10,
        return_metadata=true,
    )

    @test size(pos) == (5, 2)
    @test length(radii) == 5
    @test all(isfinite, pos)
    @test all(radii .> 0.0)
    @test meta["packing_num_triangles"] == 4

    for v in boundary
        i = Int(v) + 1
        @test norm(pos[i, :]) + radii[i] ≈ 1.0 atol=5.0e-4
    end

    tangent_edges = [
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 0),
        (0, 4),
        (1, 4),
        (2, 4),
        (3, 4),
    ]
    for (u, v) in tangent_edges
        iu = u + 1
        iv = v + 1
        @test norm(pos[iu, :] - pos[iv, :]) ≈ radii[iu] + radii[iv] atol=5.0e-3
    end
end

@testset "Circle packing can lift a disk packing to the sphere" begin
    edges = Int32[
        0 1;
        1 2;
        2 3;
        3 0;
        0 4;
        1 4;
        2 4;
        3 4;
    ]
    boundary = Int32[0, 1, 2, 3]
    triangles = Int32[
        0 1 4;
        1 2 4;
        2 3 4;
        3 0 4;
    ]

    pos3, sphere_radii, meta = compute_circle_packing_layout(
        5,
        edges,
        boundary;
        triangles=triangles,
        project_to_sphere=true,
        sphere_radius=2.0,
        maxiter=200,
        tol=1.0e-10,
        return_metadata=true,
    )

    @test size(pos3) == (5, 3)
    @test length(sphere_radii) == 5
    @test all(isfinite, pos3)
    @test all(sphere_radii .> 0.0)
    @test meta["packing_topology"] == "sphere"
    @test meta["sphere_projection"] == "inverse_stereographic"
    @test meta["sphere_radius"] == 2.0

    norms = vec(sqrt.(sum(pos3 .^ 2; dims=2)))
    @test all(abs.(norms .- 2.0) .<= 1.0e-6)

    sphere_geom = meta["sphere_circle_geometry"]
    @test size(sphere_geom["centers"]) == (5, 3)
    @test size(sphere_geom["normals"]) == (5, 3)
    @test length(sphere_geom["radii"]) == 5
    @test sphere_geom["sphere_radius"] == 2.0
end

@testset "Sphere projection scale is preserved in metadata" begin
    edges = Int32[
        0 1;
        1 2;
        2 3;
        3 0;
        0 4;
        1 4;
        2 4;
        3 4;
    ]
    boundary = Int32[0, 1, 2, 3]
    triangles = Int32[
        0 1 4;
        1 2 4;
        2 3 4;
        3 0 4;
    ]

    pos_default, _, meta_default = compute_circle_packing_layout(
        5,
        edges,
        boundary;
        triangles=triangles,
        project_to_sphere=true,
        maxiter=120,
        tol=1.0e-6,
        return_metadata=true,
    )

    pos_scaled, _, meta_scaled = compute_circle_packing_layout(
        5,
        edges,
        boundary;
        triangles=triangles,
        project_to_sphere=true,
        sphere_projection_scale=1.5,
        maxiter=120,
        tol=1.0e-6,
        return_metadata=true,
    )

    @test meta_default["sphere_projection_scale"] == 1.0
    @test meta_scaled["sphere_projection_scale"] == 1.5
    @test meta_scaled["sphere_circle_geometry"]["projection_scale"] == meta_scaled["sphere_effective_projection_scale"]
    @test all(isfinite, pos_default)
    @test all(isfinite, pos_scaled)
end

@testset "Sphere barycenter balance records consistent summary metadata" begin
    edges = Int32[
        0 1;
        1 2;
        2 3;
        3 0;
        0 4;
        1 4;
        2 4;
        3 4;
    ]
    boundary = Int32[0, 1, 2, 3]
    triangles = Int32[
        0 1 4;
        1 2 4;
        2 3 4;
        3 0 4;
    ]

    _, _, meta_balanced = compute_circle_packing_layout(
        5,
        edges,
        boundary;
        triangles=triangles,
        project_to_sphere=true,
        maxiter=120,
        tol=1.0e-6,
        return_metadata=true,
    )

    entries = DecoratedRandomPlanarMaps._sphere_balance_entries(
        meta_balanced["sphere_circle_geometry"];
        include_outer=true,
        outer_geometry=get(meta_balanced["sphere_circle_geometry"], "outer_geometry", nothing),
        outer_projection_scale=float(meta_balanced["sphere_effective_projection_scale"]),
    )
    barycenter = DecoratedRandomPlanarMaps._sphere_balance_barycenter(entries)
    barycenter_norm = norm(barycenter)

    @test meta_balanced["sphere_balance_barycenter_norm_after"] <= meta_balanced["sphere_balance_barycenter_norm_before"] + 1.0e-12
    @test barycenter_norm ≈ meta_balanced["sphere_balance_barycenter_norm_after"] atol=1.0e-10
    @test barycenter_norm <= 1.0e-8

    balanced_summary = DecoratedRandomPlanarMaps._sphere_balance_summary(
        meta_balanced["sphere_circle_geometry"];
        include_outer=true,
        outer_geometry=get(meta_balanced["sphere_circle_geometry"], "outer_geometry", nothing),
        outer_projection_scale=float(meta_balanced["sphere_effective_projection_scale"]),
    )
    @test meta_balanced["sphere_balance_method"] == "affine_barycenter_newton"
    @test balanced_summary.objective ≈ meta_balanced["sphere_balance_objective_after"] atol=1.0e-10
    @test balanced_summary.max_fraction ≈ meta_balanced["sphere_balance_max_cap_fraction_after"] atol=1.0e-10
    @test balanced_summary.min_fraction ≈ meta_balanced["sphere_balance_min_cap_fraction_after"] atol=1.0e-10
    @test balanced_summary.outer_fraction ≈ meta_balanced["sphere_balance_outer_cap_fraction_after"] atol=1.0e-10

    expected_outer = DecoratedRandomPlanarMaps._sphere_circle_geometry_from_disk(
        reshape(-Float64.(meta_balanced["sphere_balance_translation"]), 1, 2),
        Float64[1.0];
        sphere_radius=float(meta_balanced["sphere_radius"]),
        projection_scale=float(meta_balanced["sphere_effective_projection_scale"]),
    )
    actual_outer = meta_balanced["sphere_outer_circle_geometry"]
    @test actual_outer["centers"] ≈ expected_outer["centers"] atol=1.0e-10
    @test actual_outer["normals"] ≈ expected_outer["normals"] atol=1.0e-10
    @test actual_outer["offsets"] ≈ expected_outer["offsets"] atol=1.0e-10
end

@testset "Sphere barycenter balance preserves bounded-disk side information" begin
    m = build_spanning_tree_map(; faces=300, seed=712)
    problem = prepare_layout_problem(
        m;
        dimension=3,
        options=Dict("engine" => "circle_packing"),
    )

    _, _, meta = compute_circle_packing_layout(
        problem.num_vertices,
        problem.edges,
        problem.packing_boundary_vertices;
        triangles=problem.packing_triangles,
        triangle_edge_ids=problem.packing_triangle_edge_ids,
        project_to_sphere=true,
        return_metadata=true,
    )

    geom = meta["sphere_circle_geometry"]
    normals = geom["normals"]
    offsets = geom["offsets"]
    @test size(normals, 1) == length(offsets)
    @test all(normals[i, 3] < offsets[i] - 1.0e-10 for i in 1:size(normals, 1))
end

@testset "Obsolete outer-vertex preview options are rejected" begin
    m = build_spanning_tree_map(; faces=20, seed=1)
    @test_throws ArgumentError prepare_layout_problem(
        m;
        dimension=3,
        options=Dict(
            "engine" => "circle_packing",
            "outer_vertex_selection" => "first_valid",
        ),
    )
    @test_throws ArgumentError prepare_layout_problem(
        m;
        dimension=3,
        options=Dict(
            "engine" => "circle_packing",
            "outer_vertex_preview_candidates" => 4,
        ),
    )
end

@testset "Circle packing fan triangulates polygonal faces" begin
    edges = Int32[
        0 1;
        1 2;
        2 3;
        3 0;
    ]
    boundary = Int32[0, 1, 2, 3]
    faces = [Int32[0, 1, 2, 3]]

    pos, radii, meta = compute_circle_packing_layout(
        4,
        edges,
        boundary;
        faces=faces,
        maxiter=120,
        tol=1.0e-6,
        return_metadata=true,
    )

    @test size(pos) == (4, 2)
    @test length(radii) == 4
    @test all(isfinite, pos)
    @test all(radii .> 0.0)
    @test meta["triangulation_source"] == "faces_fan"
    @test meta["packing_num_edges"] == 5

    mktempdir() do out_dir
        svg_path = joinpath(out_dir, "preview.svg")
        export_svg_preview(
            nothing,
            pos,
            svg_path;
            edge_groups=Dict("generic" => edges),
            faces=faces,
            circle_radii=radii,
            metadata=meta,
        )

        export_web_binaries(
            nothing,
            pos,
            out_dir;
            edge_groups=Dict("generic" => edges),
            faces=faces,
            circle_radii=radii,
            metadata=meta,
            write_index_html=false,
        )

        @test occursin("<circle", read(svg_path, String))
        @test isfile(joinpath(out_dir, "data_radii.bin"))
    end
end

@testset "Circle packing on generated Schnyder map" begin
    m = generate_schnyder_map(; size_faces=8, seed=1)
    problem = prepare_layout_problem(m; dimension=2)

    pos, radii, meta = compute_circle_packing_layout(
        problem.num_vertices,
        problem.edges,
        problem.boundary_vertices;
        faces=problem.faces,
        maxiter=120,
        tol=1.0e-6,
        return_metadata=true,
    )

    @test size(pos) == (problem.num_vertices, 2)
    @test length(radii) == problem.num_vertices
    @test all(isfinite, pos)
    @test all(radii .> 0.0)
    @test meta["triangulation_source"] == "faces_fan"
end

@testset "Uniform circle packing is rejected explicitly" begin
    m = generate_uniform_map(; faces=100, seed=2)

    err2 = @test_throws ArgumentError prepare_layout_problem(m; dimension=2, options=Dict("engine" => "circle_packing"))
    @test occursin("uniform maps do not currently support circle_packing", sprint(showerror, err2.value))

    err3 = @test_throws ArgumentError prepare_layout_problem(m; dimension=3, options=Dict("engine" => "circle_packing"))
    @test occursin("uniform maps do not currently support circle_packing", sprint(showerror, err3.value))
end

@testset "Circle packing pipeline smoke test" begin
    cfg = Dict{String,Any}(
        "model" => Dict{String,Any}("type" => "schnyder", "faces" => 8, "seed" => 1),
        "layout" => Dict{String,Any}("dimension" => 2, "engine" => "circle_packing", "maxiter" => 50),
        "output" => Dict{String,Any}("show" => false),
    )

    timings = run_pipeline(cfg)
    @test timings isa TimingRecorder
end

@testset "Sphere circle packing pipeline smoke test" begin
    cfg = Dict{String,Any}(
        "model" => Dict{String,Any}("type" => "schnyder", "faces" => 8, "seed" => 1),
        "layout" => Dict{String,Any}("dimension" => 3, "engine" => "circle_packing", "maxiter" => 50),
        "output" => Dict{String,Any}("show" => false),
    )

    timings = run_pipeline(cfg)
    @test timings isa TimingRecorder
end

@testset "Schnyder sphere render records the appended outer circle index" begin
    m = generate_schnyder_map(; size_faces=8, seed=1)
    problem = prepare_layout_problem(m; dimension=3)
    pos3, _, meta = compute_circle_packing_layout(
        problem.num_vertices,
        problem.edges,
        problem.packing_boundary_vertices;
        faces=problem.packing_faces,
        triangles=problem.packing_triangles,
        triangle_edge_ids=problem.packing_triangle_edge_ids,
        project_to_sphere=true,
        maxiter=120,
        tol=1.0e-6,
        return_metadata=true,
    )
    layout_meta = merge(copy(problem.metadata), meta)
    render_pos, _, render_geometry = DecoratedRandomPlanarMaps._augment_sphere_circle_render(
        pos3,
        problem.edge_groups,
        meta["sphere_circle_geometry"],
        layout_meta,
        problem.render_map_data,
    )

    outer_vertex = Int(layout_meta["circle_packing_render_outer_vertex"])
    @test layout_meta["sphere_circle_geometry_includes_outer"] == true
    @test size(render_pos, 1) == problem.num_vertices + 1
    @test outer_vertex == problem.num_vertices
    @test render_pos[outer_vertex + 1, :] ≈ (-float(layout_meta["sphere_radius"]) .* render_geometry["normals"][outer_vertex + 1, :]) atol=1.0e-8
end

@testset "Circle packing preserves gasket triangle multiplicity" begin
    m = build_hc_map_from_word("hhcHHcCF")
    problem = prepare_layout_problem(m; dimension=2, options=Dict("hc_boundary_mode" => "h_gasket"))

    pos, radii, meta = compute_circle_packing_layout(
        problem.num_vertices,
        problem.edges,
        problem.boundary_vertices;
        triangles=problem.surface_triangles,
        triangle_edge_ids=problem.surface_triangle_edge_ids,
        maxiter=120,
        tol=1.0e-6,
        return_metadata=true,
    )

    @test size(pos) == (problem.num_vertices, 2)
    @test length(radii) == problem.num_vertices
    @test all(isfinite, pos)
    @test all(radii .> 0.0)
    @test meta["triangle_multiplicity_preserved"] == true
    @test meta["packing_num_triangles"] == 3
    @test problem.metadata["gasket_num_triangles"] == 3
    @test problem.metadata["gasket_removed_parallel_strip_vertices"] == 1

    for v in problem.boundary_vertices
        i = Int(v) + 1
        @test norm(pos[i, :]) + radii[i] ≈ 1.0 atol=5.0e-3
    end
end

@testset "Circle packing reports disconnected local sectors clearly" begin
    edges = Int32[
        0 1;
        1 2;
        2 0;
        0 3;
        1 3;
        2 3;
        3 4;
        4 5;
        5 3;
        3 4;
        4 5;
        5 3;
    ]
    boundary = Int32[0, 1, 2]
    triangles = Int32[
        3 0 1;
        3 1 2;
        3 2 0;
        3 4 5;
        3 5 4;
    ]

    err = @test_throws ArgumentError compute_circle_packing_layout(
        6,
        edges,
        boundary;
        triangles=triangles,
        maxiter=20,
        tol=1.0e-6,
    )
    @test occursin("disconnected sectors", sprint(showerror, err.value))
end

@testset "Circle packing succeeds after degenerate gasket cleanup" begin
    m = build_hc_map_from_word("cchCCcFF")
    problem = prepare_layout_problem(m; dimension=2, options=Dict("hc_boundary_mode" => "c_gasket"))

    @test problem.metadata["gasket_cleanup_mode"] == "fk_word_small_holes"
    @test haskey(problem.metadata, "gasket_clean_fk_word")

    pos, radii, meta = compute_circle_packing_layout(
        problem.num_vertices,
        problem.edges,
        problem.boundary_vertices;
        triangles=problem.surface_triangles,
        triangle_edge_ids=problem.surface_triangle_edge_ids,
        maxiter=120,
        tol=1.0e-6,
        return_metadata=true,
    )

    @test size(pos) == (problem.num_vertices, 2)
    @test length(radii) == problem.num_vertices
    @test all(isfinite, pos)
    @test all(radii .> 0.0)
    @test meta["packing_num_triangles"] == 3
end

@testset "Circle packing uses cleaned FK gasket word for 2D FK sample" begin
    m = build_hc_map(; faces=100, p=0.5, seed=2)
    problem = prepare_layout_problem(m; dimension=2, options=Dict("hc_boundary_mode" => "h_gasket"))

    @test problem.metadata["gasket_cleanup_mode"] == "fk_word_small_holes"
    @test haskey(problem.metadata, "gasket_clean_fk_word")
    @test _reduce_fk_fragment_test(problem.metadata["gasket_clean_fk_word"]) == "HH"
    @test problem.metadata["gasket_clean_reduced_word"] == "HH"
    @test _count_small_fresh_subgaskets_test(problem.metadata["gasket_clean_fk_word"]) == 0

    pos, radii, meta = compute_circle_packing_layout(
        problem.num_vertices,
        problem.edges,
        problem.boundary_vertices;
        triangles=problem.surface_triangles,
        triangle_edge_ids=problem.surface_triangle_edge_ids,
        maxiter=150,
        tol=1.0e-6,
        return_metadata=true,
    )

    @test size(pos) == (problem.num_vertices, 2)
    @test length(radii) == problem.num_vertices
    @test all(isfinite, pos)
    @test all(radii .> 0.0)
    @test meta["packing_num_triangles"] == problem.metadata["gasket_num_triangles"]
end

@testset "Spanning-tree 2D layouts are rejected explicitly" begin
    m = build_spanning_tree_map(; faces=6, seed=7)

    err_tutte = @test_throws ArgumentError prepare_layout_problem(m; dimension=2)
    @test occursin("spanning_tree maps do not currently support 2D embeddings", sprint(showerror, err_tutte.value))

    err_cp = @test_throws ArgumentError prepare_layout_problem(m; dimension=2, options=Dict("engine" => "circle_packing"))
    @test occursin("spanning_tree maps do not currently support 2D embeddings", sprint(showerror, err_cp.value))
end

@testset "Circle packing collapses parallel FK strips before packing" begin
    m = build_spanning_tree_map(; faces=6, seed=7)

    problem3 = prepare_layout_problem(m; dimension=3, options=Dict("engine" => "circle_packing"))
    @test problem3.num_vertices == 5
    @test size(get(problem3.edge_groups, "blue", Matrix{Int32}(undef, 0, 2)), 1) == 1
    @test get(problem3.metadata, "gasket_removed_parallel_strip_vertices", 0) == 2
    @test get(problem3.metadata, "gasket_collapsed_residual_parallel_edges", 0) == 1
    @test size(problem3.packing_triangles, 1) == 5
    @test size(problem3.edges, 1) == length(Set((min(problem3.edges[i, 1], problem3.edges[i, 2]), max(problem3.edges[i, 1], problem3.edges[i, 2])) for i in 1:size(problem3.edges, 1)))
end

@testset "Sphere FK cleanup prunes nontriangulated leftovers" begin
    m = build_hc_map(; faces=50, p=0.5, seed=2)
    problem = prepare_layout_problem(m; dimension=3, options=Dict("engine" => "circle_packing"))

    @test problem.num_vertices >= 4
    @test get(problem.metadata, "gasket_removed_nontriangulated_vertices", 0) > 0
    @test get(problem.metadata, "gasket_removed_nontriangulated_edges", 0) >= 0
    @test sort!(unique(vec(Int32.(problem.packing_triangles)))) == Int32.(collect(0:(problem.num_vertices - 1)))
    @test all(v -> any(problem.edges[i, 1] == v || problem.edges[i, 2] == v for i in 1:size(problem.edges, 1)), Int32.(0:(problem.num_vertices - 1)))

    pos3, _, meta = compute_circle_packing_layout(
        problem.num_vertices,
        problem.edges,
        problem.packing_boundary_vertices;
        triangles=problem.packing_triangles,
        triangle_edge_ids=problem.packing_triangle_edge_ids,
        project_to_sphere=true,
        return_metadata=true,
    )
    layout_meta = merge(copy(problem.metadata), meta)
    render_pos, _, render_geometry = DecoratedRandomPlanarMaps._augment_sphere_circle_render(
        pos3,
        problem.edge_groups,
        meta["sphere_circle_geometry"],
        layout_meta,
        problem.render_map_data,
    )

    finite_circle_rows = count(i -> isfinite(render_geometry["radii"][i]), eachindex(render_geometry["radii"]))
    finite_pos_rows = count(i -> all(isfinite, render_pos[i, :]), 1:size(render_pos, 1))
    @test finite_circle_rows == problem.num_vertices + 1
    @test finite_pos_rows == problem.num_vertices + 1
end

@testset "Web viewer includes circle-packing overlays and defaults" begin
    viewer_source = read(joinpath(@__DIR__, "..", "assets", "index.html"), String)
    @test occursin("circle_radii", viewer_source)
    @test occursin("addLegendToggle('circles', sphereCircleOverlay", viewer_source)
    @test occursin("if (faceFile && !sphereProjection)", viewer_source)
    @test occursin("normalizePolylineToSphere", viewer_source)
    @test occursin("buildSphereCapSurfaceGeometry", viewer_source)
    @test occursin("buildSphereCircleSegmentsFromCaps", viewer_source)
    @test occursin("buildGeodesicEdgeArcPositions", viewer_source)
    @test occursin("buildIndexedEdgeArcPositions", viewer_source)
    @test occursin("buildPositionSegmentArcPositions", viewer_source)
end

@testset "Web export writes spherical circle overlays" begin
    edges = Int32[
        0 1;
        1 2;
        2 3;
        3 0;
        0 4;
        1 4;
        2 4;
        3 4;
    ]
    boundary = Int32[0, 1, 2, 3]
    triangles = Int32[
        0 1 4;
        1 2 4;
        2 3 4;
        3 0 4;
    ]

    pos3, _, meta = compute_circle_packing_layout(
        5,
        edges,
        boundary;
        triangles=triangles,
        project_to_sphere=true,
        maxiter=120,
        tol=1.0e-6,
        return_metadata=true,
    )

    mktempdir() do out_dir
        export_web_binaries(
            nothing,
            pos3,
            out_dir;
            edge_groups=Dict("generic" => edges),
            triangles=triangles,
            sphere_circle_geometry=meta["sphere_circle_geometry"],
            metadata=Dict("model" => "generic", "dimension" => 3),
            write_index_html=false,
        )

        meta_text = read(joinpath(out_dir, "web_meta.json"), String)
        @test isfile(joinpath(out_dir, "data_sphere_circles.bin"))
        @test isfile(joinpath(out_dir, "data_sphere_circle_caps.bin"))
        @test occursin("data_sphere_circles.bin", meta_text)
        @test occursin("data_sphere_circle_caps.bin", meta_text)
        @test occursin("\"name\":\"circles\"", meta_text)
    end
end

@testset "Sphere web export keeps FK-like exploration data" begin
    m = build_spanning_tree_map(; faces=20, seed=2)
problem = prepare_layout_problem(m; dimension=3, options=Dict("engine" => "circle_packing"))
    pos3, _, meta = compute_circle_packing_layout(
        problem.num_vertices,
        problem.edges,
        problem.packing_boundary_vertices;
        triangles=problem.packing_triangles,
        triangle_edge_ids=problem.packing_triangle_edge_ids,
        project_to_sphere=true,
        sphere_projection_scale=1.4,
        maxiter=120,
        tol=1.0e-6,
        return_metadata=true,
    )
    layout_meta = merge(copy(problem.metadata), meta)
    render_pos, render_groups, render_geometry = DecoratedRandomPlanarMaps._augment_sphere_circle_render(
        pos3,
        problem.edge_groups,
        meta["sphere_circle_geometry"],
        layout_meta,
        problem.render_map_data,
    )

    expected_groups = DecoratedRandomPlanarMaps._remap_edge_groups_for_render(
        problem.edge_groups,
        layout_meta["circle_packing_packing_to_render_vertices"],
    )
    for (name, arr) in DecoratedRandomPlanarMaps._recover_outer_vertex_edge_groups(
        problem.render_map_data,
        Int32(layout_meta["circle_packing_removed_outer_vertex"]),
    )
        expected_groups = DecoratedRandomPlanarMaps._append_edge_group_edges(expected_groups, name, arr)
    end
    @test render_groups == expected_groups

    expected_triangles = DecoratedRandomPlanarMaps._augment_sphere_render_triangles(
        problem.surface_triangles,
        size(render_pos, 1),
        layout_meta,
    )
    @test maximum(expected_triangles) < size(render_pos, 1)
    outer_vertex = Int(layout_meta["circle_packing_render_outer_vertex"])
    boundary = Int.(layout_meta["circle_packing_outer_boundary_vertices"])
    expected_fan = Set(
        Set([outer_vertex, boundary[i], boundary[mod1(i + 1, length(boundary))]])
        for i in eachindex(boundary)
    )
    actual_fan = Set(
        Set(Int.(expected_triangles[i, :]))
        for i in 1:size(expected_triangles, 1)
        if outer_vertex in expected_triangles[i, :]
    )
    @test expected_fan ⊆ actual_fan
    @test render_pos[outer_vertex + 1, :] ≈ (-float(layout_meta["sphere_radius"]) .* render_geometry["normals"][outer_vertex + 1, :]) atol=1.0e-8

    mktempdir() do out_dir
        export_web_binaries(
            problem.render_map_data,
            render_pos,
            out_dir;
            edge_groups=render_groups,
            triangles=problem.surface_triangles,
            sphere_circle_geometry=render_geometry,
            metadata=merge(Dict("model" => "spanning_tree", "dimension" => 3), layout_meta),
            write_index_html=false,
        )

        meta_text = read(joinpath(out_dir, "web_meta.json"), String)
        @test occursin("data_exploration.bin", meta_text)
        @test isfile(joinpath(out_dir, "data_exploration.bin"))
        @test occursin("sphere_circle_geometry_includes_outer", meta_text)
    end
end

@testset "Restored outer sphere circle uses the balanced normalized frame" begin
    m = build_spanning_tree_map(; faces=20, seed=2)
problem = prepare_layout_problem(m; dimension=3, options=Dict("engine" => "circle_packing"))
    pos3, _, meta = compute_circle_packing_layout(
        problem.num_vertices,
        problem.edges,
        problem.packing_boundary_vertices;
        triangles=problem.packing_triangles,
        triangle_edge_ids=problem.packing_triangle_edge_ids,
        project_to_sphere=true,
        sphere_projection_scale=1.4,
        maxiter=120,
        tol=1.0e-6,
        return_metadata=true,
    )
    layout_meta = merge(copy(problem.metadata), meta)
    _, _, render_geometry = DecoratedRandomPlanarMaps._augment_sphere_circle_render(
        pos3,
        problem.edge_groups,
        meta["sphere_circle_geometry"],
        layout_meta,
        problem.render_map_data,
    )

    outer_vertex = Int(layout_meta["circle_packing_render_outer_vertex"])
    expected_geometry = layout_meta["sphere_outer_circle_geometry"]

    @test render_geometry["centers"][outer_vertex + 1, :] ≈ expected_geometry["centers"][1, :] atol=1.0e-8
    @test render_geometry["normals"][outer_vertex + 1, :] ≈ expected_geometry["normals"][1, :] atol=1.0e-8
    @test render_geometry["radii"][outer_vertex + 1] ≈ expected_geometry["radii"][1] atol=1.0e-8
    @test render_geometry["offsets"][outer_vertex + 1] ≈ expected_geometry["offsets"][1] atol=1.0e-8
    render_pos, _, render_geometry = DecoratedRandomPlanarMaps._augment_sphere_circle_render(
        pos3,
        problem.edge_groups,
        meta["sphere_circle_geometry"],
        layout_meta,
        problem.render_map_data,
    )
    @test render_pos[outer_vertex + 1, :] ≈ (-float(layout_meta["sphere_radius"]) .* render_geometry["normals"][outer_vertex + 1, :]) atol=1.0e-8
end

@testset "Sphere render remaps local packed edges back to render-space" begin
    m = build_spanning_tree_map(; faces=6, seed=7)
problem = prepare_layout_problem(m; dimension=3, options=Dict("engine" => "circle_packing"))
    pos3, _, meta = compute_circle_packing_layout(
        problem.num_vertices,
        problem.edges,
        problem.packing_boundary_vertices;
        triangles=problem.packing_triangles,
        triangle_edge_ids=problem.packing_triangle_edge_ids,
        project_to_sphere=true,
        return_metadata=true,
    )
    layout_meta = merge(copy(problem.metadata), meta)
    render_pos, render_groups, _ = DecoratedRandomPlanarMaps._augment_sphere_circle_render(
        pos3,
        problem.edge_groups,
        meta["sphere_circle_geometry"],
        layout_meta,
        problem.render_map_data,
    )

    undirected(arr) = Set((min(row[1], row[2]), max(row[1], row[2])) for row in eachrow(arr))

    @test undirected(render_groups["blue"]) == Set((
        (Int32(1), Int32(5)),
    ))
    @test undirected(render_groups["red"]) == Set((
        (Int32(0), Int32(2)),
        (Int32(2), Int32(3)),
        (Int32(3), Int32(4)),
    ))
    @test undirected(render_groups["generic"]) == Set((
        (Int32(0), Int32(1)),
        (Int32(0), Int32(5)),
        (Int32(1), Int32(2)),
        (Int32(1), Int32(3)),
        (Int32(1), Int32(4)),
        (Int32(2), Int32(5)),
        (Int32(3), Int32(5)),
        (Int32(4), Int32(5)),
    ))
    expected_render_geometry = DecoratedRandomPlanarMaps._expand_sphere_circle_geometry_for_render(
        meta["sphere_circle_geometry"],
        Int(layout_meta["circle_packing_render_vertex_count"]),
        layout_meta["circle_packing_packing_to_render_vertices"],
        Int(layout_meta["circle_packing_render_outer_vertex"]),
        layout_meta,
    )
    @test render_pos[Int(layout_meta["circle_packing_render_outer_vertex"]) + 1, :] ≈ (-float(layout_meta["sphere_radius"]) .* expected_render_geometry["normals"][Int(layout_meta["circle_packing_render_outer_vertex"]) + 1, :]) atol=1.0e-8

    render_triangles = DecoratedRandomPlanarMaps._augment_sphere_render_triangles(
        problem.surface_triangles,
        8,
        layout_meta,
    )
    @test Set(Tuple.(eachrow(render_triangles))) == Set((
        (Int32(2), Int32(1), Int32(3)),
        (Int32(3), Int32(1), Int32(4)),
        (Int32(4), Int32(1), Int32(5)),
        (Int32(5), Int32(3), Int32(4)),
        (Int32(5), Int32(2), Int32(3)),
        (Int32(0), Int32(1), Int32(5)),
        (Int32(0), Int32(5), Int32(2)),
        (Int32(0), Int32(2), Int32(1)),
    ))

    exploration = DecoratedRandomPlanarMaps.exploration_segment_points(
        m,
        render_pos;
        edge_groups=render_groups,
        triangles=render_triangles,
        metadata=layout_meta,
    )
    @test size(exploration, 1) == 16
    @test exploration[1, :] ≈ exploration[end, :]
end

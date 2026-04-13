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

@testset "Circle packing pipeline smoke test" begin
    cfg = Dict{String,Any}(
        "model" => Dict{String,Any}("type" => "schnyder", "faces" => 8, "seed" => 1),
        "layout" => Dict{String,Any}("dimension" => 2, "engine" => "circle_packing", "maxiter" => 50),
        "output" => Dict{String,Any}("show" => false),
    )

    timings = run_pipeline(cfg)
    @test timings isa TimingRecorder
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
    @test meta["packing_num_triangles"] == 5
    @test problem.metadata["gasket_num_triangles"] == 5

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

@testset "Circle packing on spanning-tree face boundary preserves triangle incidences" begin
    m = build_spanning_tree_map(; faces=100, seed=2)
    problem = prepare_layout_problem(m; dimension=2)

    @test problem.metadata["boundary_mode"] == "face"
    @test problem.surface_triangles !== nothing
    @test problem.surface_triangle_edge_ids !== nothing
    @test size(problem.surface_triangles) == size(problem.surface_triangle_edge_ids)

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
    @test meta["neighbor_order_source"] == "triangle_edge_ids"
    @test meta["triangle_multiplicity_preserved"] == true
    @test meta["packing_num_triangles"] == size(problem.surface_triangles, 1)
end

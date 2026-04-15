using Test
using JSON3
using Random
using DecoratedRandomPlanarMaps

@testset "Half-plane meandric layout preparation" begin
    m = build_half_plane_meandric_system(; order=6, boundary_half_length=3, seed=2)
    groups = grouped_edges(m)
    problem = prepare_layout_problem(m; dimension=2)

    @test num_vertices(m) == 12
    @test length(m.boundary_vertices) == 6
    @test Set(["generic", "upper", "lower"]) ⊆ Set(keys(groups))
    @test problem.boundary_vertices !== nothing
    @test problem.boundary_positions === nothing
    @test problem.metadata["model"] == "half_plane_meandric"

    pos, meta = compute_tutte_layout(
        problem.num_vertices,
        problem.edges,
        problem.boundary_vertices;
        boundary_positions=problem.boundary_positions,
        solver="direct",
        seed=7,
        return_metadata=true,
    )

    @test size(pos) == (problem.num_vertices, 2)
    @test all(isfinite, pos)
    @test meta["boundary_positions_mode"] == "harmonic_measure"
end

@testset "Uniform meandric uses auxiliary SFDP graph" begin
    m = build_uniform_meandric_system(; order=6, seed=3)
    groups = grouped_edges(m)
    problem = prepare_layout_problem(m; dimension=3)
    base_aux_edges = DecoratedRandomPlanarMaps.auxiliary_layout_edges(m; drop_loops=true)

    @test num_vertices(m) == 12
    @test any(startswith(String(name), "zz_loop_") for name in keys(groups))
    @test problem.num_vertices > num_vertices(m)
    @test problem.render_vertex_indices !== nothing
    @test length(problem.render_vertex_indices) == num_vertices(m)
    @test problem.metadata["model"] == "uniform_meandric"
    @test problem.metadata["layout_graph"] == "glued_trees_augmented"
    @test problem.metadata["layout_graph_base"] == "glued_trees"
    @test Int(problem.metadata["layout_base_edge_count"]) == size(base_aux_edges, 1)
    @test Int(problem.metadata["layout_augmented_edge_count"]) == size(problem.edges, 1)
    @test size(problem.edges, 1) > size(base_aux_edges, 1)
    @test Int(problem.metadata["layout_temporary_triangulation_edge_count"]) > 0

    pos = compute_sfdp_layout(problem.num_vertices, problem.edges; seed=5, iterations=40, scale=1.0)
    render_pos = pos[Int.(problem.render_vertex_indices) .+ 1, :]

    @test size(pos) == (problem.num_vertices, 3)
    @test size(render_pos) == (num_vertices(m), 3)
    @test all(isfinite, pos)
    @test all(isfinite, render_pos)
end

@testset "Closed meandric faces alternate loop and strand edges" begin
    m = build_uniform_meandric_system(; order=6, seed=3)

    @test length(m.faces) == 2 * (Int(m.order) + 1)
    @test length(m.face_edge_ids) == length(m.faces)
    @test length(m.face_edge_group) == length(m.edge_group) + 1

    generic_incidence = zeros(Int, length(m.upper_adj) - 1)
    for (face, edge_cycle) in zip(m.faces, m.face_edge_ids)
        @test length(face) == length(edge_cycle)
        @test !isempty(face)
        for i in eachindex(edge_cycle)
            edge_id = Int(edge_cycle[i]) + 1
            group = m.face_edge_group[edge_id]
            next_group = m.face_edge_group[Int(edge_cycle[mod1(i + 1, length(edge_cycle))]) + 1]
            @test group == "generic" ? next_group != "generic" : next_group == "generic"
            if group == "generic" && Int(edge_cycle[i]) < length(generic_incidence)
                generic_incidence[Int(edge_cycle[i]) + 1] += 1
            end
        end
    end

    @test all(generic_incidence .== 2)
end

@testset "Meandric loop highlights follow descending component size" begin
    half_plane = build_half_plane_meandric_system(; order=30, boundary_half_length=4, seed=1)
    half_groups = grouped_edges(half_plane)
    half_sizes = [size(half_groups[name], 1) for name in sort(filter(name -> startswith(String(name), "zz_loop_"), collect(keys(half_groups))))]

    closed = build_uniform_meandric_system(; order=6, seed=3)
    closed_groups = grouped_edges(closed)
    closed_sizes = [size(closed_groups[name], 1) for name in sort(filter(name -> startswith(String(name), "zz_loop_"), collect(keys(closed_groups))))]

    @test !isempty(half_sizes)
    @test issorted(half_sizes; rev=true)
    @test !isempty(closed_sizes)
    @test issorted(closed_sizes; rev=true)
end

@testset "Fast meander local rewires match full recomputation" begin
    rng = MersenneTwister(1234)

    for order in 2:8
        for _ in 1:30
            upper_steps, upper_adj = DecoratedRandomPlanarMaps.sample_uniform_dyck_path_with_adjacency(order, rng)
            lower_steps, lower_adj = DecoratedRandomPlanarMaps.sample_uniform_dyck_path_with_adjacency(order, rng)
            current_loops = length(DecoratedRandomPlanarMaps._closed_meandric_component_sizes(upper_adj, lower_adj))

            for update_upper in (true, false)
                base_steps = update_upper ? upper_steps : lower_steps
                base_adj = update_upper ? upper_adj : lower_adj

                for idx in 1:(2 * order - 1)
                    effect = DecoratedRandomPlanarMaps._dyck_adjacent_swap_effect(base_steps, base_adj, idx)

                    slow_steps = copy(base_steps)
                    slow_steps[idx], slow_steps[idx + 1] = slow_steps[idx + 1], slow_steps[idx]
                    slow_valid = true
                    slow_adj = Int32[]
                    try
                        slow_adj = DecoratedRandomPlanarMaps._pair_bracket_steps(slow_steps; allow_unmatched=false)
                    catch err
                        err isa ArgumentError || rethrow(err)
                        slow_valid = false
                    end

                    is_nontrivial_swap = base_steps[idx] != base_steps[idx + 1]
                    @test (effect !== nothing) == (slow_valid && is_nontrivial_swap)
                    effect === nothing && continue

                    fast_steps = copy(base_steps)
                    fast_adj = copy(base_adj)
                    fast_pool = DecoratedRandomPlanarMaps._build_dyck_swap_index_pool(fast_steps)
                    DecoratedRandomPlanarMaps._apply_dyck_adjacent_swap!(fast_steps, fast_adj, fast_pool, idx, effect)

                    @test fast_steps == slow_steps
                    @test fast_adj == slow_adj

                    same_component = DecoratedRandomPlanarMaps._meandric_vertices_share_component(
                        upper_adj,
                        lower_adj,
                        effect.test_u,
                        effect.test_v,
                    )
                    delta_loops = same_component ? 1 : -1
                    proposal_loops = update_upper ?
                        length(DecoratedRandomPlanarMaps._closed_meandric_component_sizes(fast_adj, lower_adj)) :
                        length(DecoratedRandomPlanarMaps._closed_meandric_component_sizes(upper_adj, fast_adj))

                    @test current_loops + delta_loops == proposal_loops
                end
            end
        end
    end
end

@testset "Uniform meander sampler returns a single loop" begin
    m = build_uniform_meander(
        ;
        order=6,
        seed=9,
        sweeps_per_temperature=3,
        search_sweeps=12,
        mixing_sweeps=3,
        restarts=12,
    )
    problem = prepare_layout_problem(m; dimension=3)

    @test num_vertices(m) == 12
    @test length(m.component_sizes) == 1
    @test sum(m.component_sizes) == num_vertices(m)
    @test problem.metadata["model"] == "uniform_meander"
    @test problem.metadata["num_loops"] == 1
    @test problem.metadata["layout_graph"] == "glued_trees_augmented"
    @test Int(problem.metadata["layout_temporary_triangulation_edge_count"]) > 0
    @test get(m.sampler_metadata, "sampler", nothing) == "tempered_mcmc"

    pos = compute_sfdp_layout(problem.num_vertices, problem.edges; seed=4, iterations=30, scale=1.0)
    render_pos = pos[Int.(problem.render_vertex_indices) .+ 1, :]

    @test size(pos) == (problem.num_vertices, 3)
    @test size(render_pos) == (num_vertices(m), 3)
    @test all(isfinite, render_pos)
end

@testset "Half-plane meandric pipeline smoke test" begin
    cfg = Dict{String,Any}(
        "model" => Dict{String,Any}(
            "type" => "half_plane_meandric",
            "pairs" => 6,
            "boundary_half_length" => 3,
            "seed" => 2,
        ),
        "layout" => Dict{String,Any}(
            "dimension" => 2,
            "engine" => "tutte",
            "solver" => "direct",
        ),
        "output" => Dict{String,Any}("show" => false),
    )

    timings = run_pipeline(cfg)
    @test timings isa TimingRecorder
end

@testset "Uniform meandric pipeline projects auxiliary layout back to strands" begin
    mktempdir() do out_dir
        cfg = Dict{String,Any}(
            "model" => Dict{String,Any}(
                "type" => "uniform_meandric",
                "pairs" => 6,
                "seed" => 3,
            ),
            "layout" => Dict{String,Any}(
                "dimension" => 3,
                "engine" => "sfdp",
                "iterations" => 40,
            ),
            "output" => Dict{String,Any}(
                "show" => false,
                "export_web" => out_dir,
            ),
        )

        timings = run_pipeline(cfg)
        meta = JSON3.read(read(joinpath(out_dir, "web_meta.json"), String))

        @test timings isa TimingRecorder
        @test String(meta["model"]) == "uniform_meandric"
        @test Int(meta["num_vertices"]) == 12
        @test Int(meta["layout_vertex_count"]) > Int(meta["render_vertex_count"])
        @test String(meta["render_vertex_projection"]) == "subset"
        @test String(meta["files"]["faces"]) == "data_faces.bin"
        @test any(String(group["name"]) == "loop 1" for group in meta["edge_groups"])
        @test count(group -> String(group["name"]) == "loop fade", meta["edge_groups"]) <= 1
    end
end

@testset "Uniform meander pipeline smoke test" begin
    mktempdir() do out_dir
        cfg = Dict{String,Any}(
            "model" => Dict{String,Any}(
                "type" => "uniform_meander",
                "pairs" => 6,
                "seed" => 9,
                "sweeps_per_temperature" => 3,
                "search_sweeps" => 12,
                "mixing_sweeps" => 3,
                "restarts" => 12,
            ),
            "layout" => Dict{String,Any}(
                "dimension" => 3,
                "engine" => "sfdp",
                "iterations" => 30,
            ),
            "output" => Dict{String,Any}(
                "show" => false,
                "export_web" => out_dir,
            ),
        )

        timings = run_pipeline(cfg)
        meta = JSON3.read(read(joinpath(out_dir, "web_meta.json"), String))

        @test timings isa TimingRecorder
        @test String(meta["model"]) == "uniform_meander"
        @test Int(meta["num_loops"]) == 1
        @test String(meta["files"]["faces"]) == "data_faces.bin"
    end
end

@testset "Closed meandric 3D web export matches half-plane palette" begin
    mktempdir() do out_dir
        m = build_uniform_meandric_system(; order=20, seed=17)
        problem = prepare_layout_problem(m; dimension=3)
        pos = compute_sfdp_layout(problem.num_vertices, problem.edges; seed=17, iterations=40, scale=1.0)
        render_pos = pos[Int.(problem.render_vertex_indices) .+ 1, :]

        export_web_binaries(
            problem.render_map_data,
            render_pos .* 10.0,
            out_dir;
            edge_groups=problem.edge_groups,
            faces=problem.faces,
            triangles=problem.surface_triangles,
            metadata=problem.metadata,
        )

        meta = JSON3.read(read(joinpath(out_dir, "web_meta.json"), String))
        groups = Dict(String(group["name"]) => group for group in meta["edge_groups"])

        @test String(groups["upper"]["color"]) == "#1f7a3a"
        @test String(groups["lower"]["color"]) == "#1f7a3a"
        @test String(groups["loop 1"]["color"]) == "#CC0000"
        @test String(groups["loop 2"]["color"]) == "#FFA500"
        @test String(groups["loop 3"]["color"]) == "#B3B300"
        @test String(groups["loop 4"]["color"]) == "#008000"
        @test String(groups["loop 5"]["color"]) == "#000080"
        if haskey(groups, "loop fade")
            @test String(groups["loop fade"]["kind"]) == "colored_index_pairs"
        end
    end
end

@testset "Meander 3D web export uses the same loop palette" begin
    mktempdir() do out_dir
        m = build_uniform_meander(
            ;
            order=20,
            seed=17,
            sweeps_per_temperature=4,
            search_sweeps=16,
            mixing_sweeps=4,
            restarts=12,
        )
        problem = prepare_layout_problem(m; dimension=3)
        pos = compute_sfdp_layout(problem.num_vertices, problem.edges; seed=17, iterations=40, scale=1.0)
        render_pos = pos[Int.(problem.render_vertex_indices) .+ 1, :]

        export_web_binaries(
            problem.render_map_data,
            render_pos .* 10.0,
            out_dir;
            edge_groups=problem.edge_groups,
            faces=problem.faces,
            triangles=problem.surface_triangles,
            metadata=problem.metadata,
        )

        meta = JSON3.read(read(joinpath(out_dir, "web_meta.json"), String))
        groups = Dict(String(group["name"]) => group for group in meta["edge_groups"])

        @test String(groups["upper"]["color"]) == "#1f7a3a"
        @test String(groups["lower"]["color"]) == "#1f7a3a"
        @test String(groups["loop 1"]["color"]) == "#CC0000"
        @test !haskey(groups, "loop fade")
    end
end

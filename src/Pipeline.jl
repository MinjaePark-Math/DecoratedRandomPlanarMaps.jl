# YAML-configurable pipeline helpers.

function _normalize_config_value(x)
    if x isa AbstractDict
        return Dict{String,Any}(string(k) => _normalize_config_value(v) for (k, v) in pairs(x))
    elseif x isa AbstractVector
        return Any[_normalize_config_value(v) for v in x]
    else
        return x
    end
end

function _as_string_dict(x)
    if x isa AbstractDict
        return Dict{String,Any}(string(k) => _normalize_config_value(v) for (k, v) in pairs(x))
    elseif x === nothing
        return Dict{String,Any}()
    else
        throw(ArgumentError("expected a mapping, got $(typeof(x))"))
    end
end

function _optional_string(x)
    if x === nothing
        return nothing
    end
    s = strip(string(x))
    return isempty(s) ? nothing : s
end

_optional_int(x) = x === nothing ? nothing : Int(x)
_optional_float(x) = x === nothing ? nothing : float(x)

function _optional_symbol(x)
    x === nothing && return nothing
    raw = lowercase(strip(string(x)))
    raw = replace(raw, '-' => '_')
    return Symbol(raw)
end

function _as_bool(x)
    if x isa Bool
        return x
    elseif x isa Integer
        return x != 0
    else
        s = lowercase(strip(string(x)))
        if s in ("1", "true", "yes", "y", "on")
            return true
        elseif s in ("0", "false", "no", "n", "off")
            return false
        end
    end
    throw(ArgumentError("could not interpret $(repr(x)) as a boolean"))
end

function _canonical_model_type(model_type)
    raw = lowercase(strip(string(model_type)))
    raw = replace(raw, '-' => '_')
    if raw in ("uniform", "quadrangulation", "uniform_quadrangulation")
        return "uniform"
    elseif raw in ("schnyder", "schnyder_wood", "schnyder_map")
        return "schnyder"
    elseif raw in ("fk", "hc", "fk_map", "hc_map", "fk_decorated_map", "hamburger_cheeseburger")
        return "fk"
    elseif raw in ("spanning_tree", "spanningtree", "tree", "tree_decorated_map")
        return "spanning_tree"
    else
        valid = join(["uniform", "schnyder", "fk", "spanning_tree"], ", ")
        throw(ArgumentError("unknown model type $(repr(model_type)); expected one of: $valid"))
    end
end

function load_config(path::AbstractString)
    data = YAML.load_file(string(path))
    data === nothing && return Dict{String,Any}()
    data isa AbstractDict || throw(ArgumentError("top-level config must be a mapping"))
    return _as_string_dict(data)
end

function build_map_from_config(model_cfg)
    cfg = _as_string_dict(model_cfg)
    raw_model_type = get(cfg, "type", nothing)
    raw_model_type === nothing && throw(ArgumentError("model.type is required"))
    model_type = _canonical_model_type(raw_model_type)
    seed = Int(get(cfg, "seed", 1))

    faces = if haskey(cfg, "faces")
        Int(cfg["faces"])
    else
        nothing
    end

    pairs = if haskey(cfg, "pairs")
        Int(cfg["pairs"])
    elseif haskey(cfg, "order")
        Int(cfg["order"])
    else
        faces
    end

    if model_type == "uniform"
        faces === nothing && throw(ArgumentError("model.faces is required"))
        return generate_uniform_map(; faces=faces, seed=seed)
    elseif model_type == "schnyder"
        faces === nothing && throw(ArgumentError("model.faces is required"))
        return generate_schnyder_map(; size_faces=faces, seed=seed)
    elseif model_type == "spanning_tree"
        faces === nothing && throw(ArgumentError("model.faces is required"))
        return build_spanning_tree_map(; faces=faces, seed=seed)
    else
        p = get(cfg, "p", nothing)
        q = get(cfg, "q", nothing)
        sampling_method = _optional_symbol(get(cfg, "sampling_method", nothing))
        pool_size = _optional_int(get(cfg, "pool_size", nothing))
        mh_steps = _optional_int(get(cfg, "mh_steps", nothing))
        exact_pilot_samples = Int(get(cfg, "exact_pilot_samples", 64))
        exact_min_acceptance = float(get(cfg, "exact_min_acceptance", 1.0e-4))
        max_exact_tries = Int(get(cfg, "max_exact_tries", 1_000_000))
        return build_fk_map(
            ;
            faces=faces,
            p=p,
            q=q,
            seed=seed,
            sampling_method=sampling_method === nothing ? :auto : sampling_method,
            pool_size=pool_size,
            mh_steps=mh_steps,
            exact_pilot_samples=exact_pilot_samples,
            exact_min_acceptance=exact_min_acceptance,
            max_exact_tries=max_exact_tries,
        )
    end
end

function _merged_metadata(parts...)
    out = Dict{String,Any}()
    for part in parts
        part === nothing && continue
        for (k, v) in pairs(part)
            out[string(k)] = v
        end
    end
    return out
end



function _model_metadata_from_config(model_type::AbstractString, model_cfg)
    cfg = _as_string_dict(model_cfg)
    md = Dict{String,Any}(
        "model" => string(model_type),
        "faces" => Int(get(cfg, "faces", 0)),
        "seed" => Int(get(cfg, "seed", 1)),
    )
    if model_type == "fk" || model_type == "spanning_tree"
        if haskey(cfg, "q") && get(cfg, "q", nothing) !== nothing
            md["q"] = float(cfg["q"])
        elseif haskey(cfg, "p") && get(cfg, "p", nothing) !== nothing
            md["p"] = float(cfg["p"])
            md["q"] = fk_q_from_p(float(cfg["p"]))
        elseif model_type == "spanning_tree"
            md["p"] = 0.0
            md["q"] = 0.0
        end
    end
    return md
end

function _default_preview_path(export_web, filename::AbstractString)
    export_web !== nothing && return joinpath(export_web, filename)
    return joinpath(pwd(), filename)
end

function _load_optional_package_in_main(pkg::AbstractString)
    sym = Symbol(pkg)
    isdefined(Main, sym) && return nothing
    try
        Core.eval(Main, :(using $sym))
        return nothing
    catch err
        return err
    end
end

function _ensure_makie_extension_loaded()
    ext = Base.get_extension(@__MODULE__, :DecoratedRandomPlanarMapsMakieExt)
    ext !== nothing && return ext

    missing = String[]
    errors = String[]
    for pkg in ("GLMakie", "GeometryBasics")
        err = _load_optional_package_in_main(pkg)
        if err !== nothing
            push!(missing, pkg)
            push!(errors, "$pkg => $(sprint(showerror, err))")
        end
    end

    ext = Base.get_extension(@__MODULE__, :DecoratedRandomPlanarMapsMakieExt)
    if ext === nothing && isdefined(Base, :retry_load_extensions)
        Base.retry_load_extensions()
        ext = Base.get_extension(@__MODULE__, :DecoratedRandomPlanarMapsMakieExt)
    end
    if ext === nothing
        if !isempty(missing)
            install_cmd = "using Pkg; Pkg.add([\"GLMakie\", \"GeometryBasics\"])"
            joined_errors = join(errors, " | ")
            detail = isempty(joined_errors) ? "" : " Original error: $joined_errors"
            throw(ArgumentError("`output.show: true` requires the optional Makie backend. Install `GLMakie` and `GeometryBasics` in the active environment, then try again. Example: `$install_cmd`.$detail"))
        end
        throw(ArgumentError("GLMakie and GeometryBasics were loaded, but the Makie extension did not activate. Please ensure both packages are installed in the active environment and restart the Julia session if needed."))
    end
    return ext
end

function run_pipeline(config_path::AbstractString)
    return run_pipeline(load_config(config_path))
end

function run_pipeline(cfg)
    cfg_dict = _as_string_dict(cfg)
    timings = TimingRecorder()

    model_cfg = _as_string_dict(get(cfg_dict, "model", Dict{String,Any}()))
    layout_cfg = _as_string_dict(get(cfg_dict, "layout", Dict{String,Any}()))
    output_cfg = _as_string_dict(get(cfg_dict, "output", Dict{String,Any}()))

    model_type = _canonical_model_type(get(model_cfg, "type", "uniform"))
    model_metadata = _model_metadata_from_config(model_type, model_cfg)

    println("Generating $(model_type) map...")
    map_data = track!(timings, "generate_map") do
        build_map_from_config(model_cfg)
    end
    println("   vertices=$(num_vertices(map_data)) faces=$(num_faces(map_data))")

    dimension = Int(get(layout_cfg, "dimension", 3))
    layout_seed = Int(get(layout_cfg, "seed", get(model_cfg, "seed", 7)))
    layout_options = _as_string_dict(get(layout_cfg, "options", Dict{String,Any}()))

    println("Preparing $(dimension)D layout problem...")
    problem = track!(timings, "prepare_layout") do
        prepare_layout_problem(
            map_data;
            dimension=dimension,
            boundary_scale=get(layout_cfg, "boundary_scale", nothing),
            options=layout_options,
        )
    end

    merged_problem_metadata = _merged_metadata(model_metadata, problem.metadata)

    boundary_count = problem.boundary_vertices === nothing ? 0 : length(problem.boundary_vertices)
    extra = dimension == 2 ? " boundary_vertices=$(boundary_count)" : ""
    println("   layout_edges=$(size(problem.edges, 1))$(extra)")

    base_pos = nothing
    circle_radii = nothing
    layout_metadata = Dict{String,Any}()
    engine_title = dimension == 2 ? "Tutte" : "SFDP"

    if dimension == 3
        engine = lowercase(strip(string(get(layout_cfg, "engine", "sfdp"))))
        engine == "sfdp" || throw(ArgumentError("unsupported 3D layout engine $(repr(engine)); only 'sfdp' is currently implemented"))
        layout_normalize_scale = float(get(layout_cfg, "normalize_scale", 1.0))
        sfdp_K = _optional_float(get(layout_cfg, "K", get(layout_options, "K", nothing)))
        sfdp_repulsiveforce = _optional_float(get(layout_cfg, "repulsiveforce", get(layout_options, "repulsiveforce", nothing)))
        sfdp_iterations = _optional_int(get(layout_cfg, "iterations", get(layout_options, "iterations", nothing)))
        sfdp_overlap = _optional_string(get(layout_cfg, "overlap", get(layout_options, "overlap", nothing)))

        println("Computing SFDP layout...")
        base_pos = track!(timings, "compute_layout") do
            compute_sfdp_layout(
                problem.num_vertices,
                problem.edges;
                scale=layout_normalize_scale,
                seed=layout_seed,
                K=sfdp_K,
                repulsiveforce=sfdp_repulsiveforce,
                iterations=sfdp_iterations,
                overlap=sfdp_overlap,
                options=layout_options,
            )
        end
        layout_metadata = Dict("engine" => "sfdp", "dimension" => 3)
        sfdp_K !== nothing && (layout_metadata["K"] = sfdp_K)
        sfdp_repulsiveforce !== nothing && (layout_metadata["repulsiveforce"] = sfdp_repulsiveforce)
        sfdp_iterations !== nothing && (layout_metadata["iterations"] = sfdp_iterations)
        sfdp_overlap !== nothing && (layout_metadata["overlap"] = sfdp_overlap)
    else
        engine = replace(lowercase(strip(string(get(layout_cfg, "engine", "tutte")))), '-' => '_')

        if engine == "tutte"
            engine_title = "Tutte"
            println("Computing Tutte embedding...")
            base_pos, tutte_meta = track!(timings, "compute_layout") do
                compute_tutte_layout(
                    problem.num_vertices,
                    problem.edges,
                    problem.boundary_vertices;
                    boundary_positions=problem.boundary_positions,
                    source_vertex=get(layout_cfg, "source_vertex", nothing),
                    seed=layout_seed,
                    radius=get(layout_cfg, "radius", nothing),
                    solver=string(get(layout_cfg, "solver", "auto")),
                    tol=float(get(layout_cfg, "tol", 1.0e-9)),
                    maxiter=_optional_int(get(layout_cfg, "maxiter", nothing)),
                    solver_options=layout_options,
                    return_metadata=true,
                )
            end

            layout_metadata = _merged_metadata(Dict("engine" => "tutte", "dimension" => 2), tutte_meta)
            if problem.boundary_vertices !== nothing
                embedding_info = get(layout_metadata, "embedding_info", Dict{String,Any}())
                solver_name = get(embedding_info, "solver", "unknown")
                boundary_mode = get(layout_metadata, "boundary_positions_mode", "unknown")
                println("   solver=$(solver_name) boundary_mode=$(boundary_mode)")
            end
        elseif engine == "circle_packing"
            engine_title = "Circle Packing"
            cp_maxiter = Int(get(layout_cfg, "maxiter", get(layout_cfg, "iterations", get(layout_options, "maxiter", get(layout_options, "iterations", 200)))))
            cp_tol = float(get(layout_cfg, "tol", get(layout_options, "tol", 1.0e-8)))
            cp_relaxation = float(get(layout_cfg, "relaxation", get(layout_options, "relaxation", 1.0)))
            cp_initial_radius = float(get(layout_cfg, "initial_radius", get(layout_options, "initial_radius", 0.5)))
            cp_min_radius = float(get(layout_cfg, "min_radius", get(layout_options, "min_radius", 1.0e-8)))

            println("Computing circle packing layout...")
            base_pos, circle_radii, packing_meta = track!(timings, "compute_layout") do
                compute_circle_packing_layout(
                    problem.num_vertices,
                    problem.edges,
                    problem.boundary_vertices;
                    faces=problem.faces,
                    triangles=problem.surface_triangles,
                    triangle_edge_ids=problem.surface_triangle_edge_ids,
                    maxiter=cp_maxiter,
                    tol=cp_tol,
                    relaxation=cp_relaxation,
                    initial_radius=cp_initial_radius,
                    min_radius=cp_min_radius,
                    return_metadata=true,
                )
            end

            layout_metadata = _merged_metadata(Dict("engine" => "circle_packing", "dimension" => 2), packing_meta)
            println("   converged=$(get(layout_metadata, "converged", false)) iterations=$(get(layout_metadata, "iterations", 0))")
        else
            throw(ArgumentError("unsupported 2D layout engine $(repr(engine)); expected 'tutte' or 'circle_packing'"))
        end
    end

    render_source = problem.render_map_data

    if problem.render_vertex_indices !== nothing
        render_idx = Int.(problem.render_vertex_indices) .+ 1
        minimum(render_idx) >= 1 || throw(ArgumentError("render_vertex_indices must be >= 0"))
        maximum(render_idx) <= size(base_pos, 1) || throw(ArgumentError("render_vertex_indices exceed the computed layout size"))
        base_pos = base_pos[render_idx, :]
        if circle_radii !== nothing
            circle_radii = circle_radii[render_idx]
        end
        layout_metadata["layout_vertex_count"] = get(layout_metadata, "layout_vertex_count", problem.num_vertices)
        layout_metadata["render_vertex_count"] = size(base_pos, 1)
        layout_metadata["render_vertex_projection"] = "subset"
    end

    export_web = _optional_string(get(output_cfg, "export_web", nothing))
    if export_web !== nothing
        web_scale = float(get(output_cfg, "web_scale", dimension == 3 ? 10.0 : 1.0))
        println("Exporting web binaries (scale=$(web_scale))...")
        track!(timings, "export_web") do
            export_web_binaries(
                render_source,
                base_pos .* web_scale,
                export_web;
                edge_groups=problem.edge_groups,
                faces=problem.faces,
                triangles=problem.surface_triangles,
                circle_radii=circle_radii === nothing ? nothing : circle_radii .* web_scale,
                metadata=_merged_metadata(merged_problem_metadata, layout_metadata),
            )
        end
    end

    export_stl = _optional_string(get(output_cfg, "export_stl", nothing))
    if export_stl !== nothing
        dimension == 3 || throw(ArgumentError("STL export is only supported for 3D layouts"))
        stl_scale = float(get(output_cfg, "stl_scale", 100.0))
        println("Exporting STL (scale=$(stl_scale))...")
        track!(timings, "export_stl") do
            export_stl_binary(map_data, base_pos .* stl_scale, export_stl)
        end
    end

    preview_svg = _optional_string(get(output_cfg, "preview_svg", nothing))
    if dimension == 2 && preview_svg !== nothing
        preview_scale = float(get(output_cfg, "show_scale", 1.0))
        println("Writing SVG preview (scale=$(preview_scale))...")
        track!(timings, "preview_svg") do
            export_svg_preview(
                render_source,
                base_pos .* preview_scale,
                preview_svg;
                edge_groups=problem.edge_groups,
                faces=problem.faces,
                triangles=problem.surface_triangles,
                circle_radii=circle_radii === nothing ? nothing : circle_radii .* preview_scale,
                title=default_render_title(render_source; metadata=_merged_metadata(merged_problem_metadata, layout_metadata), dimension=2, engine=engine_title),
                metadata=_merged_metadata(merged_problem_metadata, layout_metadata),
            )
        end
    end

    show_output = _as_bool(get(output_cfg, "show", true))
    if show_output
        if dimension == 2
            show_scale = float(get(output_cfg, "show_scale", 1.0))
            _ensure_makie_extension_loaded()
            println("Opening Makie viewer (scale=$(show_scale))...")
            track!(timings, "show_viewer") do
                Base.invokelatest(
                    render_makie_2d,
                    render_source,
                    base_pos .* show_scale,
                    ;
                    edge_groups=problem.edge_groups,
                    faces=problem.faces,
                    triangles=problem.surface_triangles,
                    circle_radii=circle_radii === nothing ? nothing : circle_radii .* show_scale,
                    title=default_render_title(render_source; metadata=_merged_metadata(merged_problem_metadata, layout_metadata), dimension=2, engine=engine_title),
                    metadata=_merged_metadata(merged_problem_metadata, layout_metadata),
                )
            end
        else
            show_scale = float(get(output_cfg, "show_scale", 10.0))
            _ensure_makie_extension_loaded()
            println("Opening Makie viewer (scale=$(show_scale))...")
            track!(timings, "show_viewer") do
                Base.invokelatest(
                    render_makie_3d,
                    render_source,
                    base_pos .* show_scale,
                    ;
                    edge_groups=problem.edge_groups,
                    faces=problem.faces,
                    triangles=problem.surface_triangles,
                    title=default_render_title(render_source; metadata=_merged_metadata(merged_problem_metadata, layout_metadata), dimension=3, engine=engine_title),
                    metadata=_merged_metadata(merged_problem_metadata, layout_metadata),
                )
            end
        end
    end

    println("Timing summary:")
    for line in pretty_lines(timings)
        println("   $(line)")
    end

    timings_json = _optional_string(get(output_cfg, "timings_json", nothing))
    if timings_json !== nothing
        path = write_json(timings, timings_json)
        println("Saved timings to $(path)")
    end

    return timings
end

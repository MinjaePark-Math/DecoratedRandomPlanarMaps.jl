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
    elseif raw in ("half_plane_meandric", "half_plane_meandric_system", "halfplane_meandric", "halfplane_meandric_system")
        return "half_plane_meandric"
    elseif raw in ("uniform_meandric", "uniform_meandric_system", "meandric", "meandric_system")
        return "uniform_meandric"
    elseif raw in ("uniform_meander", "meander")
        return "uniform_meander"
    else
        valid = join(["uniform", "schnyder", "fk", "spanning_tree", "half_plane_meandric", "uniform_meandric", "uniform_meander"], ", ")
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
    elseif model_type == "half_plane_meandric"
        pairs === nothing && throw(ArgumentError("model.pairs or model.order is required"))
        return build_half_plane_meandric_system(
            ;
            order=pairs,
            boundary_half_length=get(cfg, "boundary_half_length", nothing),
            seed=seed,
        )
    elseif model_type == "uniform_meandric"
        pairs === nothing && throw(ArgumentError("model.pairs or model.order is required"))
        return build_uniform_meandric_system(; order=pairs, seed=seed)
    elseif model_type == "uniform_meander"
        pairs === nothing && throw(ArgumentError("model.pairs or model.order is required"))
        temper_schedule = if haskey(cfg, "temper_schedule") && get(cfg, "temper_schedule", nothing) !== nothing
            Float64[float(v) for v in cfg["temper_schedule"]]
        else
            nothing
        end
        return build_uniform_meander(
            ;
            order=pairs,
            seed=seed,
            temper_schedule=temper_schedule,
            sweeps_per_temperature=_optional_int(get(cfg, "sweeps_per_temperature", nothing)),
            search_sweeps=_optional_int(get(cfg, "search_sweeps", nothing)),
            mixing_sweeps=_optional_int(get(cfg, "mixing_sweeps", nothing)),
            restarts=Int(get(cfg, "restarts", 8)),
        )
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
    size_param = if haskey(cfg, "faces")
        Int(cfg["faces"])
    elseif haskey(cfg, "pairs")
        Int(cfg["pairs"])
    elseif haskey(cfg, "order")
        Int(cfg["order"])
    else
        0
    end
    md = Dict{String,Any}(
        "model" => string(model_type),
        "faces" => size_param,
        "seed" => Int(get(cfg, "seed", 1)),
    )
    size_param > 0 && (md["pairs"] = size_param)
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
    elseif model_type == "half_plane_meandric"
        if haskey(cfg, "boundary_half_length") && get(cfg, "boundary_half_length", nothing) !== nothing
            md["boundary_half_length"] = Int(cfg["boundary_half_length"])
        end
    elseif model_type == "uniform_meander"
        md["target_loops"] = 1
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

function _circle_packing_problem_inputs(problem)
    boundary_vertices = problem.packing_boundary_vertices === nothing ? problem.boundary_vertices : problem.packing_boundary_vertices
    faces = problem.packing_faces === nothing ? problem.faces : problem.packing_faces
    triangles = problem.packing_triangles === nothing ? problem.surface_triangles : problem.packing_triangles
    triangle_edge_ids = problem.packing_triangle_edge_ids === nothing ? problem.surface_triangle_edge_ids : problem.packing_triangle_edge_ids
    return boundary_vertices, faces, triangles, triangle_edge_ids
end

function _append_edge_group_edges(edge_groups::Dict{String,Matrix{Int32}}, name::AbstractString, extra_edges)
    arr = sanitize_edge_array(extra_edges)
    size(arr, 1) == 0 && return edge_groups
    out = Dict{String,Matrix{Int32}}(k => copy(v) for (k, v) in edge_groups)
    key = string(name)
    if haskey(out, key)
        out[key] = vcat(out[key], arr)
    else
        out[key] = arr
    end
    return out
end

function _recover_outer_vertex_edge_groups(map_data::FKMap, outer_vertex::Int32)
    recovered = Dict{String,Matrix{Int32}}()
    raw = raw_active_edges_by_color(map_data; drop_loops=true)
    for color in (GREEN, RED, BLUE, PURPLE, ORANGE)
        arr = get(raw, Int(color), Matrix{Int32}(undef, 0, 2))
        size(arr, 1) == 0 && continue
        keep = BitVector(undef, size(arr, 1))
        for i in 1:size(arr, 1)
            keep[i] = arr[i, 1] == outer_vertex || arr[i, 2] == outer_vertex
        end
        selected = arr[keep, :]
        size(selected, 1) == 0 && continue
        name = FK_COLOR_NAME[Int(color)]
        name == "green" && (name = "generic")
        recovered = _append_edge_group_edges(recovered, name, selected)
    end
    return recovered
end

function _replace_outer_vertex_edges_with_generic(edge_groups::Dict{String,Matrix{Int32}}, outer_vertex::Int32, boundary_vertices)
    boundary = Int32.(collect(boundary_vertices))
    isempty(boundary) && return edge_groups
    star_edges = sanitize_edge_array(hcat(fill(outer_vertex, length(boundary)), boundary))

    star_keys = Set{Tuple{Int32,Int32}}()
    for i in 1:size(star_edges, 1)
        u = star_edges[i, 1]
        v = star_edges[i, 2]
        push!(star_keys, (min(u, v), max(u, v)))
    end

    out = Dict{String,Matrix{Int32}}()
    for (name, edge_array) in edge_groups
        arr = sanitize_edge_array(edge_array)
        size(arr, 1) == 0 && continue
        keep = BitVector(undef, size(arr, 1))
        for i in 1:size(arr, 1)
            key = (min(arr[i, 1], arr[i, 2]), max(arr[i, 1], arr[i, 2]))
            keep[i] = !(key in star_keys)
        end
        filtered = arr[keep, :]
        size(filtered, 1) == 0 && continue
        out[string(name)] = filtered
    end

    generic_existing = get(out, "generic", Matrix{Int32}(undef, 0, 2))
    out["generic"] = first(collapse_undirected_edges(vcat(generic_existing, star_edges); drop_loops=true))
    return out
end

function _remap_edge_groups_for_render(edge_groups::Dict{String,Matrix{Int32}}, mapping)
    packed_to_render = Int32.(collect(mapping))
    isempty(packed_to_render) && return Dict{String,Matrix{Int32}}(k => copy(v) for (k, v) in edge_groups)

    out = Dict{String,Matrix{Int32}}()
    for (name, edge_array) in edge_groups
        arr = sanitize_edge_array(edge_array)
        if size(arr, 1) == 0
            out[string(name)] = arr
            continue
        end
        remapped = Matrix{Int32}(undef, size(arr, 1), 2)
        for i in 1:size(arr, 1)
            remapped[i, 1] = packed_to_render[Int(arr[i, 1]) + 1]
            remapped[i, 2] = packed_to_render[Int(arr[i, 2]) + 1]
        end
        out[string(name)] = remapped
    end
    return out
end

function _fan_triangles_for_boundary(boundary_vertices, outer_vertex::Integer)
    boundary = Int32.(collect(boundary_vertices))
    length(boundary) >= 3 || return Matrix{Int32}(undef, 0, 3)
    fan = Matrix{Int32}(undef, length(boundary), 3)
    for i in 1:length(boundary)
        fan[i, 1] = Int32(outer_vertex)
        fan[i, 2] = boundary[i]
        fan[i, 3] = boundary[mod1(i + 1, length(boundary))]
    end
    return sanitize_triangles(fan; drop_degenerate=true, deduplicate=true)
end

function _augment_sphere_render_triangles(triangles, render_count::Integer, layout_metadata)
    triangles === nothing && return Matrix{Int32}(undef, 0, 3)
    tri = sanitize_triangles(triangles; drop_degenerate=true, deduplicate=false)
    topo = get(layout_metadata, "packing_topology", nothing)
    topo == "sphere" || return tri

    if haskey(layout_metadata, "circle_packing_removed_outer_vertex") && haskey(layout_metadata, "circle_packing_packing_to_render_vertices")
        outer_vertex = Int(get(layout_metadata, "circle_packing_removed_outer_vertex", 0))
        mapping = Int.(collect(layout_metadata["circle_packing_packing_to_render_vertices"]))
        remapped = _reindex_triangles(tri, Dict(i - 1 => mapping[i] for i in eachindex(mapping)); deduplicate=false)
        if haskey(layout_metadata, "circle_packing_outer_boundary_vertices")
            fan = _fan_triangles_for_boundary(layout_metadata["circle_packing_outer_boundary_vertices"], outer_vertex)
            remapped = vcat(remapped, fan)
        end
        return sanitize_triangles(remapped; drop_degenerate=true, deduplicate=true)
    end

    if haskey(layout_metadata, "circle_packing_outer_boundary_vertices")
        outer_vertex = Int(render_count) - 1
        fan = _fan_triangles_for_boundary(layout_metadata["circle_packing_outer_boundary_vertices"], outer_vertex)
        return sanitize_triangles(vcat(tri, fan); drop_degenerate=true, deduplicate=true)
    end

    return tri
end

function _expand_sphere_circle_geometry_for_render(sphere_circle_geometry, render_count::Integer, mapping, outer_vertex::Integer, layout_metadata)
    sphere_circle_geometry === nothing && return nothing
    render_n = Int(render_count)
    render_n >= 0 || throw(ArgumentError("render_count must be nonnegative"))
    packed_to_render = Int.(collect(mapping))
    length(packed_to_render) == size(sphere_circle_geometry["normals"], 1) ||
        throw(ArgumentError("sphere circle geometry must align with the packed render mapping"))

    centers = fill(NaN, render_n, 3)
    normals = fill(NaN, render_n, 3)
    radii = fill(NaN, render_n)
    offsets = fill(NaN, render_n)

    for (i, render_vertex) in enumerate(packed_to_render)
        centers[render_vertex + 1, :] .= sphere_circle_geometry["centers"][i, :]
        normals[render_vertex + 1, :] .= sphere_circle_geometry["normals"][i, :]
        radii[render_vertex + 1] = sphere_circle_geometry["radii"][i]
        offsets[render_vertex + 1] = sphere_circle_geometry["offsets"][i]
    end

    outer_geometry = get(layout_metadata, "sphere_outer_circle_geometry", nothing)
    outer_geometry === nothing && haskey(sphere_circle_geometry, "outer_geometry") && (outer_geometry = sphere_circle_geometry["outer_geometry"])
    outer_geometry === nothing && throw(ArgumentError("sphere circle render requires explicit outer-circle geometry"))
    centers[Int(outer_vertex) + 1, :] .= Float64.(outer_geometry["centers"][1, :])
    normals[Int(outer_vertex) + 1, :] .= Float64.(outer_geometry["normals"][1, :])
    radii[Int(outer_vertex) + 1] = float(outer_geometry["radii"][1])
    offsets[Int(outer_vertex) + 1] = float(outer_geometry["offsets"][1])

    out = Dict{String,Any}(
        "centers" => centers,
        "normals" => normals,
        "radii" => radii,
        "offsets" => offsets,
        "sphere_radius" => float(get(layout_metadata, "sphere_radius", sphere_circle_geometry["sphere_radius"])),
        "projection_scale" => get(layout_metadata, "sphere_effective_projection_scale", get(layout_metadata, "sphere_projection_scale", get(sphere_circle_geometry, "projection_scale", 1.0))),
        "projection" => get(sphere_circle_geometry, "projection", "inverse_stereographic"),
    )
    outer_geometry !== nothing && (out["outer_geometry"] = outer_geometry)
    return out
end

function _augment_sphere_circle_render(base_pos, edge_groups, sphere_circle_geometry, layout_metadata, map_data=nothing)
    topo = get(layout_metadata, "packing_topology", nothing)
    topo == "sphere" || return base_pos, edge_groups, sphere_circle_geometry

    sphere_radius = float(get(layout_metadata, "sphere_radius", 1.0))

    function outer_vertex_position(render_geometry, outer_vertex)
        normals = Float64.(render_geometry["normals"])
        return -sphere_radius .* normals[Int(outer_vertex) + 1, :]
    end

    if haskey(layout_metadata, "circle_packing_removed_outer_vertex") && haskey(layout_metadata, "circle_packing_packing_to_render_vertices")
        outer_vertex = Int(get(layout_metadata, "circle_packing_removed_outer_vertex", 0))
        boundary = Int32.(collect(get(layout_metadata, "circle_packing_outer_boundary_vertices", Int[])))
        render_count = Int(get(layout_metadata, "circle_packing_render_vertex_count", size(base_pos, 1) + 1))
        mapping = Int.(collect(layout_metadata["circle_packing_packing_to_render_vertices"]))
        length(mapping) == size(base_pos, 1) || throw(ArgumentError("circle_packing_packing_to_render_vertices must align with the computed sphere layout"))
        remapped_edge_groups = _remap_edge_groups_for_render(edge_groups, mapping)
        render_pos = Matrix{Float64}(undef, render_count, size(base_pos, 2))
        fill!(render_pos, NaN)
        for (i, original_vertex) in enumerate(mapping)
            render_pos[original_vertex + 1, :] .= base_pos[i, :]
        end
        render_edge_groups = if map_data isa FKMap
            recovered = _recover_outer_vertex_edge_groups(map_data, Int32(outer_vertex))
            out = remapped_edge_groups
            for (name, arr) in recovered
                out = _append_edge_group_edges(out, name, arr)
            end
            out
        else
            _replace_outer_vertex_edges_with_generic(remapped_edge_groups, Int32(outer_vertex), boundary)
        end
        render_geometry = _expand_sphere_circle_geometry_for_render(
            sphere_circle_geometry,
            render_count,
            mapping,
            outer_vertex,
            layout_metadata,
        )
        render_pos[outer_vertex + 1, :] .= outer_vertex_position(render_geometry, outer_vertex)
        layout_metadata["circle_packing_render_outer_vertex"] = outer_vertex
        layout_metadata["sphere_circle_geometry_includes_outer"] = true
        return render_pos, render_edge_groups, render_geometry
    end

    if haskey(layout_metadata, "circle_packing_outer_boundary_vertices")
        boundary = Int32.(collect(layout_metadata["circle_packing_outer_boundary_vertices"]))
        outer_vertex = Int32(size(base_pos, 1))
        render_pos = vcat(base_pos, reshape(zeros(Float64, 3), 1, :))
        extra_edges = hcat(fill(outer_vertex, length(boundary)), boundary)
        render_geometry = _expand_sphere_circle_geometry_for_render(
            sphere_circle_geometry,
            size(render_pos, 1),
            collect(0:(size(base_pos, 1) - 1)),
            outer_vertex,
            layout_metadata,
        )
        render_pos[Int(outer_vertex) + 1, :] .= outer_vertex_position(render_geometry, outer_vertex)
        layout_metadata["circle_packing_render_outer_vertex"] = Int(outer_vertex)
        layout_metadata["circle_packing_render_vertex_count"] = size(render_pos, 1)
        layout_metadata["sphere_circle_geometry_includes_outer"] = true
        return render_pos, _append_edge_group_edges(edge_groups, "generic", extra_edges), render_geometry
    end

    return base_pos, edge_groups, sphere_circle_geometry
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
    layout_engine = replace(lowercase(strip(string(get(layout_cfg, "engine", dimension == 3 ? "sfdp" : "tutte")))), '-' => '_')
    layout_options = copy(layout_options)
    layout_options["engine"] = layout_engine

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
    sphere_circle_geometry = nothing
    layout_metadata = Dict{String,Any}()
    engine_title = dimension == 2 ? "Tutte" : "SFDP"

    if dimension == 3
        engine = layout_engine
        layout_normalize_scale = float(get(layout_cfg, "normalize_scale", 1.0))
        if engine == "sfdp"
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
        elseif engine == "circle_packing"
            engine_title = "Circle Packing"
            cp_boundary, cp_faces, cp_triangles, cp_triangle_edge_ids = _circle_packing_problem_inputs(problem)
            cp_boundary === nothing && throw(ArgumentError("3D circle packing requires a sphere-topology packing boundary, but this layout problem does not provide one"))
            cp_maxiter = Int(get(layout_cfg, "maxiter", get(layout_cfg, "iterations", get(layout_options, "maxiter", get(layout_options, "iterations", 200)))))
            cp_tol = float(get(layout_cfg, "tol", get(layout_options, "tol", 1.0e-8)))
            cp_relaxation = float(get(layout_cfg, "relaxation", get(layout_options, "relaxation", 1.0)))
            cp_initial_radius = float(get(layout_cfg, "initial_radius", get(layout_options, "initial_radius", 0.5)))
            cp_min_radius = float(get(layout_cfg, "min_radius", get(layout_options, "min_radius", 1.0e-8)))
            cp_projection_scale = float(get(layout_cfg, "sphere_projection_scale", get(layout_options, "sphere_projection_scale", 1.0)))

            println("Computing sphere circle packing layout...")
            base_pos, circle_radii, packing_meta = track!(timings, "compute_layout") do
                compute_circle_packing_layout(
                    problem.num_vertices,
                    problem.edges,
                    cp_boundary;
                    faces=cp_faces,
                    triangles=cp_triangles,
                    triangle_edge_ids=cp_triangle_edge_ids,
                    maxiter=cp_maxiter,
                    tol=cp_tol,
                    relaxation=cp_relaxation,
                    initial_radius=cp_initial_radius,
                    min_radius=cp_min_radius,
                    project_to_sphere=true,
                    sphere_radius=layout_normalize_scale,
                    sphere_projection_scale=cp_projection_scale,
                    return_metadata=true,
                )
            end
            sphere_circle_geometry = pop!(packing_meta, "sphere_circle_geometry", nothing)
            layout_metadata = _merged_metadata(Dict("engine" => "circle_packing", "dimension" => 3), packing_meta)
            println("   converged=$(get(layout_metadata, "converged", false)) iterations=$(get(layout_metadata, "iterations", 0))")
        else
            throw(ArgumentError("unsupported 3D layout engine $(repr(engine)); expected 'sfdp' or 'circle_packing'"))
        end
    else
        engine = layout_engine

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
            cp_boundary, cp_faces, cp_triangles, cp_triangle_edge_ids = _circle_packing_problem_inputs(problem)
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
                    cp_boundary;
                    faces=cp_faces,
                    triangles=cp_triangles,
                    triangle_edge_ids=cp_triangle_edge_ids,
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
    render_edge_groups = problem.edge_groups
    render_faces = problem.faces
    render_triangles = problem.surface_triangles
    web_render_faces = render_faces
    web_render_triangles = render_triangles

    if problem.render_vertex_indices !== nothing
        render_idx = Int.(problem.render_vertex_indices) .+ 1
        minimum(render_idx) >= 1 || throw(ArgumentError("render_vertex_indices must be >= 0"))
        maximum(render_idx) <= size(base_pos, 1) || throw(ArgumentError("render_vertex_indices exceed the computed layout size"))
        base_pos = base_pos[render_idx, :]
        if circle_radii !== nothing
            circle_radii = circle_radii[render_idx]
        end
        if sphere_circle_geometry !== nothing
            sphere_circle_geometry = subset_sphere_circle_geometry(sphere_circle_geometry, render_idx)
        end
        layout_metadata["layout_vertex_count"] = get(layout_metadata, "layout_vertex_count", problem.num_vertices)
        layout_metadata["render_vertex_count"] = size(base_pos, 1)
        layout_metadata["render_vertex_projection"] = "subset"
    end

    render_metadata = _merged_metadata(merged_problem_metadata, layout_metadata)
    render_pos = base_pos
    if sphere_circle_geometry !== nothing
        render_pos, render_edge_groups, sphere_circle_geometry = _augment_sphere_circle_render(base_pos, render_edge_groups, sphere_circle_geometry, render_metadata, render_source)
        web_render_triangles = _augment_sphere_render_triangles(web_render_triangles, size(render_pos, 1), render_metadata)
        circle_radii = sphere_circle_geometry === nothing ? circle_radii : Float64.(sphere_circle_geometry["radii"])
        render_faces = nothing
        render_triangles = Matrix{Int32}(undef, 0, 3)
    end

    export_web = _optional_string(get(output_cfg, "export_web", nothing))
    if export_web !== nothing
        web_scale = float(get(output_cfg, "web_scale", dimension == 3 ? 10.0 : 1.0))
        println("Exporting web binaries (scale=$(web_scale))...")
        track!(timings, "export_web") do
            export_web_binaries(
                render_source,
                render_pos .* web_scale,
                export_web;
                edge_groups=render_edge_groups,
                faces=web_render_faces,
                triangles=web_render_triangles,
                circle_radii=circle_radii === nothing ? nothing : circle_radii .* web_scale,
                sphere_circle_geometry=scale_sphere_circle_geometry(sphere_circle_geometry, web_scale),
                metadata=render_metadata,
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
                render_pos .* preview_scale,
                preview_svg;
                edge_groups=render_edge_groups,
                faces=render_faces,
                triangles=render_triangles,
                circle_radii=circle_radii === nothing ? nothing : circle_radii .* preview_scale,
                title=default_render_title(render_source; metadata=render_metadata, dimension=2, engine=engine_title),
                metadata=render_metadata,
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
                    render_pos .* show_scale,
                    ;
                    edge_groups=render_edge_groups,
                    faces=render_faces,
                    triangles=render_triangles,
                    circle_radii=circle_radii === nothing ? nothing : circle_radii .* show_scale,
                    title=default_render_title(render_source; metadata=render_metadata, dimension=2, engine=engine_title),
                    metadata=render_metadata,
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
                    render_pos .* show_scale,
                    ;
                    edge_groups=render_edge_groups,
                    faces=render_faces,
                    triangles=render_triangles,
                    sphere_circle_geometry=scale_sphere_circle_geometry(sphere_circle_geometry, show_scale),
                    title=default_render_title(render_source; metadata=render_metadata, dimension=3, engine=engine_title),
                    metadata=render_metadata,
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

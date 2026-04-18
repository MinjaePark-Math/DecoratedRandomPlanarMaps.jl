# Shared helpers for converting combinatorial maps into drawable geometry.

const EDGE_COLOR_HINTS = Dict(
    "generic" => "#b0b7bf",
    "green" => "#1f7a3a",
    "red" => "#c0392b",
    "blue" => "#2980b9",
    "upper" => "#1f7a3a",
    "lower" => "#1f7a3a",
    "tree_upper" => "#1f7a3a",
    "tree_lower" => "#d35400",
    "glue_upper" => "#7f8c8d",
    "glue_lower" => "#95a5a6",
    "purple" => "#8e44ad",
    "orange" => "#d35400",
    "circles" => "#f4c430",
    "outer" => "#2c3e50",
    "navy" => "#2c3e50",
    "exploration" => "#27ae60",
    # Matches resource/meandric-systems:
    # palette(cgrad(:darktest, rev=true), 7)[1:5]
    "zz_loop_1" => "#CC0000",
    "zz_loop_2" => "#FFA500",
    "zz_loop_3" => "#B3B300",
    "zz_loop_4" => "#008000",
    "zz_loop_5" => "#000080",
)

const MEANDRIC_HIGHLIGHT_GROUPS = (
    "zz_loop_1",
    "zz_loop_2",
    "zz_loop_3",
    "zz_loop_4",
    "zz_loop_5",
)

const MEANDRIC_LOOP_FADE_PREFIX = "zy_loop_fade_"
const MEANDRIC_LOOP_COLOR_CUTOFF = length(MEANDRIC_HIGHLIGHT_GROUPS)
const MEANDRIC_FADE_TARGET_RANK = 300

function _hex_channel_pair_to_int(s::AbstractString, start::Int)
    return parse(Int, s[start:(start + 1)]; base=16)
end

function _hex_to_rgb(hex::AbstractString)
    s = strip(String(hex))
    startswith(s, "#") && (s = s[2:end])
    length(s) == 6 || throw(ArgumentError("hex colors must have 6 digits"))
    return (
        _hex_channel_pair_to_int(s, 1),
        _hex_channel_pair_to_int(s, 3),
        _hex_channel_pair_to_int(s, 5),
    )
end

function _rgb_to_hex(rgb::NTuple{3,Int})
    return "#" * uppercase(string(rgb[1], base=16, pad=2) * string(rgb[2], base=16, pad=2) * string(rgb[3], base=16, pad=2))
end

function _hex_lerp(a::AbstractString, b::AbstractString, t::Real)
    ta = clamp(float(t), 0.0, 1.0)
    ar, ag, ab = _hex_to_rgb(a)
    br, bg, bb = _hex_to_rgb(b)
    rgb = (
        round(Int, (1 - ta) * ar + ta * br),
        round(Int, (1 - ta) * ag + ta * bg),
        round(Int, (1 - ta) * ab + ta * bb),
    )
    return _rgb_to_hex(rgb)
end

function _hex_blend_over(fg::AbstractString, bg::AbstractString, alpha::Real)
    a = clamp(float(alpha), 0.0, 1.0)
    fr, fg_g, fb = _hex_to_rgb(fg)
    br, bg_g, bb = _hex_to_rgb(bg)
    rgb = (
        round(Int, a * fr + (1 - a) * br),
        round(Int, a * fg_g + (1 - a) * bg_g),
        round(Int, a * fb + (1 - a) * bb),
    )
    return _rgb_to_hex(rgb)
end

function _meandric_fade_group_name(rank::Integer)
    return MEANDRIC_LOOP_FADE_PREFIX * lpad(string(Int(rank)), 4, '0')
end

function _meandric_fade_fraction(rank::Integer; cutoff::Integer=typemax(Int), start_rank::Integer=length(MEANDRIC_HIGHLIGHT_GROUPS) + 1)
    lo = max(Int(start_rank), 1)
    hi = max(lo, Int(cutoff))
    r = clamp(Int(rank), lo, hi)
    hi == lo && return 0.0
    return (r - lo) / (hi - lo)
end

function _meandric_fade_opacity(rank::Integer; cutoff::Integer=MEANDRIC_FADE_TARGET_RANK, start_rank::Integer=length(MEANDRIC_HIGHLIGHT_GROUPS) + 1)
    lo = max(Int(start_rank), 1)
    hi = max(lo, Int(cutoff))
    r = clamp(Int(rank), lo, hi)
    hi == lo && return 0.98
    t = (r - lo) / (hi - lo)
    return 0.98 - 0.82 * t
end

function _meandric_fade_width_scale(rank::Integer; cutoff::Integer=MEANDRIC_FADE_TARGET_RANK, start_rank::Integer=length(MEANDRIC_HIGHLIGHT_GROUPS) + 1)
    lo = max(Int(start_rank), 1)
    hi = max(lo, Int(cutoff))
    r = clamp(Int(rank), lo, hi)
    hi == lo && return 1.55
    t = (r - lo) / (hi - lo)
    return 1.55 - 0.45 * t
end

function _meandric_rank_opacity(rank::Integer; cutoff::Integer=length(MEANDRIC_HIGHLIGHT_GROUPS))
    cutoff_int = max(1, Int(cutoff))
    rank_int = clamp(Int(rank), 1, cutoff_int)
    cutoff_int == 1 && return 1.0
    return 1.0 - 0.08 * (rank_int - 1) / (cutoff_int - 1)
end

function _edge_group_style(name::AbstractString)
    normalized = _normalize_group_string(name)

    if haskey(EDGE_COLOR_HINTS, normalized)
        if startswith(normalized, "zz_loop_")
            loop_idx = parse(Int, split(normalized, "_")[end])
            return (
                color=EDGE_COLOR_HINTS[normalized],
                opacity=_meandric_rank_opacity(loop_idx),
                width_scale=3.9,
                force_wide=true,
            )
        elseif normalized in ("upper", "lower")
            return (
                color=EDGE_COLOR_HINTS[normalized],
                opacity=0.42,
                width_scale=1.4,
                force_wide=false,
            )
        elseif normalized == "generic"
            return (
                color=EDGE_COLOR_HINTS[normalized],
                opacity=0.8,
                width_scale=1.2,
                force_wide=false,
            )
        end
        return (
            color=EDGE_COLOR_HINTS[normalized],
            opacity=0.9,
            width_scale=(normalized in ("red", "blue", "purple", "orange", "green", "navy")) ? 1.8 : 1.0,
            force_wide=false,
        )
    end

    fade_match = match(r"^zy_loop_fade_(\d+)$", normalized)
    if fade_match !== nothing
        rank = parse(Int, fade_match.captures[1])
        t = _meandric_fade_fraction(rank; cutoff=MEANDRIC_FADE_TARGET_RANK)
        return (
            color=_hex_lerp("#87CEEB", "#800080", t),
            opacity=_meandric_fade_opacity(rank; cutoff=MEANDRIC_FADE_TARGET_RANK),
            width_scale=_meandric_fade_width_scale(rank; cutoff=MEANDRIC_FADE_TARGET_RANK),
            force_wide=true,
        )
    end

    return (
        color="#2c3e50",
        opacity=0.9,
        width_scale=1.0,
        force_wide=false,
    )
end

edge_group_color(name::AbstractString) = _edge_group_style(name).color
edge_group_opacity(name::AbstractString) = _edge_group_style(name).opacity
edge_group_width_scale(name::AbstractString) = _edge_group_style(name).width_scale
edge_group_force_wide(name::AbstractString) = _edge_group_style(name).force_wide

function sanitize_triangles(triangles; drop_degenerate::Bool=true, deduplicate::Bool=true)
    arr = Int32.(triangles)
    isempty(arr) && return _empty_triangles()
    size(arr, 2) == 3 || throw(ArgumentError("triangle array must have shape (T, 3)"))

    keep = trues(size(arr, 1))
    for i in axes(arr, 1)
        tri = arr[i, :]
        if any(tri .< 0)
            keep[i] = false
            continue
        end
        if drop_degenerate
            keep[i] = (tri[1] != tri[2]) && (tri[2] != tri[3]) && (tri[1] != tri[3])
        end
    end
    arr = arr[keep, :]
    size(arr, 1) == 0 && return _empty_triangles()

    if deduplicate
        seen = Set{NTuple{3,Int32}}()
        unique_rows = NTuple{3,Int32}[]
        sizehint!(unique_rows, size(arr, 1))
        for i in axes(arr, 1)
            tri = (arr[i, 1], arr[i, 2], arr[i, 3])
            key_sorted = sort(collect(tri))
            key = (key_sorted[1], key_sorted[2], key_sorted[3])
            key in seen && continue
            push!(seen, key)
            push!(unique_rows, tri)
        end
        arr = Matrix{Int32}(undef, length(unique_rows), 3)
        for (i, tri) in enumerate(unique_rows)
            arr[i, 1] = tri[1]
            arr[i, 2] = tri[2]
            arr[i, 3] = tri[3]
        end
    end

    return arr
end

function _iter_faces(faces)
    if faces isa AbstractMatrix
        return [Int32[v for v in row if v >= 0] for row in eachrow(Int32.(faces))]
    else
        return [Int32[Int(v) for v in face if Int(v) >= 0] for face in faces]
    end
end

function fan_triangulate_faces(faces; drop_degenerate::Bool=true, deduplicate::Bool=true)
    triangles = NTuple{3,Int32}[]
    for face in _iter_faces(faces)
        length(face) < 3 && continue
        anchor = face[1]
        for i in 2:(length(face)-1)
            push!(triangles, (anchor, face[i], face[i+1]))
        end
    end
    if isempty(triangles)
        return _empty_triangles()
    end
    arr = Matrix{Int32}(undef, length(triangles), 3)
    for (i, tri) in enumerate(triangles)
        arr[i, 1] = tri[1]
        arr[i, 2] = tri[2]
        arr[i, 3] = tri[3]
    end
    return sanitize_triangles(arr; drop_degenerate=drop_degenerate, deduplicate=deduplicate)
end

function surface_triangles(map_data=nothing; faces=nothing, triangles=nothing, drop_degenerate::Bool=true, deduplicate::Bool=true)
    if triangles !== nothing
        return sanitize_triangles(triangles; drop_degenerate=drop_degenerate, deduplicate=deduplicate)
    elseif faces !== nothing
        return fan_triangulate_faces(faces; drop_degenerate=drop_degenerate, deduplicate=deduplicate)
    elseif map_data === nothing
        return _empty_triangles()
    elseif map_data isa FKMap
        return sanitize_triangles(map_data.triangulation_faces; drop_degenerate=drop_degenerate, deduplicate=deduplicate)
    elseif map_data isa MatedCRTMap
        faces = map_data.topology == :disk && Int(map_data.outer_face_index) >= 0 ?
            [map_data.faces[i] for i in eachindex(map_data.faces) if i != Int(map_data.outer_face_index) + 1] :
            map_data.faces
        return fan_triangulate_faces(faces; drop_degenerate=drop_degenerate, deduplicate=deduplicate)
    elseif map_data isa UniformMap || map_data isa SchnyderMap
        return fan_triangulate_faces(map_data.faces; drop_degenerate=drop_degenerate, deduplicate=deduplicate)
    elseif map_data isa UniformMeandricSystemMap
        return fan_triangulate_faces(map_data.faces; drop_degenerate=drop_degenerate, deduplicate=deduplicate)
    else
        return _empty_triangles()
    end
end

_normalize_group_string(k) = lowercase(strip(string(k)))

function _looks_like_fk_group_names(names)
    normalized = Set(_normalize_group_string(name) for name in names)
    (("green" in normalized) || ("generic" in normalized)) || return false
    (("red" in normalized) || ("blue" in normalized) || ("purple" in normalized) || ("orange" in normalized)) || return false
    (("navy" in normalized) || ("outer" in normalized)) && return false
    return true
end

function _metadata_model_name(metadata)
    metadata === nothing && return nothing
    try
        if haskey(metadata, "model")
            return lowercase(strip(string(metadata["model"])))
        elseif haskey(metadata, :model)
            return lowercase(strip(string(metadata[:model])))
        end
    catch
    end
    return nothing
end



function _metadata_lookup(metadata, key::AbstractString, default=nothing)
    metadata === nothing && return default
    if metadata isa AbstractDict
        if haskey(metadata, key)
            return metadata[key]
        end
        sk = Symbol(key)
        if haskey(metadata, sk)
            return metadata[sk]
        end
    end
    return default
end

function _metadata_has(metadata, key::AbstractString)
    metadata isa AbstractDict || return false
    return haskey(metadata, key) || haskey(metadata, Symbol(key))
end

function _maybe_float_string(x; digits::Int=4)
    x === nothing && return nothing
    if x isa AbstractFloat || x isa Real
        y = round(Float64(x); digits=digits)
        if isinteger(y)
            return string(Int(round(y)))
        end
        return string(y)
    end
    return strip(string(x))
end

function fk_q_from_p(p::Real)
    pf = float(p)
    pf == 1 && return Inf
    return (2pf / (1 - pf))^2
end

function model_parameter_pairs(map_data=nothing; metadata=nothing)
    pairs = Pair{String,String}[]

    faces_val = if _metadata_has(metadata, "faces")
        _metadata_lookup(metadata, "faces")
    elseif _metadata_has(metadata, "vertices") && _metadata_lookup(metadata, "model", nothing) == "mated_crt"
        nothing
    elseif map_data !== nothing
        try
            num_faces(map_data)
        catch
            nothing
        end
    else
        nothing
    end
    faces_val === nothing || push!(pairs, "faces" => string(Int(faces_val)))

    if _metadata_lookup(metadata, "model", nothing) == "mated_crt"
        vertices_val = _metadata_lookup(metadata, "vertices", map_data isa MatedCRTMap ? num_vertices(map_data) : nothing)
        vertices_val === nothing || push!(pairs, "vertices" => string(Int(vertices_val)))
    end

    seed_val = _metadata_lookup(metadata, "seed", nothing)
    seed_val === nothing || push!(pairs, "seed" => string(Int(seed_val)))

    if map_data isa FKMap || (_metadata_lookup(metadata, "model", nothing) in ("fk", "spanning_tree", "hc"))
        q_val = _metadata_lookup(metadata, "q", nothing)
        p_val = _metadata_lookup(metadata, "p", nothing)
        variant_val = _metadata_lookup(metadata, "variant", nothing)
        if q_val === nothing
            if p_val !== nothing
                q_val = fk_q_from_p(float(p_val))
            elseif variant_val == "spanning_tree"
                q_val = 0.0
            end
        end
        q_val === nothing || push!(pairs, "q" => _maybe_float_string(q_val; digits=5))
        if p_val !== nothing
            push!(pairs, "p" => _maybe_float_string(p_val; digits=5))
        end
        bmode = _metadata_lookup(metadata, "boundary_mode", nothing)
        bmode === nothing || push!(pairs, "boundary" => string(bmode))
    elseif map_data isa MatedCRTMap || (_metadata_lookup(metadata, "model", nothing) == "mated_crt")
        gamma_val = _metadata_lookup(metadata, "gamma", map_data isa MatedCRTMap ? map_data.gamma : nothing)
        topology_val = _metadata_lookup(metadata, "topology", map_data isa MatedCRTMap ? String(map_data.topology) : nothing)
        corr_val = _metadata_lookup(metadata, "correlation", map_data isa MatedCRTMap ? map_data.brownian_correlation : nothing)
        gamma_val === nothing || push!(pairs, "gamma" => _maybe_float_string(gamma_val; digits=5))
        topology_val === nothing || push!(pairs, "topology" => string(topology_val))
        corr_val === nothing || push!(pairs, "rho" => _maybe_float_string(corr_val; digits=5))
    end

    return pairs
end

function model_parameter_string(map_data=nothing; metadata=nothing)
    pairs = model_parameter_pairs(map_data; metadata=metadata)
    isempty(pairs) && return ""
    return join(["$(k)=$(v)" for (k, v) in pairs], ", ")
end

function model_display_name(map_data=nothing; metadata=nothing)
    default_model = map_data isa HalfPlaneMeandricSystemMap ? "half_plane_meandric" :
                    map_data isa UniformMeandricSystemMap ? _meandric_model_name(map_data) :
                    map_data isa MatedCRTMap ? "mated_crt" :
                    map_data isa UniformMap ? "uniform" :
                    map_data isa SchnyderMap ? "schnyder" :
                    map_data isa FKMap ? "fk" : "map"
    model = lowercase(strip(string(something(_metadata_lookup(metadata, "model", nothing), default_model))))
    if model == "fk"
        variant_val = _metadata_lookup(metadata, "variant", nothing)
        if variant_val == "spanning_tree"
            return "spanning_tree"
        end
    end
    return model
end

function default_render_title(map_data=nothing; metadata=nothing, dimension=nothing, engine=nothing)
    base = model_display_name(map_data; metadata=metadata)
    bits = String[]
    push!(bits, base)
    if dimension !== nothing
        eng = engine === nothing ? (Int(dimension) == 2 ? "Tutte" : "SFDP") : string(engine)
        push!(bits, "$(Int(dimension))D $(eng)")
    end
    params = model_parameter_string(map_data; metadata=metadata)
    isempty(params) || push!(bits, params)
    return join(bits, " · ")
end

function should_include_exploration(map_data=nothing; edge_groups=nothing, metadata=nothing)
    map_data isa FKMap && return true
    map_data isa MatedCRTMap && return false

    model_name = _metadata_model_name(metadata)
    model_name in ("fk", "spanning_tree", "hc") && return true
    model_name == "mated_crt" && return false

    if edge_groups !== nothing
        try
            return _looks_like_fk_group_names(keys(edge_groups))
        catch
            return false
        end
    end

    return false
end

function _sanitize_edge_group_dict(edge_groups::Dict{String,<:Any}; drop_loops::Bool=true)
    out = Dict{String,Matrix{Int32}}()
    for (key, value) in edge_groups
        arr = sanitize_edge_array(value)
        if drop_loops && size(arr, 1) > 0
            mask = arr[:, 1] .!= arr[:, 2]
            arr = arr[mask, :]
        end
        size(arr, 1) == 0 && continue
        if haskey(out, key)
            out[key] = vcat(out[key], arr)
        else
            out[key] = arr
        end
    end
    return out
end

function _meandric_component_labels(map_data::HalfPlaneMeandricSystemMap)
    labels, component_sizes = _meandric_component_labels_and_sizes(map_data.upper_adj, map_data.lower_adj)
    excluded_labels = _component_labels_touching_vertices(labels, map_data.boundary_vertices)
    return labels, component_sizes, excluded_labels
end

function _meandric_component_labels(map_data::UniformMeandricSystemMap)
    labels, component_sizes = _meandric_component_labels_and_sizes(map_data.upper_adj, map_data.lower_adj)
    return labels, component_sizes, Set{Int32}()
end

function _append_meandric_highlights!(out::Dict{String,Matrix{Int32}}, map_data::Union{HalfPlaneMeandricSystemMap,UniformMeandricSystemMap})
    labels, component_sizes, excluded_labels = _meandric_component_labels(map_data)
    ranked_labels = _sorted_meandric_component_labels(component_sizes; excluded_labels=excluded_labels)
    isempty(ranked_labels) && return out

    top_count = min(length(ranked_labels), length(MEANDRIC_HIGHLIGHT_GROUPS))
    root_to_group = Dict{Int32,String}()
    highlighted = Dict{String,Vector{NTuple{2,Int32}}}()

    for i in 1:top_count
        group = MEANDRIC_HIGHLIGHT_GROUPS[i]
        root_to_group[ranked_labels[i]] = group
        highlighted[group] = NTuple{2,Int32}[]
    end

    if length(ranked_labels) > top_count
        for rank_idx in (top_count + 1):length(ranked_labels)
            group = _meandric_fade_group_name(rank_idx)
            root_to_group[ranked_labels[rank_idx]] = group
            highlighted[group] = NTuple{2,Int32}[]
        end
    end

    for i in eachindex(map_data.edge_u)
        group = _normalize_group_string(map_data.edge_group[i])
        group in ("upper", "lower") || continue
        u = map_data.edge_u[i]
        v = map_data.edge_v[i]
        labels[Int(u) + 1] == labels[Int(v) + 1] || continue
        highlight_group = get(root_to_group, labels[Int(u) + 1], nothing)
        highlight_group === nothing && continue
        push!(highlighted[highlight_group], (u, v))
    end

    for (name, edges_raw) in highlighted
        isempty(edges_raw) && continue
        arr = Matrix{Int32}(undef, length(edges_raw), 2)
        for (i, (u, v)) in enumerate(edges_raw)
            arr[i, 1] = u
            arr[i, 2] = v
        end
        out[name] = arr
    end
    return out
end

function grouped_edges(map_data=nothing; edge_groups=nothing, drop_loops::Bool=true)
    if edge_groups !== nothing
        fk_like = should_include_exploration(map_data; edge_groups=edge_groups)
        d = Dict{String,Matrix{Int32}}()
        for (k, v) in pairs(edge_groups)
            name = _normalize_group_string(k)
            if fk_like && name == "green"
                name = "generic"
            end
            arr = sanitize_edge_array(v)
            if haskey(d, name)
                d[name] = vcat(d[name], arr)
            else
                d[name] = arr
            end
        end
        out = _sanitize_edge_group_dict(d; drop_loops=drop_loops)
        if map_data isa HalfPlaneMeandricSystemMap || map_data isa UniformMeandricSystemMap
            _append_meandric_highlights!(out, map_data)
            out = _sanitize_edge_group_dict(out; drop_loops=drop_loops)
        end
        return out
    end

    if map_data === nothing
        return Dict{String,Matrix{Int32}}()
    elseif map_data isa FKMap
        out = Dict{String,Matrix{Int32}}()
        raw = raw_active_edges_by_color(map_data; drop_loops=drop_loops)
        for color in (GREEN, RED, BLUE, PURPLE, ORANGE)
            arr = get(raw, Int(color), _empty_edges())
            size(arr, 1) == 0 && continue
            name = FK_COLOR_NAME[Int(color)]
            if name == "green"
                name = "generic"
            end
            out[name] = arr
        end
        return out
    elseif map_data isa HalfPlaneMeandricSystemMap || map_data isa UniformMeandricSystemMap
        temp = Dict{String,Vector{NTuple{2,Int32}}}()
        for i in eachindex(map_data.edge_u)
            u = map_data.edge_u[i]
            v = map_data.edge_v[i]
            if drop_loops && u == v
                continue
            end
            group = _normalize_group_string(map_data.edge_group[i])
            push!(get!(temp, group, NTuple{2,Int32}[]), (u, v))
        end
        out = Dict{String,Matrix{Int32}}()
        for (name, edges_raw) in temp
            arr = Matrix{Int32}(undef, length(edges_raw), 2)
            for (i, (u, v)) in enumerate(edges_raw)
                arr[i, 1] = u
                arr[i, 2] = v
            end
            out[name] = arr
        end
        _append_meandric_highlights!(out, map_data)
        return _sanitize_edge_group_dict(out; drop_loops=drop_loops)
    elseif map_data isa MatedCRTMap
        temp = Dict{String,Vector{NTuple{2,Int32}}}()
        for i in eachindex(map_data.edge_u)
            u = map_data.edge_u[i]
            v = map_data.edge_v[i]
            if drop_loops && u == v
                continue
            end
            push!(get!(temp, "generic", NTuple{2,Int32}[]), (u, v))
            kind = map_data.edge_kind[i]
            if kind == :upper
                push!(get!(temp, "red", NTuple{2,Int32}[]), (u, v))
            elseif kind == :lower
                push!(get!(temp, "blue", NTuple{2,Int32}[]), (u, v))
            end
        end
        out = Dict{String,Matrix{Int32}}()
        for (name, pairs_raw) in temp
            arr = Matrix{Int32}(undef, length(pairs_raw), 2)
            for (i, (u, v)) in enumerate(pairs_raw)
                arr[i, 1] = u
                arr[i, 2] = v
            end
            out[name] = arr
        end
        return _sanitize_edge_group_dict(out; drop_loops=drop_loops)
    elseif map_data isa SchnyderMap
        temp = Dict{String,Vector{NTuple{2,Int32}}}()
        for i in eachindex(map_data.edge_u)
            u = map_data.edge_u[i]
            v = map_data.edge_v[i]
            if drop_loops && u == v
                continue
            end
            edge = (u, v)
            push!(get!(temp, "generic", NTuple{2,Int32}[]), edge)
            color = map_data.edge_color[i]
            if color in ("orange", "navy", "green", "outer")
                push!(get!(temp, color, NTuple{2,Int32}[]), edge)
            end
        end
        out = Dict{String,Matrix{Int32}}()
        for (k, vs) in temp
            arr = Matrix{Int32}(undef, length(vs), 2)
            for (i, (u, v)) in enumerate(vs)
                arr[i, 1] = u
                arr[i, 2] = v
            end
            out[k] = arr
        end
        return out
    elseif map_data isa UniformMap
        arr = layout_edges(map_data; drop_loops=drop_loops)
        return Dict("generic" => arr)
    else
        arr = layout_edges(map_data; drop_loops=drop_loops)
        return Dict("generic" => arr)
    end
end

function pad_positions_3d(pos)
    arr = Float32.(pos)
    ndims(arr) == 2 || throw(ArgumentError("positions must be a matrix"))
    if size(arr, 2) == 3
        return arr
    elseif size(arr, 2) == 2
        out = zeros(Float32, size(arr, 1), 3)
        out[:, 1:2] .= arr
        return out
    else
        throw(ArgumentError("positions must have shape (N, 2) or (N, 3)"))
    end
end

function subset_sphere_circle_geometry(geometry, keep_indices)
    geometry === nothing && return nothing
    idx = Int.(collect(keep_indices))
    out = Dict{String,Any}(
        "centers" => geometry["centers"][idx, :],
        "normals" => geometry["normals"][idx, :],
        "radii" => geometry["radii"][idx],
        "offsets" => geometry["offsets"][idx],
        "sphere_radius" => geometry["sphere_radius"],
        "projection_scale" => get(geometry, "projection_scale", 1.0),
        "projection" => get(geometry, "projection", "inverse_stereographic"),
    )
    if haskey(geometry, "outer_geometry")
        out["outer_geometry"] = geometry["outer_geometry"]
    end
    return out
end

function scale_sphere_circle_geometry(geometry, scale::Real)
    geometry === nothing && return nothing
    sf = float(scale)
    return Dict{String,Any}(
        "centers" => sf .* Float64.(geometry["centers"]),
        "normals" => Float64.(geometry["normals"]),
        "radii" => sf .* Float64.(geometry["radii"]),
        "offsets" => Float64.(geometry["offsets"]),
        "sphere_radius" => sf * float(geometry["sphere_radius"]),
        "projection_scale" => get(geometry, "projection_scale", 1.0),
        "projection" => get(geometry, "projection", "inverse_stereographic"),
    )
end

function sphere_circle_segment_points(geometry; segments::Integer=128)
    geometry === nothing && return zeros(Float32, 0, 3)
    seg_count = Int(segments)
    seg_count >= 8 || throw(ArgumentError("segments must be at least 8"))

    centers = Float64.(geometry["centers"])
    normals = Float64.(geometry["normals"])
    radii = Float64.(geometry["radii"])
    size(centers, 1) == size(normals, 1) == length(radii) || throw(ArgumentError("sphere circle geometry arrays must align"))

    rows = Float32[]
    sizehint!(rows, 6 * seg_count * size(centers, 1))

    for i in 1:size(centers, 1)
        radius = radii[i]
        radius > 0.0 || continue
        isfinite(radius) || continue

        normal = normals[i, :]
        nlen = norm(normal)
        nlen > 0.0 || continue
        n = normal ./ nlen
        center = centers[i, :]

        ref = abs(n[3]) < 0.9 ? [0.0, 0.0, 1.0] : [1.0, 0.0, 0.0]
        u = cross(n, ref)
        ulen = norm(u)
        ulen > 0.0 || continue
        u ./= ulen
        v = cross(n, u)

        prev = center .+ radius .* u
        for j in 1:seg_count
            θ = 2π * j / seg_count
            current = center .+ radius .* (cos(θ) .* u .+ sin(θ) .* v)
            append!(rows, Float32.(prev))
            append!(rows, Float32.(current))
            prev = current
        end
    end

    isempty(rows) && return zeros(Float32, 0, 3)
    return Matrix{Float32}(permutedims(reshape(rows, 3, :)))
end

@inline _edge_key(u::Integer, v::Integer) = Int(u) <= Int(v) ? (Int32(u), Int32(v)) : (Int32(v), Int32(u))

@inline function _row_is_finite(pos::AbstractMatrix, idx::Integer)
    @inbounds for j in axes(pos, 2)
        isfinite(pos[idx, j]) || return false
    end
    return true
end

function _triangle_edge_midpoint_segments(pos::AbstractMatrix, triangles::AbstractMatrix, generic_edges::AbstractMatrix)
    size(pos, 2) in (2, 3) || throw(ArgumentError("positions must have shape (N, 2) or (N, 3)"))
    tri = sanitize_triangles(triangles; drop_degenerate=true, deduplicate=false)
    size(tri, 1) == 0 && return zeros(Float32, 0, size(pos, 2))
    gen = sanitize_edge_array(generic_edges)
    size(gen, 1) == 0 && return zeros(Float32, 0, size(pos, 2))

    generic_set = Set{NTuple{2,Int32}}()
    sizehint!(generic_set, size(gen, 1))
    for i in 1:size(gen, 1)
        push!(generic_set, _edge_key(gen[i, 1], gen[i, 2]))
    end

    posf = Float32.(pos)
    d = size(posf, 2)
    rows = Float32[]
    sizehint!(rows, 2 * d * size(tri, 1))

    @inbounds for i in 1:size(tri, 1)
        a = tri[i, 1]
        b = tri[i, 2]
        c = tri[i, 3]
        keys = (_edge_key(a, b), _edge_key(b, c), _edge_key(c, a))
        flags = (keys[1] in generic_set, keys[2] in generic_set, keys[3] in generic_set)
        nflags = Int(flags[1]) + Int(flags[2]) + Int(flags[3])
        nflags == 2 || continue

        selected = NTuple{2,Int32}[]
        flags[1] && push!(selected, keys[1])
        flags[2] && push!(selected, keys[2])
        flags[3] && push!(selected, keys[3])
        length(selected) == 2 || continue

        (u0, v0) = selected[1]
        (u1, v1) = selected[2]
        iu0 = Int(u0) + 1
        iv0 = Int(v0) + 1
        iu1 = Int(u1) + 1
        iv1 = Int(v1) + 1

        (1 <= iu0 <= size(posf, 1) && 1 <= iv0 <= size(posf, 1) && 1 <= iu1 <= size(posf, 1) && 1 <= iv1 <= size(posf, 1)) || continue
        (_row_is_finite(posf, iu0) && _row_is_finite(posf, iv0) && _row_is_finite(posf, iu1) && _row_is_finite(posf, iv1)) || continue

        for j in 1:d
            push!(rows, 0.5f0 * (posf[iu0, j] + posf[iv0, j]))
        end
        for j in 1:d
            push!(rows, 0.5f0 * (posf[iu1, j] + posf[iv1, j]))
        end
    end

    isempty(rows) && return zeros(Float32, 0, d)
    return Matrix{Float32}(permutedims(reshape(rows, d, :)))
end

function _contracted_exploration_loop_segments_2d(pos::AbstractMatrix, loop_edges)
    loops = sanitize_edge_array(loop_edges)
    size(loops, 1) == 0 && return zeros(Float32, 0, 2)
    posf = Float32.(pos)
    rows = Float32[]

    @inbounds for i in 1:size(loops, 1)
        u = Int(loops[i, 1]) + 1
        v = Int(loops[i, 2]) + 1
        (1 <= u <= size(posf, 1) && 1 <= v <= size(posf, 1)) || continue
        (_row_is_finite(posf, u) && _row_is_finite(posf, v)) || continue

        pu = Float64.(@view posf[u, :])
        pv = Float64.(@view posf[v, :])
        mid = 0.5 .* (pu .+ pv)
        dir = pv .- pu
        len = norm(dir)
        len > 0.0 || continue
        tangent = dir ./ len
        normal = Float64[-tangent[2], tangent[1]]
        major = 0.18 * len
        minor = 0.12 * len

        prev = nothing
        for j in 0:20
            θ = 2π * j / 20
            point = mid .+ major * cos(θ) .* tangent .+ minor * sin(θ) .* normal
            if prev !== nothing
                append!(rows, Float32.(prev))
                append!(rows, Float32.(point))
            end
            prev = point
        end
    end

    isempty(rows) && return zeros(Float32, 0, 2)
    return Matrix{Float32}(permutedims(reshape(rows, 2, :)))
end

function _infer_sphere_radius_from_positions(pos::AbstractMatrix)
    size(pos, 2) == 3 || throw(ArgumentError("sphere radius inference expects 3D positions"))
    radii = Float64[]
    sizehint!(radii, size(pos, 1))
    @inbounds for i in 1:size(pos, 1)
        _row_is_finite(pos, i) || continue
        r = norm(Float64.(view(pos, i, :)))
        (isfinite(r) && r > 0.0) || continue
        push!(radii, r)
    end
    isempty(radii) && return 1.0
    return sum(radii) / length(radii)
end

function _triangle_edge_spherical_midpoint_segments(
    pos::AbstractMatrix,
    triangles::AbstractMatrix,
    generic_edges::AbstractMatrix;
    sphere_radius::Real=1.0,
)
    size(pos, 2) == 3 || throw(ArgumentError("spherical exploration segments require 3D positions"))
    tri = sanitize_triangles(triangles; drop_degenerate=true, deduplicate=false)
    size(tri, 1) == 0 && return zeros(Float32, 0, 3)
    gen = sanitize_edge_array(generic_edges)
    size(gen, 1) == 0 && return zeros(Float32, 0, 3)

    generic_set = Set{NTuple{2,Int32}}()
    sizehint!(generic_set, size(gen, 1))
    for i in 1:size(gen, 1)
        push!(generic_set, _edge_key(gen[i, 1], gen[i, 2]))
    end

    radius = float(sphere_radius)
    radius > 0.0 || throw(ArgumentError("sphere radius must be positive"))
    posf = Float32.(pos)
    rows = Float32[]
    sizehint!(rows, 6 * size(tri, 1))

    @inline function _sphere_midpoint(iu::Int, iv::Int)
        vec = Float64.(@view(posf[iu, :])) .+ Float64.(@view(posf[iv, :]))
        nrm = norm(vec)
        nrm > 0.0 || return nothing
        return Float32.(radius .* (vec ./ nrm))
    end

    @inbounds for i in 1:size(tri, 1)
        a = tri[i, 1]
        b = tri[i, 2]
        c = tri[i, 3]
        keys = (_edge_key(a, b), _edge_key(b, c), _edge_key(c, a))
        flags = (keys[1] in generic_set, keys[2] in generic_set, keys[3] in generic_set)
        nflags = Int(flags[1]) + Int(flags[2]) + Int(flags[3])
        nflags == 2 || continue

        selected = NTuple{2,Int32}[]
        flags[1] && push!(selected, keys[1])
        flags[2] && push!(selected, keys[2])
        flags[3] && push!(selected, keys[3])
        length(selected) == 2 || continue

        (u0, v0) = selected[1]
        (u1, v1) = selected[2]
        iu0 = Int(u0) + 1
        iv0 = Int(v0) + 1
        iu1 = Int(u1) + 1
        iv1 = Int(v1) + 1

        (1 <= iu0 <= size(posf, 1) && 1 <= iv0 <= size(posf, 1) && 1 <= iu1 <= size(posf, 1) && 1 <= iv1 <= size(posf, 1)) || continue
        (_row_is_finite(posf, iu0) && _row_is_finite(posf, iv0) && _row_is_finite(posf, iu1) && _row_is_finite(posf, iv1)) || continue

        mid0 = _sphere_midpoint(iu0, iv0)
        mid1 = _sphere_midpoint(iu1, iv1)
        (mid0 === nothing || mid1 === nothing) && continue
        append!(rows, mid0)
        append!(rows, mid1)
    end

    isempty(rows) && return zeros(Float32, 0, 3)
    return Matrix{Float32}(permutedims(reshape(rows, 3, :)))
end

function _contracted_exploration_loop_segments_sphere(pos::AbstractMatrix, loop_edges; sphere_radius::Real=1.0)
    loops = sanitize_edge_array(loop_edges)
    size(loops, 1) == 0 && return zeros(Float32, 0, 3)
    posf = Float32.(pos)
    radius = float(sphere_radius)
    rows = Float32[]

    @inbounds for i in 1:size(loops, 1)
        u = Int(loops[i, 1]) + 1
        v = Int(loops[i, 2]) + 1
        (1 <= u <= size(posf, 1) && 1 <= v <= size(posf, 1)) || continue
        (_row_is_finite(posf, u) && _row_is_finite(posf, v)) || continue

        pu = Float64.(@view posf[u, :])
        pv = Float64.(@view posf[v, :])
        mid = pu .+ pv
        nrm_mid = norm(mid)
        nrm_mid > 0.0 || continue
        center_dir = mid ./ nrm_mid

        edge_dir = pv .- pu
        edge_dir .-= dot(edge_dir, center_dir) .* center_dir
        edge_norm = norm(edge_dir)
        edge_norm > 0.0 || continue
        tangent = edge_dir ./ edge_norm
        bitangent = cross(center_dir, tangent)
        bitangent ./= norm(bitangent)

        chord = norm(pv .- pu)
        α = clamp(0.22 * chord / max(radius, 1.0e-8), 0.05, 0.18)

        prev = nothing
        for j in 0:24
            θ = 2π * j / 24
            vec = cos(α) .* center_dir .+ sin(α) .* (cos(θ) .* tangent .+ sin(θ) .* bitangent)
            point = Float32.(radius .* (vec ./ norm(vec)))
            if prev !== nothing
                append!(rows, prev)
                append!(rows, point)
            end
            prev = point
        end
    end

    isempty(rows) && return zeros(Float32, 0, 3)
    return Matrix{Float32}(permutedims(reshape(rows, 3, :)))
end

function _exploration_green_edge_pairs(map_data, metadata)
    pairs = try
        _metadata_lookup(metadata, "exploration_green_edge_pairs", nothing)
    catch
        nothing
    end
    if pairs !== nothing
        return sanitize_edge_array(pairs)
    elseif map_data isa FKMap
        return sanitize_edge_array(map_data.triangle_green_edges)
    end
    return Matrix{Int32}(undef, 0, 2)
end

function _exploration_vertex_maps(pos::AbstractMatrix, metadata)
    current_to_base = Int32.(collect(0:(size(pos, 1) - 1)))
    kept = try
        _metadata_lookup(metadata, "gasket_kept_vertex_indices", nothing)
    catch
        nothing
    end
    if kept !== nothing
        kept_vec = Int32.(collect(kept))
        if length(kept_vec) == size(pos, 1)
            current_to_base = kept_vec
        end
    end

    base_to_current = Dict{Int32,Int32}()
    for (i, base_vertex) in enumerate(current_to_base)
        base_to_current[base_vertex] = Int32(i - 1)
    end
    return current_to_base, base_to_current
end

function _exploration_pair_point_map(pos::AbstractMatrix, generic_edges::AbstractMatrix, metadata)
    d = size(pos, 2)
    current_to_base, _ = _exploration_vertex_maps(pos, metadata)
    posf = Float32.(pos)
    points = Dict{Tuple{Int32,Int32},Vector{Float32}}()

    sphere_mode = d == 3 && _metadata_lookup(metadata, "sphere_projection", nothing) == "inverse_stereographic"
    sphere_radius = sphere_mode ? float(something(
        try
            _metadata_lookup(metadata, "sphere_radius", nothing)
        catch
            nothing
        end,
        _infer_sphere_radius_from_positions(posf),
    )) : 1.0

    @inbounds for i in 1:size(generic_edges, 1)
        u = Int(generic_edges[i, 1]) + 1
        v = Int(generic_edges[i, 2]) + 1
        (1 <= u <= size(posf, 1) && 1 <= v <= size(posf, 1)) || continue
        (_row_is_finite(posf, u) && _row_is_finite(posf, v)) || continue

        base_u = current_to_base[u]
        base_v = current_to_base[v]
        key = _edge_key(base_u, base_v)

        point = if sphere_mode
            vec = Float64.(@view(posf[u, :])) .+ Float64.(@view(posf[v, :]))
            nrm = norm(vec)
            nrm > 0.0 || continue
            Float32.(sphere_radius .* (vec ./ nrm))
        else
            0.5f0 .* (Float32.(@view(posf[u, :])) .+ Float32.(@view(posf[v, :])))
        end
        points[key] = collect(point)
    end

    return points
end

function _exploration_walk_components(edge_pairs::AbstractMatrix)
    pairs = sanitize_edge_array(edge_pairs)
    size(pairs, 1) == 0 && return Tuple{Vector{Int32},Bool}[]

    incident = Dict{Int32,Vector{Int}}()
    for i in 1:size(pairs, 1)
        a = pairs[i, 1]
        b = pairs[i, 2]
        a == b && continue
        push!(get!(incident, a, Int[]), i)
        push!(get!(incident, b, Int[]), i)
    end

    used = falses(size(pairs, 1))
    walks = Tuple{Vector{Int32},Bool}[]

    function other_endpoint(edge_idx::Int, node::Int32)
        a = pairs[edge_idx, 1]
        b = pairs[edge_idx, 2]
        return a == node ? b : a
    end

    for seed_edge in 1:size(pairs, 1)
        used[seed_edge] && continue
        a = pairs[seed_edge, 1]
        b = pairs[seed_edge, 2]
        a == b && continue

        component_edges = Int[]
        component_nodes = Set{Int32}()
        queue = Int[seed_edge]
        while !isempty(queue)
            edge_idx = popfirst!(queue)
            used[edge_idx] && continue
            push!(component_edges, edge_idx)
            used[edge_idx] = true
            u = pairs[edge_idx, 1]
            v = pairs[edge_idx, 2]
            push!(component_nodes, u)
            push!(component_nodes, v)
            for node in (u, v)
                for next_edge in get(incident, node, Int[])
                    used[next_edge] || push!(queue, next_edge)
                end
            end
        end

        component_used = Set(component_edges)
        degree = Dict{Int32,Int}(node => 0 for node in component_nodes)
        for edge_idx in component_edges
            degree[pairs[edge_idx, 1]] += 1
            degree[pairs[edge_idx, 2]] += 1
        end

        start = let degree_one = sort!(Int32[node for (node, deg) in degree if deg == 1])
            isempty(degree_one) ? minimum(collect(component_nodes)) : degree_one[1]
        end

        walk = Int32[start]
        current = start
        previous_edge = 0
        while !isempty(component_used)
            next_edge = 0
            for edge_idx in get(incident, current, Int[])
                edge_idx in component_used || continue
                if previous_edge == 0 || edge_idx != previous_edge || length(component_used) == 1
                    next_edge = edge_idx
                    break
                end
            end
            next_edge == 0 && break
            delete!(component_used, next_edge)
            current = other_endpoint(next_edge, current)
            push!(walk, current)
            previous_edge = next_edge
        end

        closed = length(walk) >= 2 && walk[end] == walk[1]
        push!(walks, (walk, closed))
    end

    return walks
end

function _fk_exploration_segments(map_data::FKMap, pos::AbstractMatrix, generic_edges::AbstractMatrix, metadata)
    green_pairs = _exploration_green_edge_pairs(map_data, metadata)
    size(green_pairs, 1) == 0 && return zeros(Float32, 0, size(pos, 2))
    generic = sanitize_edge_array(generic_edges)
    size(generic, 1) == 0 && return zeros(Float32, 0, size(pos, 2))

    point_map = _exploration_pair_point_map(pos, generic, metadata)
    isempty(point_map) && return zeros(Float32, 0, size(pos, 2))

    rep_map = Dict{Int32,Tuple{Int32,Int32}}()
    for edge_id in unique(Int32[v for v in vec(green_pairs) if v >= 0])
        edge_idx = Int(edge_id) + 1
        edge_idx <= length(map_data.edge_u) || continue
        key = _edge_key(map_data.edge_u[edge_idx], map_data.edge_v[edge_idx])
        haskey(point_map, key) || continue
        rep_map[edge_id] = key
    end

    rows = Float32[]
    d = size(pos, 2)
    for (walk, closed) in _exploration_walk_components(green_pairs)
        mapped = Tuple{Int32,Int32}[]
        for edge_id in walk
            rep = get(rep_map, edge_id, nothing)
            rep === nothing && continue
            isempty(mapped) || mapped[end] != rep || continue
            push!(mapped, rep)
        end
        closed && length(mapped) >= 2 && mapped[1] == mapped[end] && pop!(mapped)
        length(mapped) >= 2 || continue

        for i in 1:(length(mapped) - 1)
            pa = point_map[mapped[i]]
            pb = point_map[mapped[i + 1]]
            append!(rows, pa)
            append!(rows, pb)
        end
        if closed && length(mapped) >= 3
            pa = point_map[mapped[end]]
            pb = point_map[mapped[1]]
            append!(rows, pa)
            append!(rows, pb)
        end
    end

    isempty(rows) && return zeros(Float32, 0, d)
    return Matrix{Float32}(permutedims(reshape(rows, d, :)))
end

function exploration_segment_points(map_data=nothing, pos=nothing; edge_groups=nothing, faces=nothing, triangles=nothing, metadata=nothing)
    pos === nothing && return zeros(Float32, 0, 2)
    pos_arr = Float32.(pos)
    ndims(pos_arr) == 2 || throw(ArgumentError("positions must be a matrix"))
    d = size(pos_arr, 2)
    d in (2, 3) || throw(ArgumentError("positions must have shape (N, 2) or (N, 3)"))

    should_include_exploration(map_data; edge_groups=edge_groups, metadata=metadata) || return zeros(Float32, 0, d)

    groups = grouped_edges(map_data; edge_groups=edge_groups, drop_loops=true)
    generic = get(groups, "generic", _empty_edges())
    size(generic, 1) == 0 && return zeros(Float32, 0, d)

    if map_data isa FKMap || _metadata_has(metadata, "exploration_green_edge_pairs")
        return (map_data isa FKMap) ? _fk_exploration_segments(map_data, pos_arr, generic, metadata) : zeros(Float32, 0, d)
    end

    tri = surface_triangles(map_data; faces=faces, triangles=triangles, drop_degenerate=true, deduplicate=true)
    size(tri, 1) == 0 && return zeros(Float32, 0, d)

    if d == 3 && _metadata_lookup(metadata, "sphere_projection", nothing) == "inverse_stereographic"
        sphere_radius = something(
            try
                _metadata_lookup(metadata, "sphere_radius", nothing)
            catch
                nothing
            end,
            _infer_sphere_radius_from_positions(pos_arr),
        )
        return _triangle_edge_spherical_midpoint_segments(pos_arr, tri, generic; sphere_radius=float(sphere_radius))
    end

    segs = _triangle_edge_midpoint_segments(pos_arr, tri, generic)
    return segs
end

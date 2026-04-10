# Shared helpers for converting combinatorial maps into drawable geometry.

const EDGE_COLOR_HINTS = Dict(
    "generic" => "#b0b7bf",
    "green" => "#1f7a3a",
    "red" => "#c0392b",
    "blue" => "#2980b9",
    "purple" => "#8e44ad",
    "orange" => "#d35400",
    "outer" => "#2c3e50",
    "navy" => "#2c3e50",
    "exploration" => "#27ae60",
)

_empty_triangles() = Matrix{Int32}(undef, 0, 3)
_empty_edges() = Matrix{Int32}(undef, 0, 2)

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
        for i in 2:(length(face) - 1)
            push!(triangles, (anchor, face[i], face[i + 1]))
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
    elseif map_data isa UniformMap || map_data isa SchnyderMap
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
    end

    return pairs
end

function model_parameter_string(map_data=nothing; metadata=nothing)
    pairs = model_parameter_pairs(map_data; metadata=metadata)
    isempty(pairs) && return ""
    return join(["$(k)=$(v)" for (k, v) in pairs], ", ")
end

function model_display_name(map_data=nothing; metadata=nothing)
    model = lowercase(strip(string(something(_metadata_lookup(metadata, "model", nothing), map_data isa UniformMap ? "uniform" : map_data isa SchnyderMap ? "schnyder" : map_data isa FKMap ? "fk" : "map"))))
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

    model_name = _metadata_model_name(metadata)
    model_name in ("fk", "spanning_tree", "hc") && return true

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
        return _sanitize_edge_group_dict(d; drop_loops=drop_loops)
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

function exploration_segment_points(map_data=nothing, pos=nothing; edge_groups=nothing, faces=nothing, triangles=nothing, metadata=nothing)
    pos === nothing && return zeros(Float32, 0, 2)
    pos_arr = Float32.(pos)
    ndims(pos_arr) == 2 || throw(ArgumentError("positions must be a matrix"))
    d = size(pos_arr, 2)
    d in (2, 3) || throw(ArgumentError("positions must have shape (N, 2) or (N, 3)"))

    should_include_exploration(map_data; edge_groups=edge_groups, metadata=metadata) || return zeros(Float32, 0, d)

    tri = surface_triangles(map_data; faces=faces, triangles=triangles, drop_degenerate=true, deduplicate=true)
    size(tri, 1) == 0 && return zeros(Float32, 0, d)

    groups = grouped_edges(map_data; edge_groups=edge_groups, drop_loops=true)
    generic = get(groups, "generic", _empty_edges())
    size(generic, 1) == 0 && return zeros(Float32, 0, d)

    return _triangle_edge_midpoint_segments(pos_arr, tri, generic)
end

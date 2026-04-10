# Model-specific preparation of graph layout problems.

Base.@kwdef struct LayoutProblem
    num_vertices::Int
    edges::Matrix{Int32}
    edge_groups::Dict{String,Matrix{Int32}}
    boundary_vertices::Union{Nothing,Vector{Int32}} = nothing
    boundary_positions::Union{Nothing,Matrix{Float64}} = nothing
    faces::Any = nothing
    surface_triangles::Union{Nothing,Matrix{Int32}} = nothing
    metadata::Dict{String,Any} = Dict{String,Any}()
    render_map_data::Any = nothing
end

function _regular_polygon(num_vertices::Integer, radius::Real; phase::Real=0.0)
    n = Int(num_vertices)
    n <= 0 && return zeros(Float64, 0, 2)
    angles = range(float(phase), float(phase) + 2π; length=n + 1)[1:end-1]
    return hcat(float(radius) .* cos.(angles), float(radius) .* sin.(angles))
end

equilateral_triangle(radius::Real) = _regular_polygon(3, radius; phase=π / 2)
square(radius::Real) = Float64[-radius -radius; radius -radius; radius radius; -radius radius]

function _edge_groups_for_map(map_data; drop_loops::Bool=true)
    return grouped_edges(map_data; drop_loops=drop_loops)
end

_faces_as_lists(faces) = _iter_faces(faces)

_deduplicate_undirected_edges(edges) = first(collapse_undirected_edges(edges; drop_loops=true))

function _ordered_unique_vertices(chunks...)
    out = Int32[]
    seen = Set{Int32}()
    for chunk in chunks
        for raw in Int32.(collect(chunk))
            if raw < 0 || raw in seen
                continue
            end
            push!(seen, raw)
            push!(out, raw)
        end
    end
    return out
end

function _reindex_edges(edge_array, old_to_new)
    arr = sanitize_edge_array(edge_array)
    if size(arr, 1) == 0
        return Matrix{Int32}(undef, 0, 2)
    end
    remapped = Matrix{Int32}(undef, size(arr, 1), 2)
    for i in 1:size(arr, 1)
        remapped[i, 1] = Int32(old_to_new[Int(arr[i, 1])])
        remapped[i, 2] = Int32(old_to_new[Int(arr[i, 2])])
    end
    return _deduplicate_undirected_edges(remapped)
end

function _reindex_triangles(triangles, old_to_new)
    tri = Int32.(triangles)
    size(tri, 1) == 0 && return Matrix{Int32}(undef, 0, 3)
    remapped = Matrix{Int32}(undef, size(tri, 1), 3)
    for i in 1:size(tri, 1)
        for j in 1:3
            remapped[i, j] = Int32(old_to_new[Int(tri[i, j])])
        end
    end
    return sanitize_triangles(remapped; drop_degenerate=true, deduplicate=true)
end

function _first_simple_face(faces, target_size::Integer)
    best = nothing
    best_unique = -1
    for (idx1, face) in enumerate(_faces_as_lists(faces))
        unique_in_order = Int32[]
        seen = Set{Int32}()
        for v in face
            if v in seen
                continue
            end
            push!(seen, v)
            push!(unique_in_order, v)
        end
        if length(unique_in_order) == Int(target_size)
            return (idx1 - 1, unique_in_order)
        end
        if length(unique_in_order) > best_unique
            best = (idx1 - 1, unique_in_order)
            best_unique = length(unique_in_order)
        end
    end
    return best
end

function _prepare_schnyder_problem(map_data::SchnyderMap; dimension::Int, boundary_scale::Float64)
    edges = layout_edges(map_data; drop_loops=true)
    edge_groups = _edge_groups_for_map(map_data)
    faces = _faces_as_lists(map_data.faces)

    if dimension == 3
        return LayoutProblem(
            num_vertices=map_data.nverts,
            edges=edges,
            edge_groups=edge_groups,
            faces=faces,
            metadata=Dict("model" => "schnyder", "dimension" => 3),
            render_map_data=map_data,
        )
    end

    outer_edges = get(edge_groups, "outer", Matrix{Int32}(undef, 0, 2))
    outer_vertices = size(outer_edges, 1) == 0 ? Int32[0, 1, 2] : sort(unique(vec(outer_edges)))
    length(outer_vertices) == 3 || throw(ArgumentError("Schnyder 2D layout expects exactly three outer vertices"))

    outer_set = Set(outer_vertices)
    interior_faces = [face for face in faces if Set(face) != outer_set]

    return LayoutProblem(
        num_vertices=map_data.nverts,
        edges=edges,
        edge_groups=edge_groups,
        boundary_vertices=outer_vertices,
        boundary_positions=equilateral_triangle(boundary_scale),
        faces=interior_faces,
        metadata=Dict("model" => "schnyder", "dimension" => 2, "outer_face_removed" => true),
        render_map_data=map_data,
    )
end

function _prepare_uniform_problem(map_data::UniformMap; dimension::Int, boundary_scale::Float64, outer_face_index=nothing)
    edges = layout_edges(map_data; drop_loops=true)
    edge_groups = _edge_groups_for_map(map_data)
    faces = _faces_as_lists(map_data.faces)

    if dimension == 3
        return LayoutProblem(
            num_vertices=num_vertices(map_data),
            edges=edges,
            edge_groups=edge_groups,
            faces=faces,
            metadata=Dict("model" => "uniform", "dimension" => 3),
            render_map_data=map_data,
        )
    end

    choice = if outer_face_index === nothing
        _first_simple_face(faces, 4)
    else
        idx = Int(outer_face_index)
        0 <= idx < length(faces) || throw(ArgumentError("outer_face_index is out of range"))
        (idx, unique(Int32.(faces[idx+1])))
    end

    choice === nothing && throw(ArgumentError("failed to find a quadrilateral face with four distinct boundary vertices for the 2D layout"))
    idx, boundary = choice
    length(boundary) >= 4 || throw(ArgumentError("failed to find a quadrilateral face with four distinct boundary vertices for the 2D layout"))

    boundary_vertices = Int32.(boundary[1:4])
    interior_faces = [face for (j1, face) in enumerate(faces) if (j1 - 1) != idx]

    return LayoutProblem(
        num_vertices=num_vertices(map_data),
        edges=edges,
        edge_groups=edge_groups,
        boundary_vertices=boundary_vertices,
        boundary_positions=square(boundary_scale),
        faces=interior_faces,
        metadata=Dict(
            "model" => "uniform",
            "dimension" => 2,
            "outer_face_removed" => true,
            "outer_face_index" => idx,
        ),
        render_map_data=map_data,
    )
end

_hc_variant(map_data::FKMap) = variant(map_data)

function _hc_gasket_selection(map_data::FKMap, gasket::AbstractString)
    gasket_norm = lowercase(strip(string(gasket)))
    gasket_norm in ("h", "c", "h_gasket", "c_gasket") || throw(ArgumentError("invalid gasket"))
    kind = gasket_norm[1]
    boundary_order = kind == 'h' ? "H" : "C"
    boundary_color = kind == 'h' ? "red" : "blue"
    fresh_color = kind == 'h' ? "orange" : "purple"
    prod_idx, start_step, end_step = maximal_gasket_span(map_data, string(kind))
    return string(kind), Int(prod_idx), Int(start_step), Int(end_step), boundary_order, boundary_color, fresh_color
end

function _hc_gasket_boundary(map_data::FKMap, gasket::AbstractString)
    kind, prod_idx, start_step, end_step, boundary_order, boundary_color, fresh_color = _hc_gasket_selection(map_data, gasket)

    boundary_vertices = Int32[map_data.opposite_vertex_at_production[prod_idx+1]]
    boundary_sources = Int32[prod_idx]

    order = sortperm(map_data.order_steps)
    for j in order
        order_step = Int(map_data.order_steps[j])
        if order_step <= start_step || order_step > end_step
            continue
        end
        if map_data.order_symbols[j] != boundary_order[1]
            continue
        end
        if Int(map_data.production_steps[j]) >= start_step
            continue
        end
        push!(boundary_vertices, map_data.vertex_of_prod[j])
        push!(boundary_sources, Int32(j - 1))
    end
    push!(boundary_vertices, map_data.opposite_vertex_at_order[prod_idx+1])
    push!(boundary_sources, Int32(prod_idx))

    unique_boundary = Int32[]
    unique_sources = Int32[]
    seen = Set{Int32}()
    for (v, src) in zip(boundary_vertices, boundary_sources)
        if v in seen
            continue
        end
        push!(seen, v)
        push!(unique_boundary, v)
        push!(unique_sources, src)
    end

    return unique_boundary, Dict{String,Any}(
        "gasket_type" => "$(kind)_gasket",
        "gasket_root_prod_index" => prod_idx,
        "gasket_span" => [start_step, end_step],
        "gasket_boundary_source_productions" => Int.(unique_sources),
        "gasket_boundary_vertex_color" => boundary_color,
        "gasket_boundary_edge_colors" => [boundary_color, fresh_color],
    )
end

function _hc_gasket_subgraph(map_data::FKMap, gasket::AbstractString)
    kind, prod_idx, start_step, end_step, _, boundary_color, fresh_color = _hc_gasket_selection(map_data, gasket)
    boundary_vertices_old, meta = _hc_gasket_boundary(map_data, gasket)

    keep_steps = [s for s in start_step:end_step if s != Int(map_data.order_face[prod_idx+1])]
    tri_faces_all = Int32.(map_data.triangulation_faces)
    tri_green_all = Int32.(map_data.triangle_green_edges)
    tri_diag_all = Int32.(map_data.triangle_diag_edge)

    selected_triangles_old = tri_faces_all[keep_steps.+1, :]
    selected_green = tri_green_all[keep_steps.+1, :]
    selected_diag = tri_diag_all[keep_steps.+1]

    edge_ids = unique(Int32[v for v in vcat(vec(selected_green), selected_diag) if v >= 0])
    edge_u = map_data.edge_u
    edge_v = map_data.edge_v
    edge_color = Int32.(map_data.edge_color)

    original_vertices = _ordered_unique_vertices(
        boundary_vertices_old,
        vec(selected_triangles_old),
        edge_u[edge_ids.+1],
        edge_v[edge_ids.+1],
    )
    old_to_new = Dict(Int(v) => i - 1 for (i, v) in enumerate(original_vertices))

    triangles = _reindex_triangles(selected_triangles_old, old_to_new)
    boundary_vertices = Int32[old_to_new[Int(v)] for v in boundary_vertices_old]

    edge_groups = Dict{String,Matrix{Int32}}()
    color_names = Dict(0 => "green", 1 => "red", 2 => "blue", 3 => "purple", 4 => "orange")
    for (color_id, name) in color_names
        mask = edge_color[edge_ids.+1] .== color_id
        any(mask) || continue
        pairs_old = hcat(edge_u[edge_ids[mask].+1], edge_v[edge_ids[mask].+1])
        pairs_new = _reindex_edges(pairs_old, old_to_new)
        if size(pairs_new, 1) > 0
            edge_groups[name] = pairs_new
        end
    end

    all_edges = if isempty(edge_groups)
        Matrix{Int32}(undef, 0, 2)
    else
        _deduplicate_undirected_edges(reduce(vcat, collect(values(edge_groups))))
    end

    meta["boundary_mode"] = "$(kind)_gasket"
    meta["gasket_removed_exterior_triangle_step"] = Int(map_data.order_face[prod_idx+1])
    meta["gasket_original_vertex_ids"] = Int.(original_vertices)
    meta["gasket_num_vertices"] = length(original_vertices)
    meta["gasket_num_edges"] = size(all_edges, 1)
    meta["gasket_boundary_vertex_color"] = boundary_color
    meta["gasket_boundary_edge_colors"] = [boundary_color, fresh_color]

    return length(original_vertices), all_edges, edge_groups, boundary_vertices, triangles, meta
end

function _prepare_hc_problem(map_data::FKMap; dimension::Int, boundary_scale::Float64, hc_boundary_mode=nothing)
    edges = layout_edges(map_data; drop_loops=true)
    edge_groups = _edge_groups_for_map(map_data)
    map_variant = _hc_variant(map_data)
    metadata = Dict{String,Any}("model" => "fk", "dimension" => dimension, "variant" => map_variant)

    if dimension == 3
        return LayoutProblem(
            num_vertices=num_vertices(map_data),
            edges=edges,
            edge_groups=edge_groups,
            surface_triangles=sanitize_triangles(map_data.triangulation_faces; drop_degenerate=true, deduplicate=true),
            metadata=metadata,
            render_map_data=map_data,
        )
    end

    boundary_mode = lowercase(strip(string(something(hc_boundary_mode, map_variant == "spanning_tree" ? "face" : "h_gasket"))))

    if boundary_mode in ("h", "c", "h_gasket", "c_gasket")
        nverts, gasket_edges, gasket_edge_groups, boundary_vertices, triangles, extra = _hc_gasket_subgraph(map_data, boundary_mode)
        merge!(metadata, extra)
        boundary_positions = length(boundary_vertices) >= 3 ? nothing : _regular_polygon(length(boundary_vertices), max(boundary_scale, 1.0))
        if boundary_positions !== nothing
            metadata["gasket_boundary_fallback"] = "regular_polygon"
        end
        return LayoutProblem(
            num_vertices=nverts,
            edges=gasket_edges,
            edge_groups=gasket_edge_groups,
            boundary_vertices=boundary_vertices,
            boundary_positions=boundary_positions,
            surface_triangles=triangles,
            metadata=metadata,
            render_map_data=nothing,
        )
    end

    tri_faces = sanitize_triangles(map_data.triangulation_faces; drop_degenerate=true, deduplicate=true)
    choice = _first_simple_face(tri_faces, 3)
    choice === nothing && throw(ArgumentError("failed to find a simple FK triangular face for 2D boundary conditions"))
    idx, boundary = choice
    length(boundary) >= 3 || throw(ArgumentError("failed to find a simple FK triangular face for 2D boundary conditions"))
    boundary_vertices = Int32.(boundary[1:3])
    boundary_positions = equilateral_triangle(boundary_scale)

    if size(tri_faces, 1) > 0
        keep = trues(size(tri_faces, 1))
        keep[idx+1] = false
        tri_faces = tri_faces[keep, :]
    end

    metadata["boundary_mode"] = "face"
    metadata["outer_face_removed"] = true
    metadata["outer_triangle_index"] = idx

    return LayoutProblem(
        num_vertices=num_vertices(map_data),
        edges=edges,
        edge_groups=edge_groups,
        boundary_vertices=boundary_vertices,
        boundary_positions=boundary_positions,
        surface_triangles=tri_faces,
        metadata=metadata,
        render_map_data=map_data,
    )
end

function prepare_layout_problem(map_data; dimension::Integer, boundary_scale=nothing, options=Dict{String,Any}())
    dim = Int(dimension)
    dim in (2, 3) || throw(ArgumentError("dimension must be 2 or 3"))
    scale = boundary_scale === nothing ? max(sqrt(max(num_vertices(map_data), 1)), 2.0) : float(boundary_scale)
    opts = Dict{String,Any}(string(k) => v for (k, v) in pairs(options))

    if map_data isa SchnyderMap
        return _prepare_schnyder_problem(map_data; dimension=dim, boundary_scale=scale)
    elseif map_data isa FKMap
        return _prepare_hc_problem(map_data; dimension=dim, boundary_scale=scale, hc_boundary_mode=get(opts, "hc_boundary_mode", nothing))
    elseif map_data isa UniformMap
        return _prepare_uniform_problem(map_data; dimension=dim, boundary_scale=scale, outer_face_index=get(opts, "outer_face_index", nothing))
    else
        edges = layout_edges(map_data; drop_loops=true)
        return LayoutProblem(
            num_vertices=num_vertices(map_data),
            edges=edges,
            edge_groups=Dict("generic" => edges),
            metadata=Dict("model" => "generic", "dimension" => dim),
            render_map_data=map_data,
        )
    end
end

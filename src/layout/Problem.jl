# Model-specific preparation of graph layout problems.

Base.@kwdef struct LayoutProblem
    num_vertices::Int
    edges::Matrix{Int32}
    edge_groups::Dict{String,Matrix{Int32}}
    boundary_vertices::Union{Nothing,Vector{Int32}} = nothing
    boundary_positions::Union{Nothing,Matrix{Float64}} = nothing
    faces::Any = nothing
    surface_triangles::Union{Nothing,Matrix{Int32}} = nothing
    surface_triangle_edge_ids::Union{Nothing,Matrix{Int32}} = nothing
    packing_boundary_vertices::Union{Nothing,Vector{Int32}} = nothing
    packing_faces::Any = nothing
    packing_triangles::Union{Nothing,Matrix{Int32}} = nothing
    packing_triangle_edge_ids::Union{Nothing,Matrix{Int32}} = nothing
    metadata::Dict{String,Any} = Dict{String,Any}()
    render_map_data::Any = nothing
    render_vertex_indices::Union{Nothing,Vector{Int32}} = nothing
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

function _reindex_edges(edge_array, old_to_new; deduplicate::Bool=true)
    arr = sanitize_edge_array(edge_array)
    if size(arr, 1) == 0
        return Matrix{Int32}(undef, 0, 2)
    end
    remapped = Matrix{Int32}(undef, size(arr, 1), 2)
    for i in 1:size(arr, 1)
        remapped[i, 1] = Int32(old_to_new[Int(arr[i, 1])])
        remapped[i, 2] = Int32(old_to_new[Int(arr[i, 2])])
    end
    return deduplicate ? _deduplicate_undirected_edges(remapped) : remapped
end

function _reindex_triangles(triangles, old_to_new; deduplicate::Bool=true)
    tri = Int32.(triangles)
    size(tri, 1) == 0 && return Matrix{Int32}(undef, 0, 3)
    remapped = Matrix{Int32}(undef, size(tri, 1), 3)
    for i in 1:size(tri, 1)
        for j in 1:3
            remapped[i, j] = Int32(old_to_new[Int(tri[i, j])])
        end
    end
    return sanitize_triangles(remapped; drop_degenerate=true, deduplicate=deduplicate)
end

function _triangle_side_edge_ids(triangles, green_edges, diag_edges, edge_u, edge_v)
    tri = Int32.(triangles)
    greens = Int32.(green_edges)
    diags = Int32.(diag_edges)
    size(tri, 1) == 0 && return Matrix{Int32}(undef, 0, 3)

    side_edges = Matrix{Int32}(undef, size(tri, 1), 3)
    for i in 1:size(tri, 1)
        a = tri[i, 1]
        b = tri[i, 2]
        c = tri[i, 3]
        side_pairs = (
            (min(a, b), max(a, b)),
            (min(b, c), max(b, c)),
            (min(c, a), max(c, a)),
        )
        triangle_edges = Int32[greens[i, 1], greens[i, 2], diags[i]]
        assigned = falses(3)

        for edge_id in triangle_edges
            edge_id >= 0 || throw(ArgumentError("triangle edge ids must be nonnegative"))
            pair = (
                min(edge_u[Int(edge_id)+1], edge_v[Int(edge_id)+1]),
                max(edge_u[Int(edge_id)+1], edge_v[Int(edge_id)+1]),
            )

            matched = 0
            for j in 1:3
                if !assigned[j] && side_pairs[j] == pair
                    matched = j
                    break
                end
            end
            matched != 0 || throw(ArgumentError("failed to match an FK triangle edge id to a triangle side"))
            side_edges[i, matched] = edge_id
            assigned[matched] = true
        end

        all(assigned) || throw(ArgumentError("incomplete FK triangle edge-id assignment"))
    end

    return side_edges
end

function _fan_triangulate_faces_with_edge_ids(faces, edge_u, edge_v; drop_degenerate::Bool=true)
    actual_pair_to_id = Dict{Tuple{Int32,Int32},Int32}()
    for i in eachindex(edge_u)
        u = Int32(edge_u[i])
        v = Int32(edge_v[i])
        u == v && continue
        pair = (min(u, v), max(u, v))
        haskey(actual_pair_to_id, pair) || (actual_pair_to_id[pair] = Int32(i - 1))
    end

    triangles = NTuple{3,Int32}[]
    triangle_edge_ids = NTuple{3,Int32}[]
    next_synthetic = Int32(length(edge_u))

    for face in _iter_faces(faces)
        length(face) < 3 && continue
        anchor = Int32(face[1])
        synthetic_pair_to_id = Dict{Tuple{Int32,Int32},Int32}()

        local function side_id(u::Int32, v::Int32)
            pair = (min(u, v), max(u, v))
            if haskey(actual_pair_to_id, pair)
                return actual_pair_to_id[pair]
            end
            if haskey(synthetic_pair_to_id, pair)
                return synthetic_pair_to_id[pair]
            end
            eid = next_synthetic
            next_synthetic += 1
            synthetic_pair_to_id[pair] = eid
            return eid
        end

        for i in 2:(length(face) - 1)
            a = anchor
            b = Int32(face[i])
            c = Int32(face[i + 1])
            tri = (a, b, c)
            if drop_degenerate && !_triangle_has_distinct_vertices(tri)
                continue
            end
            push!(triangles, tri)
            push!(triangle_edge_ids, (
                side_id(a, b),
                side_id(b, c),
                side_id(c, a),
            ))
        end
    end

    if isempty(triangles)
        return Matrix{Int32}(undef, 0, 3), Matrix{Int32}(undef, 0, 3)
    end

    tri = Matrix{Int32}(undef, length(triangles), 3)
    tri_ids = Matrix{Int32}(undef, length(triangle_edge_ids), 3)
    for i in eachindex(triangles)
        tri[i, 1] = triangles[i][1]
        tri[i, 2] = triangles[i][2]
        tri[i, 3] = triangles[i][3]
        tri_ids[i, 1] = triangle_edge_ids[i][1]
        tri_ids[i, 2] = triangle_edge_ids[i][2]
        tri_ids[i, 3] = triangle_edge_ids[i][3]
    end
    return sanitize_triangles(tri; drop_degenerate=drop_degenerate, deduplicate=false), tri_ids
end

function _fan_triangulate_faces_with_face_edge_ids(faces, face_edge_ids; next_synthetic::Integer=0, drop_degenerate::Bool=true)
    triangles = NTuple{3,Int32}[]
    triangle_edge_ids = NTuple{3,Int32}[]
    next_diag = Int32(next_synthetic)

    face_iter = collect(_iter_faces(faces))
    length(face_iter) == length(face_edge_ids) || throw(ArgumentError("face_edge_ids must align with faces"))

    for (face, edge_cycle_raw) in zip(face_iter, face_edge_ids)
        edge_cycle = Int32.(edge_cycle_raw)
        length(face) == length(edge_cycle) || throw(ArgumentError("each face_edge_ids cycle must align with the face vertex cycle"))
        length(face) < 3 && continue
        anchor = Int32(face[1])
        diag_ids = Dict{Int,Int32}()

        for i in 2:(length(face) - 1)
            a = anchor
            b = Int32(face[i])
            c = Int32(face[i + 1])
            tri = (a, b, c)
            if drop_degenerate && !_triangle_has_distinct_vertices(tri)
                continue
            end

            ab = if i == 2
                edge_cycle[1]
            else
                get!(diag_ids, i, begin
                    eid = next_diag
                    next_diag += 1
                    eid
                end)
            end
            bc = edge_cycle[i]
            ca = if i == length(face) - 1
                edge_cycle[end]
            else
                get!(diag_ids, i + 1, begin
                    eid = next_diag
                    next_diag += 1
                    eid
                end)
            end

            push!(triangles, tri)
            push!(triangle_edge_ids, (ab, bc, ca))
        end
    end

    if isempty(triangles)
        return Matrix{Int32}(undef, 0, 3), Matrix{Int32}(undef, 0, 3)
    end

    tri = Matrix{Int32}(undef, length(triangles), 3)
    tri_ids = Matrix{Int32}(undef, length(triangle_edge_ids), 3)
    for i in eachindex(triangles)
        tri[i, 1] = triangles[i][1]
        tri[i, 2] = triangles[i][2]
        tri[i, 3] = triangles[i][3]
        tri_ids[i, 1] = triangle_edge_ids[i][1]
        tri_ids[i, 2] = triangle_edge_ids[i][2]
        tri_ids[i, 3] = triangle_edge_ids[i][3]
    end
    return sanitize_triangles(tri; drop_degenerate=drop_degenerate, deduplicate=false), tri_ids
end

function _undirected_edge_key_set(edges)
    arr = sanitize_edge_array(edges)
    keys = Set{Tuple{Int32,Int32}}()
    for i in 1:size(arr, 1)
        u = arr[i, 1]
        v = arr[i, 2]
        u == v && continue
        push!(keys, (min(u, v), max(u, v)))
    end
    return keys
end

function _map_render_edges_to_auxiliary(edge_array, render_vertex_indices)
    arr = sanitize_edge_array(edge_array)
    size(arr, 1) == 0 && return Matrix{Int32}(undef, 0, 2)

    interface = Int32.(collect(render_vertex_indices))
    mapped = Matrix{Int32}(undef, size(arr, 1), 2)
    for i in 1:size(arr, 1)
        mapped[i, 1] = interface[Int(arr[i, 1]) + 1]
        mapped[i, 2] = interface[Int(arr[i, 2]) + 1]
    end
    return _deduplicate_undirected_edges(mapped)
end

function _augment_closed_meandric_sfdp_graph(map_data::UniformMeandricSystemMap, triangles)
    base_edges = auxiliary_layout_edges(map_data; drop_loops=true)
    render_support_edges = _map_render_edges_to_auxiliary(layout_edges(map_data; drop_loops=true), map_data.render_vertex_indices)
    triangulation_edges = _map_render_edges_to_auxiliary(_triangle_edge_array(triangles), map_data.render_vertex_indices)

    final_edges = first(collapse_undirected_edges(vcat(base_edges, render_support_edges, triangulation_edges); drop_loops=true))

    base_keys = _undirected_edge_key_set(base_edges)
    support_keys = _undirected_edge_key_set(render_support_edges)
    triangulation_keys = _undirected_edge_key_set(triangulation_edges)
    final_keys = _undirected_edge_key_set(final_edges)

    triangulation_only = Set{Tuple{Int32,Int32}}()
    for edge in triangulation_keys
        (edge in base_keys || edge in support_keys) && continue
        push!(triangulation_only, edge)
    end

    metadata = Dict{String,Any}(
        "layout_graph" => "glued_trees_augmented",
        "layout_graph_base" => "glued_trees",
        "layout_graph_support" => "render_edges",
        "layout_graph_augmentation" => "triangulated_faces",
        "layout_base_edge_count" => length(base_keys),
        "layout_render_support_edge_count" => length(support_keys),
        "layout_triangulation_edge_count" => length(triangulation_keys),
        "layout_temporary_triangulation_edge_count" => length(triangulation_only),
        "layout_augmented_edge_count" => length(final_keys),
    )
    return final_edges, metadata
end

function _pinched_quad_collapsed_edge(face)
    verts = Int32.(collect(face))
    length(verts) == 4 || return nothing
    length(unique(verts)) == 3 || return nothing

    for shift in 0:3
        idx = ntuple(i -> mod1(i + shift, 4), 4)
        rot = verts[[idx...]]
        if rot[1] == rot[3] && rot[1] != rot[2] && rot[1] != rot[4] && rot[2] != rot[4]
            return Int32[rot[1], rot[4]]
        end
    end

    return nothing
end

function _filter_collapsible_quad_faces_with_edge_ids(faces, face_edge_ids)
    face_iter = [Int32.(collect(face)) for face in _iter_faces(faces)]
    length(face_iter) == length(face_edge_ids) || throw(ArgumentError("face_edge_ids must align with faces"))

    kept_faces = Vector{Vector{Int32}}()
    kept_face_edge_ids = Vector{Vector{Int32}}()
    keep_mask = falses(length(face_iter))
    collapsed_edges = Matrix{Int32}(undef, 0, 2)

    collapsed_pairs = Tuple{Int32,Int32}[]
    for i in eachindex(face_iter)
        face = face_iter[i]
        edge_cycle = Int32.(collect(face_edge_ids[i]))
        length(face) == length(edge_cycle) || throw(ArgumentError("each face_edge_ids cycle must align with the face vertex cycle"))

        collapsed = _pinched_quad_collapsed_edge(face)
        if collapsed === nothing
            keep_mask[i] = true
            push!(kept_faces, face)
            push!(kept_face_edge_ids, edge_cycle)
        else
            push!(collapsed_pairs, (collapsed[1], collapsed[2]))
        end
    end

    if !isempty(collapsed_pairs)
        collapsed_edges = Matrix{Int32}(undef, length(collapsed_pairs), 2)
        for (i, (u, v)) in enumerate(collapsed_pairs)
            collapsed_edges[i, 1] = u
            collapsed_edges[i, 2] = v
        end
        collapsed_edges = _deduplicate_undirected_edges(collapsed_edges)
    end

    return kept_faces, kept_face_edge_ids, keep_mask, collapsed_edges
end

function _surface_triangle_data(map_data::FKMap)
    tri_faces_all = Int32.(map_data.triangulation_faces)
    tri_green_all = Int32.(map_data.triangle_green_edges)
    tri_diag_all = Int32.(map_data.triangle_diag_edge)

    nrows = size(tri_faces_all, 1)
    size(tri_green_all, 1) == nrows || throw(ArgumentError("triangle_green_edges must align with triangulation_faces"))
    length(tri_diag_all) == nrows || throw(ArgumentError("triangle_diag_edge must align with triangulation_faces"))

    valid = BitVector(undef, nrows)
    for i in 1:nrows
        valid[i] = _triangle_has_distinct_vertices(view(tri_faces_all, i, :))
    end

    tri_faces = tri_faces_all[valid, :]
    tri_green = tri_green_all[valid, :]
    tri_diag = tri_diag_all[valid]
    triangle_edge_ids = _triangle_side_edge_ids(tri_faces, tri_green, tri_diag, map_data.edge_u, map_data.edge_v)
    return tri_faces, triangle_edge_ids
end

@inline function _triangle_has_distinct_vertices(tri)
    return (tri[1] != tri[2]) && (tri[2] != tri[3]) && (tri[1] != tri[3])
end

function _drop_self_loops(edge_array)
    arr = sanitize_edge_array(edge_array)
    size(arr, 1) == 0 && return arr
    keep = arr[:, 1] .!= arr[:, 2]
    return arr[keep, :]
end

function _gasket_cleanup_mask(
    map_data::FKMap,
    selected_triangles_old,
    selected_green,
    selected_diag,
)
    nrows = size(selected_triangles_old, 1)
    valid = BitVector(undef, nrows)
    supported = Dict{Int32,Int}()

    for i in 1:nrows
        tri = view(selected_triangles_old, i, :)
        valid[i] = _triangle_has_distinct_vertices(tri)
        valid[i] || continue
        for e in selected_green[i, :]
            e < 0 && continue
            supported[e] = get(supported, e, 0) + 1
        end
        d = selected_diag[i]
        d < 0 && continue
        supported[d] = get(supported, d, 0) + 1
    end

    drop_edge_ids = Set{Int32}()
    for i in 1:nrows
        valid[i] && continue

        g1 = selected_green[i, 1]
        g2 = selected_green[i, 2]
        if g1 >= 0 || g2 >= 0
            candidate = Int32(-1)
            if g1 >= 0 && g2 >= 0 && g1 != g2
                s1 = get(supported, g1, 0)
                s2 = get(supported, g2, 0)
                candidate = if s1 < s2
                    g1
                elseif s2 < s1
                    g2
                else
                    max(g1, g2)
                end
            else
                candidate = max(g1, g2)
            end
            candidate >= 0 && push!(drop_edge_ids, candidate)
        end

        d = selected_diag[i]
        if d >= 0 && map_data.edge_u[Int(d)+1] == map_data.edge_v[Int(d)+1]
            push!(drop_edge_ids, d)
        end
    end

    return valid, drop_edge_ids
end

abstract type _CurrentGasketItem end

struct _CurrentGasketPlain <: _CurrentGasketItem
    ch::Char
    step::Int
end

mutable struct _CurrentGasketHole <: _CurrentGasketItem
    kind::Char
    start_step::Int
    end_step::Int
    children::Vector{_CurrentGasketItem}
end

mutable struct _GasketLayoutVertex
    is_hamburger::Bool
    darts::Vector{Int32}
end

mutable struct _GasketLayoutState
    vertices::Vector{_GasketLayoutVertex}
    dart_tail::Vector{Int32}
    dart_head::Vector{Int32}
    dart_twin::Vector{Int32}
    dart_edge::Vector{Int32}
    edge_u::Vector{Int32}
    edge_v::Vector{Int32}
    edge_kind::Vector{String}
    current_v_h::Int32
    current_v_c::Int32
    prev_v_h::Vector{Int32}
    prev_v_c::Vector{Int32}
    prev_v_f::Vector{Int32}
    boundary_cycle::Vector{Int32}
end

_current_gasket_start_char(kind::Char) = kind == 'h' ? 'c' : 'h'
_current_gasket_order_char(kind::Char) = kind == 'h' ? 'H' : 'C'
_current_gasket_child_kind(kind::Char) = kind == 'h' ? 'c' : 'h'

function _current_gasket_plain_allowed(kind::Char)
    return kind == 'h' ? Set(('c', 'C')) : Set(('h', 'H'))
end

function _gasket_pairing_data(map_data::FKMap)
    phi = Dict{Int,Int}()
    h_h_steps = Set{Int}()
    c_c_steps = Set{Int}()
    f_h_steps = Set{Int}()
    f_c_steps = Set{Int}()

    for i in eachindex(map_data.production_steps)
        prod_step = Int(map_data.production_steps[i])
        order_step = Int(map_data.order_steps[i])
        phi[order_step] = prod_step
        if map_data.order_symbols[i] == 'H'
            push!(h_h_steps, order_step)
        elseif map_data.order_symbols[i] == 'C'
            push!(c_c_steps, order_step)
        elseif map_data.order_symbols[i] == 'F'
            if map_data.fulfilled_by[i] == 'h'
                push!(f_h_steps, order_step)
            elseif map_data.fulfilled_by[i] == 'c'
                push!(f_c_steps, order_step)
            end
        end
    end

    return Dict{String,Any}(
        "phi" => phi,
        "h_h_steps" => h_h_steps,
        "c_c_steps" => c_c_steps,
        "f_h_steps" => f_h_steps,
        "f_c_steps" => f_c_steps,
    )
end

function _current_gasket_step_indices(start_step::Int, end_step::Int, kind::Char, pair_data)
    phi = pair_data["phi"]
    h_h_steps = pair_data["h_h_steps"]
    c_c_steps = pair_data["c_c_steps"]
    f_h_steps = pair_data["f_h_steps"]
    f_c_steps = pair_data["f_c_steps"]

    gasket_idx = collect(start_step+1:end_step-1)
    if kind == 'h'
        for f in sort!(collect(intersect(f_c_steps, Set(gasket_idx))))
            setdiff!(gasket_idx, phi[f]:f)
        end
        for h_h in sort!(collect(intersect(h_h_steps, Set(gasket_idx))))
            setdiff!(gasket_idx, [phi[h_h], h_h])
        end
    else
        for f in sort!(collect(intersect(f_h_steps, Set(gasket_idx))))
            setdiff!(gasket_idx, phi[f]:f)
        end
        for c_c in sort!(collect(intersect(c_c_steps, Set(gasket_idx))))
            setdiff!(gasket_idx, [phi[c_c], c_c])
        end
    end
    reverse!(gasket_idx)
    return gasket_idx
end

function _build_current_gasket_tree_children(map_data::FKMap, pair_data, kind::Char, reversed_steps::Vector{Int})
    phi = pair_data["phi"]
    allowed = _current_gasket_plain_allowed(kind)
    child_kind = _current_gasket_child_kind(kind)
    work = copy(reversed_steps)
    children = _CurrentGasketItem[]

    for step in copy(reversed_steps)
        step in work || continue
        ch = map_data.word[step+1]
        if ch == 'F'
            start = phi[step]
            inner = Int[q for q in work if start < q < step]
            child_children = _build_current_gasket_tree_children(map_data, pair_data, child_kind, inner)
            push!(children, _CurrentGasketHole(child_kind, start, step, child_children))
            setdiff!(work, start:step-1)
        elseif ch in allowed
            push!(children, _CurrentGasketPlain(ch, step))
        end
    end

    reverse!(children)
    return children
end

function _build_current_gasket_tree(map_data::FKMap, kind::Char, start_step::Int, end_step::Int)
    pair_data = _gasket_pairing_data(map_data)
    reversed_steps = _current_gasket_step_indices(start_step, end_step, kind, pair_data)
    children = _build_current_gasket_tree_children(map_data, pair_data, kind, reversed_steps)
    return _CurrentGasketHole(kind, start_step, end_step, children)
end

_current_gasket_original_size(hole::_CurrentGasketHole) = length(hole.children) + 1

function _current_gasket_plain_perimeter(ch::Char)
    if ch == 'H' || ch == 'C'
        return 1
    elseif ch == 'h' || ch == 'c'
        return -1
    end
    throw(ArgumentError("unexpected gasket symbol $(repr(ch))"))
end

function _current_gasket_hole_perimeter(hole::_CurrentGasketHole)
    result = 0
    for child in hole.children
        if child isa _CurrentGasketPlain
            result += _current_gasket_plain_perimeter((child::_CurrentGasketPlain).ch)
        else
            result += _current_gasket_hole_perimeter(child::_CurrentGasketHole)
        end
    end
    return result
end

function _current_gasket_word(item::_CurrentGasketItem)
    if item isa _CurrentGasketPlain
        return string(item.ch)
    end
    hole = item::_CurrentGasketHole
    return string(_current_gasket_start_char(hole.kind)) * join(_current_gasket_word(child) for child in hole.children) * "F"
end

function _count_small_current_gasket_holes(item::_CurrentGasketItem)
    item isa _CurrentGasketPlain && return (0, 0)
    hole = item::_CurrentGasketHole
    size_one = 0
    size_two = 0
    for child in hole.children
        child isa _CurrentGasketHole || continue
        perimeter = _current_gasket_hole_perimeter(child)
        if perimeter == 0
            size_one += 1
        elseif perimeter == 1
            size_two += 1
        end
        child_one, child_two = _count_small_current_gasket_holes(child)
        size_one += child_one
        size_two += child_two
    end
    return size_one, size_two
end

function _collapse_small_current_gasket_children(item::_CurrentGasketItem)
    if item isa _CurrentGasketPlain
        return _CurrentGasketItem[item]
    end

    hole = item::_CurrentGasketHole
    perimeter = _current_gasket_hole_perimeter(hole)
    if perimeter == 0
        return _CurrentGasketItem[]
    elseif perimeter == 1
        return _CurrentGasketItem[_CurrentGasketPlain(_current_gasket_order_char(hole.kind), hole.start_step)]
    end

    cleaned_children = _CurrentGasketItem[]
    for child in hole.children
        append!(cleaned_children, _collapse_small_current_gasket_children(child))
    end
    return _CurrentGasketItem[_CurrentGasketHole(hole.kind, hole.start_step, hole.end_step, cleaned_children)]
end

function _clean_current_gasket_root(root::_CurrentGasketHole)
    cleaned_children = _CurrentGasketItem[]
    for child in root.children
        append!(cleaned_children, _collapse_small_current_gasket_children(child))
    end
    return _CurrentGasketHole(root.kind, root.start_step, root.end_step, cleaned_children)
end

function _current_gasket_boundary_order_steps(map_data::FKMap, hole::_CurrentGasketHole)
    boundary_order = hole.kind == 'h' ? 'H' : 'C'
    boundary_steps = Dict{Int,Char}()
    order = sortperm(map_data.order_steps)
    for j in order
        order_step = Int(map_data.order_steps[j])
        if order_step <= hole.start_step || order_step > hole.end_step
            continue
        end
        if map_data.order_symbols[j] != boundary_order
            continue
        end
        if Int(map_data.production_steps[j]) >= hole.start_step
            continue
        end
        boundary_steps[order_step] = boundary_order
    end
    return boundary_steps
end

function _fresh_order_pair_map(map_data::FKMap)
    pair_map = Dict{Int,Int}()
    for i in eachindex(map_data.production_steps)
        map_data.order_symbols[i] == 'F' || continue
        pair_map[Int(map_data.production_steps[i])] = Int(map_data.order_steps[i])
    end
    return pair_map
end

function _clean_current_gasket_fk_span(
    map_data::FKMap,
    start_step::Int,
    end_step::Int,
    pair_map::Dict{Int,Int};
    collapse_self::Bool=false,
)
    haskey(pair_map, start_step) || throw(ArgumentError("expected an FK fresh-order span starting at step $start_step"))
    pair_map[start_step] == end_step || throw(ArgumentError("FK fresh-order span endpoints do not match the expected gasket span"))

    pieces = String[]
    pos = start_step
    while pos <= end_step
        if pos == start_step || pos == end_step
            push!(pieces, string(map_data.word[pos+1]))
            pos += 1
            continue
        end

        child_end = get(pair_map, pos, -1)
        if child_end != -1 && child_end <= end_step
            push!(pieces, _clean_current_gasket_fk_span(map_data, pos, child_end, pair_map; collapse_self=true))
            pos = child_end + 1
        else
            push!(pieces, string(map_data.word[pos+1]))
            pos += 1
        end
    end

    cleaned = join(pieces)
    if collapse_self
        reduced = _reduce_fk_fragment(cleaned)
        if length(reduced) <= 1
            return reduced
        end
    end
    return cleaned
end

function _clean_current_gasket_fk_subword(map_data::FKMap, hole::_CurrentGasketHole; collapse_self::Bool=false)
    pair_map = _fresh_order_pair_map(map_data)
    cleaned = _clean_current_gasket_fk_span(
        map_data,
        hole.start_step,
        hole.end_step,
        pair_map;
        collapse_self=collapse_self,
    )
    if collapse_self
        reduced = _reduce_fk_fragment(cleaned)
        if length(reduced) <= 1
            return _reduce_fk_fragment(cleaned)
        end
    end
    return cleaned
end

function _decorate_current_gasket_tree(map_data::FKMap, hole::_CurrentGasketHole)
    hole_children = Dict{Int,_CurrentGasketHole}()
    plain_steps = Dict{Int,_CurrentGasketPlain}()
    for child in hole.children
        if child isa _CurrentGasketHole
            decorated_child = _decorate_current_gasket_tree(map_data, child)
            hole_children[decorated_child.start_step] = decorated_child
        else
            plain = child::_CurrentGasketPlain
            plain_steps[plain.step] = plain
        end
    end
    boundary_steps = _current_gasket_boundary_order_steps(map_data, hole)

    children = _CurrentGasketItem[]
    pos = hole.start_step + 1
    while pos < hole.end_step
        child = get(hole_children, pos, nothing)
        if child !== nothing
            push!(children, child)
            pos = child.end_step + 1
            continue
        end

        plain = get(plain_steps, pos, nothing)
        if plain !== nothing
            push!(children, plain)
            pos += 1
            continue
        end

        boundary_ch = get(boundary_steps, pos, '\0')
        if boundary_ch != '\0'
            push!(children, _CurrentGasketPlain(boundary_ch, pos))
            pos += 1
            continue
        end

        pos += 1
    end

    return _CurrentGasketHole(hole.kind, hole.start_step, hole.end_step, children)
end

function _reduce_fk_fragment(word::AbstractString)
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
                next_stack[top_any+1] = idx
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
            top_any == -1 ? ('\0', -1) : (prod_symbols[top_any+1], top_any)
        else
            throw(ArgumentError("unexpected FK symbol $(repr(ch))"))
        end

        if idx == -1
            push!(unmatched_orders, ch)
            continue
        end

        ps = prev_stack[idx+1]
        ns = next_stack[idx+1]
        if ps != -1
            next_stack[ps+1] = ns
        end
        if ns != -1
            prev_stack[ns+1] = ps
        else
            top_any = ps
        end
        top_by_type[target_type] = prev_same[idx+1]
    end

    unmatched_burgers = Char[]
    idx = top_any
    while idx != -1 && prev_stack[idx+1] != -1
        idx = prev_stack[idx+1]
    end
    while idx != -1
        push!(unmatched_burgers, prod_symbols[idx+1])
        idx = next_stack[idx+1]
    end

    return String(vcat(unmatched_orders, unmatched_burgers))
end

_matching_fk_production(ch::Char) = ch == 'H' ? 'h' : ch == 'C' ? 'c' : throw(ArgumentError("expected H or C in reduced FK fragment"))

function _build_cleaned_fk_gasket_map(cleaned_subword::AbstractString)
    reduced = _reduce_fk_fragment(cleaned_subword)
    prefix = isempty(reduced) ? "" : String([_matching_fk_production(ch) for ch in reduced])
    balanced_word = prefix * String(cleaned_subword)
    return build_hc_map_from_word(balanced_word), reduced, balanced_word
end

function _new_gasket_layout_state(root_kind::Char)
    root_is_hamburger = root_kind == 'c'
    return _GasketLayoutState(
        [_GasketLayoutVertex(!root_is_hamburger, Int32[])],
        Int32[],
        Int32[],
        Int32[],
        Int32[],
        Int32[],
        Int32[],
        String[],
        !root_is_hamburger ? Int32(0) : Int32(-1),
        root_is_hamburger ? Int32(0) : Int32(-1),
        Int32[],
        Int32[],
        Int32[],
        Int32[],
    )
end

function _gasket_add_vertex!(state::_GasketLayoutState, is_hamburger::Bool)
    push!(state.vertices, _GasketLayoutVertex(is_hamburger, Int32[]))
    return Int32(length(state.vertices) - 1)
end

_gasket_vertex_degree(state::_GasketLayoutState, v::Int32) = length(state.vertices[Int(v)+1].darts)

function _gasket_add_edge!(state::_GasketLayoutState, v::Int32, w::Int32; insert_at::Int=0, boundary::Bool=false)
    (v < 0 || w < 0 || v == w) && return Int32(-1)

    kind = if boundary
        "boundary"
    elseif state.vertices[Int(v)+1].is_hamburger && state.vertices[Int(w)+1].is_hamburger
        "red"
    elseif !state.vertices[Int(v)+1].is_hamburger && !state.vertices[Int(w)+1].is_hamburger
        "blue"
    else
        "green"
    end

    push!(state.edge_u, v)
    push!(state.edge_v, w)
    push!(state.edge_kind, kind)
    edge_id = Int32(length(state.edge_u) - 1)

    d_vw = Int32(length(state.dart_tail) + 1)
    d_wv = Int32(length(state.dart_tail) + 2)
    append!(state.dart_tail, [v, w])
    append!(state.dart_head, [w, v])
    append!(state.dart_twin, [d_wv, d_vw])
    append!(state.dart_edge, [edge_id, edge_id])

    push!(state.vertices[Int(v)+1].darts, d_vw)
    if insert_at <= 0
        push!(state.vertices[Int(w)+1].darts, d_wv)
    else
        splice!(state.vertices[Int(w)+1].darts, insert_at:insert_at-1, d_wv)
    end
    return edge_id
end

function _process_current_gasket_plain!(state::_GasketLayoutState, ch::Char)
    if ch == 'h'
        pred = state.current_v_h
        opp = state.current_v_c
        push!(state.prev_v_h, pred)
        nxt = _gasket_add_vertex!(state, true)
        _gasket_add_edge!(state, pred, nxt)
        _gasket_add_edge!(state, opp, nxt)
        state.current_v_h = nxt
    elseif ch == 'c'
        pred = state.current_v_c
        opp = state.current_v_h
        push!(state.prev_v_c, pred)
        nxt = _gasket_add_vertex!(state, false)
        _gasket_add_edge!(state, pred, nxt)
        _gasket_add_edge!(state, opp, nxt)
        state.current_v_c = nxt
    elseif ch == 'H'
        if isempty(state.prev_v_h)
            push!(state.boundary_cycle, state.current_v_h)
            nxt = _gasket_add_vertex!(state, true)
            _gasket_add_edge!(state, state.current_v_h, nxt; boundary=true)
            _gasket_add_edge!(state, state.current_v_c, nxt)
            state.current_v_h = nxt
        else
            state.current_v_h = pop!(state.prev_v_h)
            _gasket_add_edge!(state, state.current_v_c, state.current_v_h)
        end
    elseif ch == 'C'
        if isempty(state.prev_v_c)
            push!(state.boundary_cycle, state.current_v_c)
            nxt = _gasket_add_vertex!(state, false)
            _gasket_add_edge!(state, state.current_v_c, nxt; boundary=true)
            _gasket_add_edge!(state, state.current_v_h, nxt)
            state.current_v_c = nxt
        else
            state.current_v_c = pop!(state.prev_v_c)
            _gasket_add_edge!(state, state.current_v_h, state.current_v_c)
        end
    else
        throw(ArgumentError("unexpected gasket symbol $(repr(ch))"))
    end
end

function _enter_current_gasket_hole!(state::_GasketLayoutState, kind::Char; is_root::Bool)
    if kind == 'c'
        if is_root
            nxt = _gasket_add_vertex!(state, true)
            _gasket_add_edge!(state, state.current_v_c, nxt)
            state.current_v_h = nxt
            return 1
        end
        push!(state.prev_v_f, state.current_v_c)
        push!(state.prev_v_h, state.current_v_h)
        nxt = _gasket_add_vertex!(state, true)
        _gasket_add_edge!(state, state.current_v_c, nxt)
        prev_deg = _gasket_vertex_degree(state, state.current_v_c)
        state.current_v_h = nxt
        return prev_deg
    end

    if is_root
        nxt = _gasket_add_vertex!(state, false)
        _gasket_add_edge!(state, state.current_v_h, nxt)
        state.current_v_c = nxt
        return 1
    end
    push!(state.prev_v_f, state.current_v_h)
    push!(state.prev_v_c, state.current_v_c)
    nxt = _gasket_add_vertex!(state, false)
    _gasket_add_edge!(state, state.current_v_h, nxt)
    prev_deg = _gasket_vertex_degree(state, state.current_v_h)
    state.current_v_c = nxt
    return prev_deg
end

function _leave_current_gasket_hole!(state::_GasketLayoutState, kind::Char, prev_deg::Int; is_root::Bool)
    if kind == 'c'
        if is_root
            push!(state.boundary_cycle, state.current_v_c)
            _gasket_add_edge!(state, state.current_v_c, Int32(0); boundary=true)
            return
        end
        v_f = pop!(state.prev_v_f)
        state.current_v_h = pop!(state.prev_v_h)
        _gasket_add_edge!(state, state.current_v_c, v_f; insert_at=prev_deg)
        _gasket_add_edge!(state, state.current_v_c, state.current_v_h)
        return
    end

    if is_root
        push!(state.boundary_cycle, state.current_v_h)
        _gasket_add_edge!(state, state.current_v_h, Int32(0); boundary=true)
        return
    end
    v_f = pop!(state.prev_v_f)
    state.current_v_c = pop!(state.prev_v_c)
    _gasket_add_edge!(state, state.current_v_h, v_f; insert_at=prev_deg)
    _gasket_add_edge!(state, state.current_v_h, state.current_v_c)
end

function _walk_current_gasket_tree!(state::_GasketLayoutState, hole::_CurrentGasketHole; is_root::Bool=false)
    prev_deg = _enter_current_gasket_hole!(state, hole.kind; is_root=is_root)
    for child in hole.children
        if child isa _CurrentGasketPlain
            _process_current_gasket_plain!(state, child.ch)
        else
            _walk_current_gasket_tree!(state, child; is_root=false)
        end
    end
    _leave_current_gasket_hole!(state, hole.kind, prev_deg; is_root=is_root)
    return state
end

function _gasket_prev_darts(state::_GasketLayoutState)
    prev_dart = fill(Int32(0), length(state.dart_tail))
    for vertex in state.vertices
        darts = vertex.darts
        isempty(darts) && continue
        for i in eachindex(darts)
            prev_dart[Int(darts[i])] = darts[mod1(i - 1, length(darts))]
        end
    end
    return prev_dart
end

function _extract_gasket_faces(state::_GasketLayoutState)
    prev_dart = _gasket_prev_darts(state)
    seen = falses(length(state.dart_tail))
    faces = Vector{Vector{Int32}}()
    face_edges = Vector{Vector{Int32}}()

    for dart in eachindex(state.dart_tail)
        seen[dart] && continue
        face = Int32[]
        edges = Int32[]
        cur = dart
        while !seen[cur]
            seen[cur] = true
            push!(face, state.dart_tail[cur])
            push!(edges, state.dart_edge[cur])
            cur = Int(prev_dart[Int(state.dart_twin[cur])])
        end
        push!(faces, face)
        push!(face_edges, edges)
    end

    return faces, face_edges
end

function _select_outer_gasket_face(state::_GasketLayoutState, face_edges)
    best_idx = 1
    best_score = (-1, -1)
    for (i, edges) in enumerate(face_edges)
        boundary_count = count(e -> state.edge_kind[Int(e)+1] == "boundary", edges)
        score = (boundary_count, length(edges))
        if score > best_score
            best_idx = i
            best_score = score
        end
    end
    return best_idx
end

function _face_vertices_in_order(face::Vector{Int32})
    out = Int32[]
    seen = Set{Int32}()
    for v in face
        if v in seen
            continue
        end
        push!(seen, v)
        push!(out, v)
    end
    return out
end

function _cleaned_current_gasket_layout(kind::Char, root::_CurrentGasketHole, boundary_color::String, fresh_color::String)
    state = _new_gasket_layout_state(kind)
    _walk_current_gasket_tree!(state, root; is_root=true)
    faces, face_edges = _extract_gasket_faces(state)
    outer_idx = _select_outer_gasket_face(state, face_edges)
    boundary_vertices = isempty(state.boundary_cycle) ? _face_vertices_in_order(faces[outer_idx]) : copy(state.boundary_cycle)
    interior_faces = [faces[i] for i in eachindex(faces) if i != outer_idx]
    triangles = fan_triangulate_faces(interior_faces; drop_degenerate=true, deduplicate=false)

    edge_groups = Dict{String,Matrix{Int32}}(
        "green" => Matrix{Int32}(undef, 0, 2),
        "red" => Matrix{Int32}(undef, 0, 2),
        "blue" => Matrix{Int32}(undef, 0, 2),
        "orange" => Matrix{Int32}(undef, 0, 2),
        "purple" => Matrix{Int32}(undef, 0, 2),
    )

    cycle_pairs = Tuple{Int32,Int32}[]
    for i in eachindex(boundary_vertices)
        u = boundary_vertices[i]
        v = boundary_vertices[mod1(i + 1, length(boundary_vertices))]
        push!(cycle_pairs, (min(u, v), max(u, v)))
    end
    closing_pair = isempty(cycle_pairs) ? nothing : cycle_pairs[end]

    grouped = Dict(
        "green" => Tuple{Int32,Int32}[],
        "red" => Tuple{Int32,Int32}[],
        "blue" => Tuple{Int32,Int32}[],
        "orange" => Tuple{Int32,Int32}[],
        "purple" => Tuple{Int32,Int32}[],
    )

    for i in eachindex(state.edge_u)
        u = state.edge_u[i]
        v = state.edge_v[i]
        pair = (min(u, v), max(u, v))
        group = state.edge_kind[i]
        if group == "boundary"
            target = pair == closing_pair ? fresh_color : boundary_color
            push!(grouped[target], pair)
        else
            push!(grouped[group], pair)
        end
    end

    for (name, pairs) in grouped
        isempty(pairs) && continue
        arr = Matrix{Int32}(undef, length(pairs), 2)
        for (i, (u, v)) in enumerate(pairs)
            arr[i, 1] = u
            arr[i, 2] = v
        end
        edge_groups[name] = arr
    end

    all_edges = sanitize_edge_array(reduce(vcat, [arr for arr in values(edge_groups) if size(arr, 1) > 0]; init=Matrix{Int32}(undef, 0, 2)))
    return (
        length(state.vertices),
        all_edges,
        Dict(name => arr for (name, arr) in edge_groups if size(arr, 1) > 0),
        boundary_vertices,
        interior_faces,
        triangles,
        Dict{String,Any}(
            "gasket_cleanup_mode" => "tree_small_holes",
            "gasket_num_vertices" => length(state.vertices),
            "gasket_num_edges" => size(all_edges, 1),
            "gasket_num_triangles" => size(triangles, 1),
        ),
    )
end

function _fk_local_euler_path(local_adj_component::Dict{Int32,Vector{Tuple{Int32,Int}}}, start::Int32, vertex::Integer)
    wedge_ids = sort!(unique(Int[wedge_id for adjlist in values(local_adj_component) for (_, wedge_id) in adjlist]))
    used = falses(isempty(wedge_ids) ? 0 : maximum(wedge_ids))
    stack = Int32[start]
    path = Int32[]

    while !isempty(stack)
        current = stack[end]
        adjlist = get(local_adj_component, current, Tuple{Int32,Int}[])
        while !isempty(adjlist) && used[last(adjlist)[2]]
            pop!(adjlist)
        end

        if isempty(adjlist)
            push!(path, current)
            pop!(stack)
        else
            next_edge, wedge_id = pop!(adjlist)
            used[wedge_id] = true
            push!(stack, next_edge)
        end
    end

    all(used[wedge_ids]) || throw(ArgumentError("failed to recover an FK gasket cleanup path around vertex $vertex"))
    reverse!(path)
    return path
end

function _fk_gasket_bubble_vertices(num_vertices::Integer, edges, boundary_vertices, triangles, triangle_edge_ids)
    triangle_edge_ids === nothing && return Int32[]
    n = Int(num_vertices)
    tri = Int32.(triangles)
    tri_edges = Int32.(triangle_edge_ids)
    size(tri, 1) == size(tri_edges, 1) || throw(ArgumentError("triangle_edge_ids must align with FK gasket triangles"))

    adjacency = [Int32[] for _ in 1:n]
    for i in 1:size(edges, 1)
        u = Int(edges[i, 1]) + 1
        v = Int(edges[i, 2]) + 1
        u == v && continue
        push!(adjacency[u], Int32(v - 1))
        push!(adjacency[v], Int32(u - 1))
    end

    boundary_set = Set(Int32.(boundary_vertices))
    to_remove = Set{Int32}()

    for vertex in Int32.(0:(n-1))
        local_adj = Dict{Int32,Vector{Tuple{Int32,Int}}}()
        local_deg = Dict{Int32,Int}()
        local_neighbor = Dict{Int32,Int32}()
        wedge_count = 0

        function register_neighbor!(edge_id::Int32, neighbor::Int32)
            existing = get(local_neighbor, edge_id, Int32(-1))
            if existing != Int32(-1) && existing != neighbor
                throw(ArgumentError("triangle_edge_ids are inconsistent with the FK gasket triangulation around vertex $vertex"))
            end
            local_neighbor[edge_id] = neighbor
        end

        function add_local_wedge!(left_edge::Int32, left_neighbor::Int32, right_edge::Int32, right_neighbor::Int32)
            wedge_count += 1
            register_neighbor!(left_edge, left_neighbor)
            register_neighbor!(right_edge, right_neighbor)
            push!(get!(local_adj, left_edge, Tuple{Int32,Int}[]), (right_edge, wedge_count))
            push!(get!(local_adj, right_edge, Tuple{Int32,Int}[]), (left_edge, wedge_count))
            local_deg[left_edge] = get(local_deg, left_edge, 0) + 1
            local_deg[right_edge] = get(local_deg, right_edge, 0) + 1
        end

        for i in 1:size(tri, 1)
            a = tri[i, 1]
            b = tri[i, 2]
            c = tri[i, 3]
            eab = tri_edges[i, 1]
            ebc = tri_edges[i, 2]
            eca = tri_edges[i, 3]
            if a == vertex
                add_local_wedge!(eca, c, eab, b)
            end
            if b == vertex
                add_local_wedge!(eab, a, ebc, c)
            end
            if c == vertex
                add_local_wedge!(ebc, b, eca, a)
            end
        end

        visited = Set{Int32}()
        for start in keys(local_adj)
            start in visited && continue
            component = Int32[]
            queue = Int32[start]
            push!(visited, start)
            while !isempty(queue)
                current = popfirst!(queue)
                push!(component, current)
                for (next_edge, _) in get(local_adj, current, Tuple{Int32,Int}[])
                    if !(next_edge in visited)
                        push!(visited, next_edge)
                        push!(queue, next_edge)
                    end
                end
            end

            comp_set = Set(component)
            comp_adj = Dict{Int32,Vector{Tuple{Int32,Int}}}()
            odd_edges = Int32[]
            for edge_id in component
                values = [(next_edge, wedge_id) for (next_edge, wedge_id) in get(local_adj, edge_id, Tuple{Int32,Int}[]) if next_edge in comp_set]
                comp_adj[edge_id] = values
                isodd(length(values)) && push!(odd_edges, edge_id)
            end
            length(odd_edges) == 2 || continue

            shared_neighbor = local_neighbor[odd_edges[1]]
            shared_neighbor == local_neighbor[odd_edges[2]] || continue

            edge_path = _fk_local_euler_path(comp_adj, odd_edges[1], vertex)
            order = Int32[local_neighbor[edge_id] for edge_id in edge_path]
            seeds = Set(Int32[v for v in order[2:end-1] if v != shared_neighbor])
            isempty(seeds) && continue

            component_vertices = Set{Int32}()
            search = Int32[collect(seeds)...]
            valid = true
            while !isempty(search)
                current = popfirst!(search)
                current in component_vertices && continue
                if current == vertex || current == shared_neighbor
                    continue
                end
                if current in boundary_set
                    valid = false
                    break
                end
                push!(component_vertices, current)
                for next_vertex in adjacency[Int(current)+1]
                    if next_vertex != vertex && next_vertex != shared_neighbor && !(next_vertex in component_vertices)
                        push!(search, next_vertex)
                    end
                end
            end

            valid || continue
            union!(to_remove, component_vertices)
        end
    end

    return sort!(collect(to_remove))
end

function _boundary_sector_bubble_vertices(num_vertices::Integer, edges, boundary_vertices, triangles, triangle_edge_ids)
    triangle_edge_ids === nothing && return Int32[]
    n = Int(num_vertices)
    tri = Int32.(triangles)
    tri_edges = Int32.(triangle_edge_ids)
    size(tri, 1) == size(tri_edges, 1) || throw(ArgumentError("triangle_edge_ids must align with triangles"))

    boundary = Int32.(collect(boundary_vertices))
    length(boundary) >= 3 || return Int32[]
    boundary_set = Set(boundary)
    boundary_pos = Dict(v => i for (i, v) in enumerate(boundary))

    adjacency = [Int32[] for _ in 1:n]
    for i in 1:size(edges, 1)
        u = Int(edges[i, 1]) + 1
        v = Int(edges[i, 2]) + 1
        u == v && continue
        push!(adjacency[u], Int32(v - 1))
        push!(adjacency[v], Int32(u - 1))
    end

    to_remove = Set{Int32}()

    for vertex in boundary
        local_adj = Dict{Int32,Vector{Tuple{Int32,Int}}}()
        local_neighbor = Dict{Int32,Int32}()
        wedge_count = 0

        function register_neighbor!(edge_id::Int32, neighbor::Int32)
            existing = get(local_neighbor, edge_id, Int32(-1))
            if existing != Int32(-1) && existing != neighbor
                throw(ArgumentError("triangle_edge_ids are inconsistent with the triangulation around boundary vertex $vertex"))
            end
            local_neighbor[edge_id] = neighbor
        end

        function add_local_wedge!(left_edge::Int32, left_neighbor::Int32, right_edge::Int32, right_neighbor::Int32)
            wedge_count += 1
            register_neighbor!(left_edge, left_neighbor)
            register_neighbor!(right_edge, right_neighbor)
            push!(get!(local_adj, left_edge, Tuple{Int32,Int}[]), (right_edge, wedge_count))
            push!(get!(local_adj, right_edge, Tuple{Int32,Int}[]), (left_edge, wedge_count))
        end

        for i in 1:size(tri, 1)
            a = tri[i, 1]
            b = tri[i, 2]
            c = tri[i, 3]
            eab = tri_edges[i, 1]
            ebc = tri_edges[i, 2]
            eca = tri_edges[i, 3]
            if a == vertex
                add_local_wedge!(eca, c, eab, b)
            elseif b == vertex
                add_local_wedge!(eab, a, ebc, c)
            elseif c == vertex
                add_local_wedge!(ebc, b, eca, a)
            end
        end

        isempty(local_adj) && continue
        idx = boundary_pos[vertex]
        prev_boundary = boundary[idx == 1 ? length(boundary) : idx - 1]
        next_boundary = boundary[idx == length(boundary) ? 1 : idx + 1]

        components = Vector{Vector{Int32}}()
        visited = Set{Int32}()
        for start in keys(local_adj)
            start in visited && continue
            queue = Int32[start]
            push!(visited, start)
            component = Int32[]
            while !isempty(queue)
                current = popfirst!(queue)
                push!(component, current)
                for (next_edge, _) in get(local_adj, current, Tuple{Int32,Int}[])
                    if !(next_edge in visited)
                        push!(visited, next_edge)
                        push!(queue, next_edge)
                    end
                end
            end
            push!(components, component)
        end

        main_components = Set{Int}()
        for (i, component) in enumerate(components)
            neighbors = Set(Int32[local_neighbor[eid] for eid in component])
            if prev_boundary in neighbors || next_boundary in neighbors
                push!(main_components, i)
            end
        end

        for (i, component) in enumerate(components)
            i in main_components && continue
            seeds = Set(Int32[local_neighbor[eid] for eid in component if local_neighbor[eid] != prev_boundary && local_neighbor[eid] != next_boundary])
            isempty(seeds) && continue

            component_vertices = Set{Int32}()
            search = Int32[collect(seeds)...]
            valid = true
            while !isempty(search)
                current = popfirst!(search)
                current in component_vertices && continue
                if current == vertex
                    continue
                end
                if current in boundary_set
                    valid = false
                    break
                end
                push!(component_vertices, current)
                for next_vertex in adjacency[Int(current) + 1]
                    if next_vertex != vertex && !(next_vertex in component_vertices)
                        push!(search, next_vertex)
                    end
                end
            end

            valid || continue
            union!(to_remove, component_vertices)
        end
    end

    return sort!(collect(to_remove))
end

function _recover_boundary_cycle_from_triangles(triangles, triangle_edge_ids)
    triangle_edge_ids === nothing && return nothing
    tri = Int32.(triangles)
    tri_ids = Int32.(triangle_edge_ids)
    size(tri, 1) == size(tri_ids, 1) || throw(ArgumentError("triangle_edge_ids must align with triangle rows"))
    size(tri_ids, 2) == 3 || throw(ArgumentError("triangle_edge_ids must have three columns"))

    edge_pair = Dict{Int32,Tuple{Int32,Int32}}()
    incidence = Dict{Int32,Int}()
    for i in 1:size(tri, 1)
        a = tri[i, 1]
        b = tri[i, 2]
        c = tri[i, 3]
        pairs = (
            (min(a, b), max(a, b)),
            (min(b, c), max(b, c)),
            (min(c, a), max(c, a)),
        )
        for j in 1:3
            eid = tri_ids[i, j]
            pair = pairs[j]
            existing = get(edge_pair, eid, pair)
            existing == pair || return nothing
            edge_pair[eid] = pair
            incidence[eid] = get(incidence, eid, 0) + 1
        end
    end

    boundary_ids = sort!(Int32[eid for (eid, count) in incidence if count == 1])
    isempty(boundary_ids) && return nothing
    adjacency = Dict{Int32,Vector{Int32}}()
    for eid in boundary_ids
        u, v = edge_pair[eid]
        push!(get!(adjacency, u, Int32[]), v)
        push!(get!(adjacency, v, Int32[]), u)
    end
    all(length(nbrs) == 2 for nbrs in values(adjacency)) || return nothing

    start = minimum(collect(keys(adjacency)))
    order = Int32[start]
    prev = Int32(-1)
    current = start
    while true
        nbrs = adjacency[current]
        next_vertex = prev == Int32(-1) ? minimum(nbrs) : (nbrs[1] == prev ? nbrs[2] : nbrs[1])
        next_vertex == start && break
        push!(order, next_vertex)
        prev, current = current, next_vertex
        length(order) > length(adjacency) && return nothing
    end
    length(order) == length(adjacency) || return nothing
    return order
end

function _fk_cleanup_reindex_edges(edge_array, removed::Set{Int32}, old_to_new)
    arr = sanitize_edge_array(edge_array)
    keep = [i for i in 1:size(arr, 1) if !(arr[i, 1] in removed) && !(arr[i, 2] in removed)]
    isempty(keep) && return Matrix{Int32}(undef, 0, 2)
    return _drop_self_loops(_reindex_edges(arr[keep, :], old_to_new; deduplicate=false))
end

function _fk_cleanup_reindex_edge_groups(edge_groups, removed::Set{Int32}, old_to_new)
    next_groups = Dict{String,Matrix{Int32}}()
    for (name, arr) in pairs(edge_groups)
        remapped = _fk_cleanup_reindex_edges(arr, removed, old_to_new)
        size(remapped, 1) > 0 && (next_groups[string(name)] = remapped)
    end
    return next_groups
end

function _deduplicate_edge_pairs_preserve_order(edge_array)
    arr = sanitize_edge_array(edge_array)
    size(arr, 1) == 0 && return Matrix{Int32}(undef, 0, 2)

    pairs = Tuple{Int32,Int32}[]
    seen = Set{Tuple{Int32,Int32}}()
    for i in 1:size(arr, 1)
        u = arr[i, 1]
        v = arr[i, 2]
        u == v && continue
        pair = (min(u, v), max(u, v))
        pair in seen && continue
        push!(seen, pair)
        push!(pairs, pair)
    end

    out = Matrix{Int32}(undef, length(pairs), 2)
    for (i, (u, v)) in enumerate(pairs)
        out[i, 1] = u
        out[i, 2] = v
    end
    return out
end

function _canonical_triangle_edge_ids_by_pair(triangles, triangle_edge_ids)
    triangle_edge_ids === nothing && return nothing

    tri = Int32.(triangles)
    tri_ids = Int32.(triangle_edge_ids)
    size(tri, 1) == size(tri_ids, 1) || throw(ArgumentError("triangle_edge_ids must align with the triangle rows"))
    size(tri_ids, 2) == 3 || throw(ArgumentError("triangle_edge_ids must have three columns"))

    pair_to_id = Dict{Tuple{Int32,Int32},Int32}()
    next_ids = Matrix{Int32}(undef, size(tri_ids, 1), 3)

    for i in 1:size(tri, 1)
        a = tri[i, 1]
        b = tri[i, 2]
        c = tri[i, 3]
        pairs = (
            (min(a, b), max(a, b)),
            (min(b, c), max(b, c)),
            (min(c, a), max(c, a)),
        )
        for j in 1:3
            next_ids[i, j] = get!(pair_to_id, pairs[j], Int32(length(pair_to_id)))
        end
    end

    return next_ids
end

function _collapse_residual_parallel_edges(edges, edge_groups, triangles, triangle_edge_ids)
    collapsed_edges = _deduplicate_edge_pairs_preserve_order(edges)
    collapsed_groups = Dict{String,Matrix{Int32}}()
    for (name, arr) in pairs(edge_groups)
        remapped = _deduplicate_edge_pairs_preserve_order(arr)
        size(remapped, 1) > 0 && (collapsed_groups[string(name)] = remapped)
    end
    collapsed_triangle_edge_ids = _canonical_triangle_edge_ids_by_pair(triangles, triangle_edge_ids)
    return collapsed_edges, collapsed_groups, collapsed_triangle_edge_ids
end

function _triangle_supported_edge_pairs(triangles)
    tri = Int32.(triangles)
    allowed = Set{Tuple{Int32,Int32}}()
    for i in 1:size(tri, 1)
        a = tri[i, 1]
        b = tri[i, 2]
        c = tri[i, 3]
        push!(allowed, (min(a, b), max(a, b)))
        push!(allowed, (min(b, c), max(b, c)))
        push!(allowed, (min(c, a), max(c, a)))
    end
    return allowed
end

function _prune_to_triangle_support(num_vertices::Integer, edges, edge_groups, boundary_vertices, triangles, triangle_edge_ids, current_to_original)
    tri = Int32.(triangles)
    size(tri, 1) == 0 && return (
        Int(num_vertices),
        sanitize_edge_array(edges),
        Dict{String,Matrix{Int32}}(string(k) => sanitize_edge_array(v) for (k, v) in pairs(edge_groups)),
        Int32.(boundary_vertices),
        tri,
        triangle_edge_ids === nothing ? nothing : Int32.(triangle_edge_ids),
        Int32.(current_to_original),
        0,
        0,
    )

    current_n = Int(num_vertices)
    current_edges = sanitize_edge_array(edges)
    current_edge_groups = Dict{String,Matrix{Int32}}(string(k) => sanitize_edge_array(v) for (k, v) in pairs(edge_groups))
    current_boundary = Int32.(boundary_vertices)
    current_triangle_edge_ids = triangle_edge_ids === nothing ? nothing : Int32.(triangle_edge_ids)
    current_to_original = Int32.(current_to_original)
    edge_count_before = size(current_edges, 1)

    active_vertices = sort!(unique(vec(tri)))
    active_set = Set(active_vertices)
    removed_vertices = Int32[v for v in 0:(current_n - 1) if !(Int32(v) in active_set)]
    if !isempty(removed_vertices)
        removed_set = Set(removed_vertices)
        any(v -> v in removed_set, current_boundary) &&
            throw(ArgumentError("triangle-support pruning removed a boundary vertex from the FK circle-packing graph"))
        old_to_new = Dict(Int(v) => i - 1 for (i, v) in enumerate(active_vertices))
        current_edges = _fk_cleanup_reindex_edges(current_edges, removed_set, old_to_new)
        current_edge_groups = _fk_cleanup_reindex_edge_groups(current_edge_groups, removed_set, old_to_new)
        current_boundary = Int32[old_to_new[Int(v)] for v in current_boundary]
        tri = _reindex_triangles(tri, old_to_new; deduplicate=false)
        current_to_original = current_to_original[Int.(active_vertices) .+ 1]
        current_n = length(active_vertices)
    end

    return (
        current_n,
        current_edges,
        current_edge_groups,
        current_boundary,
        tri,
        current_triangle_edge_ids,
        current_to_original,
        length(removed_vertices),
        max(edge_count_before - size(current_edges, 1), 0),
    )
end

function _fk_parallel_strip_vertices(num_vertices::Integer, edges, boundary_vertices)
    components = _fk_parallel_strip_components(num_vertices, edges, boundary_vertices)
    isempty(components) && return Int32[]
    return sort!(unique(vcat([comp.vertices for comp in components]...)))
end

function _fk_parallel_strip_components(num_vertices::Integer, edges, boundary_vertices)
    current_edges = sanitize_edge_array(edges)
    size(current_edges, 1) == 0 && return NamedTuple{(:vertices, :attachment),Tuple{Vector{Int32},Tuple{Int32,Int32}}}[]

    n = Int(num_vertices)
    boundary_set = Set(Int32.(boundary_vertices))
    pair_multiplicity = Dict{Tuple{Int32,Int32},Int}()
    adjacency = [Set{Int32}() for _ in 1:n]

    for i in 1:size(current_edges, 1)
        u = current_edges[i, 1]
        v = current_edges[i, 2]
        key = (min(u, v), max(u, v))
        pair_multiplicity[key] = get(pair_multiplicity, key, 0) + 1
        push!(adjacency[Int(u)+1], v)
        push!(adjacency[Int(v)+1], u)
    end

    components = NamedTuple{(:vertices, :attachment),Tuple{Vector{Int32},Tuple{Int32,Int32}}}[]
    for ((u, v), multiplicity) in pair_multiplicity
        multiplicity < 2 && continue

        reduced_adjacency = [Set{Int32}() for _ in 1:n]
        for i in 1:size(current_edges, 1)
            a = current_edges[i, 1]
            b = current_edges[i, 2]
            if a == u || a == v || b == u || b == v
                continue
            end
            push!(reduced_adjacency[Int(a)+1], b)
            push!(reduced_adjacency[Int(b)+1], a)
        end

        seen = Set{Int32}([u, v])
        for start in Int32(0):Int32(n - 1)
            start in seen && continue
            queue = Int32[start]
            push!(seen, start)
            component = Set{Int32}()
            touches_boundary = false
            touches_u = false
            touches_v = false

            while !isempty(queue)
                current = popfirst!(queue)
                push!(component, current)
                current in boundary_set && (touches_boundary = true)
                (u in adjacency[Int(current)+1]) && (touches_u = true)
                (v in adjacency[Int(current)+1]) && (touches_v = true)

                for next_vertex in reduced_adjacency[Int(current)+1]
                    if !(next_vertex in seen)
                        push!(seen, next_vertex)
                        push!(queue, next_vertex)
                    end
                end
            end

            isempty(component) && continue
            (!touches_boundary && touches_u && touches_v) || continue
            push!(components, (vertices=sort!(collect(component)), attachment=(u, v)))
        end
    end

    return components
end

function _cleanup_fk_gasket_layout_graph(num_vertices::Integer, edges, edge_groups, boundary_vertices, triangles, triangle_edge_ids)
    current_n = Int(num_vertices)
    current_edges = sanitize_edge_array(edges)
    current_edge_groups = Dict{String,Matrix{Int32}}(string(k) => sanitize_edge_array(v) for (k, v) in pairs(edge_groups))
    current_boundary = Int32.(boundary_vertices)
    current_triangles = Int32.(triangles)
    current_triangle_edge_ids = triangle_edge_ids === nothing ? nothing : Int32.(triangle_edge_ids)
    current_to_original = Int32.(collect(0:(current_n - 1)))

    removed_bubble_total = Set{Int32}()
    removed_parallel_total = Set{Int32}()
    rounds = 0
    bubble_rounds = 0
    parallel_rounds = 0

    while true
        bubble_vertices = current_triangle_edge_ids === nothing ? Int32[] : _fk_gasket_bubble_vertices(
            current_n,
            current_edges,
            current_boundary,
            current_triangles,
            current_triangle_edge_ids,
        )
        boundary_sector_vertices = current_triangle_edge_ids === nothing ? Int32[] : _boundary_sector_bubble_vertices(
            current_n,
            current_edges,
            current_boundary,
            current_triangles,
            current_triangle_edge_ids,
        )
        parallel_components = _fk_parallel_strip_components(current_n, current_edges, current_boundary)
        parallel_vertices = isempty(parallel_components) ? Int32[] : sort!(unique(vcat([comp.vertices for comp in parallel_components]...)))

        removed_vertices = sort!(unique(vcat(bubble_vertices, boundary_sector_vertices, parallel_vertices)))
        isempty(removed_vertices) && break
        rounds += 1
        if !isempty(vcat(bubble_vertices, boundary_sector_vertices))
            bubble_rounds += 1
            union!(removed_bubble_total, bubble_vertices)
            union!(removed_bubble_total, boundary_sector_vertices)
        end
        if !isempty(parallel_vertices)
            parallel_rounds += 1
            union!(removed_parallel_total, parallel_vertices)
        end
        removed_set = Set(removed_vertices)

        keep_vertices = Int32[v for v in 0:(current_n-1) if !(Int32(v) in removed_set)]
        old_to_new = Dict(Int(v) => i - 1 for (i, v) in enumerate(keep_vertices))
        current_to_original = current_to_original[Int.(keep_vertices) .+ 1]

        current_edges = _fk_cleanup_reindex_edges(current_edges, removed_set, old_to_new)
        current_edge_groups = _fk_cleanup_reindex_edge_groups(current_edge_groups, removed_set, old_to_new)

        keep_triangles = [i for i in 1:size(current_triangles, 1) if !(current_triangles[i, 1] in removed_set) && !(current_triangles[i, 2] in removed_set) && !(current_triangles[i, 3] in removed_set)]
        current_triangle_edge_ids = current_triangle_edge_ids === nothing ? nothing : current_triangle_edge_ids[keep_triangles, :]
        current_triangles = _reindex_triangles(current_triangles[keep_triangles, :], old_to_new; deduplicate=false)
        current_boundary = Int32[old_to_new[Int(v)] for v in current_boundary if !(v in removed_set)]
        current_n = length(keep_vertices)
    end

    edge_count_before_collapse = size(current_edges, 1)
    current_edges, current_edge_groups, current_triangle_edge_ids = _collapse_residual_parallel_edges(
        current_edges,
        current_edge_groups,
        current_triangles,
        current_triangle_edge_ids,
    )
    current_n,
    current_edges,
    current_edge_groups,
    current_boundary,
    current_triangles,
    current_triangle_edge_ids,
    current_to_original,
    removed_unsupported_vertices,
    removed_unsupported_edges = _prune_to_triangle_support(
        current_n,
        current_edges,
        current_edge_groups,
        current_boundary,
        current_triangles,
        current_triangle_edge_ids,
        current_to_original,
    )
    recovered_boundary = _recover_boundary_cycle_from_triangles(current_triangles, current_triangle_edge_ids)
    recovered_boundary !== nothing && (current_boundary = Int32.(recovered_boundary))

    removed_total = union(removed_bubble_total, removed_parallel_total)
    return (
        current_n,
        current_edges,
        current_edge_groups,
        current_boundary,
        current_triangles,
        current_triangle_edge_ids,
        Dict{String,Any}(
            "gasket_removed_bubble_vertices" => length(removed_bubble_total),
            "gasket_removed_bubble_rounds" => bubble_rounds,
            "gasket_removed_parallel_strip_vertices" => length(removed_parallel_total),
            "gasket_removed_parallel_strip_rounds" => parallel_rounds,
            "gasket_geometric_cleanup" => rounds > 0,
            "circle_packing_removed_parallel_strip_vertices" => length(removed_parallel_total),
            "circle_packing_removed_parallel_strip_rounds" => parallel_rounds,
            "circle_packing_geometric_cleanup" => rounds > 0,
            "gasket_num_vertices" => current_n,
            "gasket_num_edges" => size(current_edges, 1),
            "gasket_num_collapsed_edges" => size(_deduplicate_undirected_edges(current_edges), 1),
            "gasket_collapsed_residual_parallel_edges" => max(edge_count_before_collapse - size(current_edges, 1), 0),
            "gasket_removed_nontriangulated_vertices" => removed_unsupported_vertices,
            "gasket_removed_nontriangulated_edges" => removed_unsupported_edges,
            "gasket_num_triangles" => size(current_triangles, 1),
            "gasket_num_collapsed_triangles" => size(sanitize_triangles(current_triangles; drop_degenerate=true, deduplicate=true), 1),
            "gasket_kept_vertex_indices" => Int.(current_to_original),
        ),
    )
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

function _drop_triangle_row(triangles, idx::Integer)
    tri = Int32.(triangles)
    size(tri, 1) == 0 && return tri
    keep = trues(size(tri, 1))
    keep[Int(idx) + 1] = false
    return tri[keep, :]
end

function _schnyder_outer_vertices(edge_groups::Dict{String,Matrix{Int32}})
    outer_edges = get(edge_groups, "outer", Matrix{Int32}(undef, 0, 2))
    outer_vertices = size(outer_edges, 1) == 0 ? Int32[0, 1, 2] : sort(unique(vec(outer_edges)))
    length(outer_vertices) == 3 || throw(ArgumentError("Schnyder layout expects exactly three outer vertices"))
    return outer_vertices
end

function _uniform_disk_face_data(map_data::UniformMap; outer_face_index=nothing)
    faces = [Int32.(collect(face)) for face in _faces_as_lists(map_data.faces)]
    kept_faces, kept_face_edge_ids, keep_mask, collapsed_quad_edges =
        _filter_collapsible_quad_faces_with_edge_ids(faces, map_data.face_edge_ids)
    choice = if outer_face_index === nothing
        _first_simple_face(faces, 4)
    else
        idx = Int(outer_face_index)
        0 <= idx < length(faces) || throw(ArgumentError("outer_face_index is out of range"))
        (idx, unique(faces[idx + 1]))
    end

    choice === nothing && throw(ArgumentError("failed to find a quadrilateral face with four distinct vertices for uniform circle packing"))
    face_idx, boundary = choice
    length(boundary) >= 4 || throw(ArgumentError("failed to find a quadrilateral face with four distinct vertices for uniform circle packing"))
    keep_mask[face_idx + 1] || throw(ArgumentError("outer_face_index refers to a collapsed degenerate quadrilateral; choose a simple quadrilateral face instead"))

    keep = [j1 for j1 in eachindex(faces) if (j1 - 1) != face_idx && keep_mask[j1]]
    interior_faces = [faces[j1] for j1 in keep]
    interior_face_edge_ids = [map_data.face_edge_ids[j1] for j1 in keep]
    triangles, triangle_edge_ids = _fan_triangulate_faces_with_face_edge_ids(interior_faces, interior_face_edge_ids; next_synthetic=num_edges(map_data), drop_degenerate=true)

    meta = Dict{String,Any}(
        "outer_face_removed" => true,
        "outer_face_index" => face_idx,
        "circle_packing_collapsed_degenerate_quads" => count(!, keep_mask),
        "circle_packing_collapsed_degenerate_quad_edges" => Int.(collapsed_quad_edges),
    )
    return Int32.(boundary[1:4]), interior_faces, triangles, triangle_edge_ids, meta
end

function _sphere_circle_packing_hc_face_data(map_data::FKMap)
    tri_faces, triangle_edge_ids = _surface_triangle_data(map_data)
    choice = _first_simple_face(tri_faces, 3)
    choice === nothing && throw(ArgumentError("failed to find a simple FK triangular face for sphere circle packing"))
    idx, boundary = choice
    length(boundary) >= 3 || throw(ArgumentError("failed to find a simple FK triangular face for sphere circle packing"))

    return Int32.(boundary[1:3]), _drop_triangle_row(tri_faces, idx), _drop_triangle_row(triangle_edge_ids, idx), Dict{String,Any}(
        "circle_packing_outer_triangle_index" => idx,
        "circle_packing_outer_boundary_vertices" => Int.(boundary[1:3]),
    )
end

function _vertex_link_order_from_triangle_edge_ids(num_vertices::Integer, triangles, triangle_edge_ids, vertex::Integer)
    n = Int(num_vertices)
    v = Int32(vertex)
    0 <= v < n || throw(ArgumentError("vertex must lie in [0, num_vertices)"))

    tri = sanitize_triangles(triangles; drop_degenerate=true, deduplicate=false)
    tri_edge = Int32.(triangle_edge_ids)
    size(tri_edge, 1) == size(tri, 1) || throw(ArgumentError("triangle_edge_ids must align with the triangle rows"))
    size(tri_edge, 2) == 3 || throw(ArgumentError("triangle_edge_ids must have three columns ordered as (ab, bc, ca)"))

    local_adj = Dict{Int32,Vector{Tuple{Int32,Int}}}()
    local_deg = Dict{Int32,Int}()
    local_neighbor = Dict{Int32,Int32}()
    local_wedge_count = 0

    function register_neighbor!(edge_id::Int32, neighbor::Int32)
        existing = get(local_neighbor, edge_id, Int32(-1))
        if existing != Int32(-1) && existing != neighbor
            throw(ArgumentError("triangle_edge_ids are inconsistent with the triangulation around vertex $vertex"))
        end
        local_neighbor[edge_id] = neighbor
    end

    function add_local_wedge!(left_edge::Int32, left_neighbor::Int32, right_edge::Int32, right_neighbor::Int32)
        register_neighbor!(left_edge, left_neighbor)
        register_neighbor!(right_edge, right_neighbor)
        local_wedge_count += 1
        wedge_id = local_wedge_count
        push!(get!(local_adj, left_edge, Tuple{Int32,Int}[]), (right_edge, wedge_id))
        push!(get!(local_adj, right_edge, Tuple{Int32,Int}[]), (left_edge, wedge_id))
        local_deg[left_edge] = get(local_deg, left_edge, 0) + 1
        local_deg[right_edge] = get(local_deg, right_edge, 0) + 1
    end

    for i in 1:size(tri, 1)
        a = tri[i, 1]
        b = tri[i, 2]
        c = tri[i, 3]
        eab = tri_edge[i, 1]
        ebc = tri_edge[i, 2]
        eca = tri_edge[i, 3]
        if a == v
            add_local_wedge!(eca, c, eab, b)
        elseif b == v
            add_local_wedge!(eab, a, ebc, c)
        elseif c == v
            add_local_wedge!(ebc, b, eca, a)
        end
    end

    isempty(local_adj) && return Int32[]

    for values in values(local_adj)
        sort!(values; by=x -> (x[1], x[2]))
    end

    odd_edges = sort!(Int32[eid for (eid, d) in local_deg if isodd(d)])
    visited_edges = Set{Int32}()
    component_count = 0
    for start in keys(local_adj)
        start in visited_edges && continue
        component_count += 1
        queue = Int32[start]
        push!(visited_edges, start)
        while !isempty(queue)
            current = popfirst!(queue)
            for (next_edge, _) in get(local_adj, current, Tuple{Int32,Int}[])
                if !(next_edge in visited_edges)
                    push!(visited_edges, next_edge)
                    push!(queue, next_edge)
                end
            end
        end
    end
    component_count == 1 || throw(ArgumentError("chosen FK outer vertex $vertex induces $component_count disconnected local sectors"))

    function eulerian_trail(start::Int32)
        wedge_ids = sort!(unique(Int[wedge_id for values in values(local_adj) for (_, wedge_id) in values]))
        used = falses(isempty(wedge_ids) ? 0 : maximum(wedge_ids))
        stack = Int32[start]
        path = Int32[]

        while !isempty(stack)
            current = stack[end]
            adjlist = get(local_adj, current, Tuple{Int32,Int}[])
            while !isempty(adjlist) && used[last(adjlist)[2]]
                pop!(adjlist)
            end

            if isempty(adjlist)
                push!(path, current)
                pop!(stack)
            else
                next_edge, wedge_id = pop!(adjlist)
                used[wedge_id] = true
                push!(stack, next_edge)
            end
        end

        all(used[wedge_ids]) || throw(ArgumentError("failed to recover the neighbor ordering around FK outer vertex $vertex"))
        reverse!(path)
        return path
    end

    if isempty(odd_edges)
        start_edge = minimum(collect(keys(local_adj)))
        edge_cycle = eulerian_trail(start_edge)
        length(edge_cycle) >= 2 && edge_cycle[1] == edge_cycle[end] || throw(ArgumentError("failed to close the FK outer-vertex link cycle"))
        order_edges = edge_cycle[1:end-1]
        length(order_edges) == local_wedge_count || throw(ArgumentError("FK outer-vertex link cycle has the wrong length"))
        return Int32[local_neighbor[eid] for eid in order_edges]
    elseif length(odd_edges) == 2 && local_neighbor[odd_edges[1]] == local_neighbor[odd_edges[2]]
        edge_path = eulerian_trail(odd_edges[1])
        edge_path[end] == odd_edges[2] || throw(ArgumentError("failed to recover the FK repeated-neighbor outer-vertex link"))
        length(edge_path) == local_wedge_count + 1 || throw(ArgumentError("FK outer-vertex link path has the wrong length"))
        return Int32[local_neighbor[eid] for eid in edge_path]
    end

    throw(ArgumentError("chosen FK outer vertex $vertex induces odd link endpoints $(Tuple(Int.(odd_edges)))"))
end

function _sphere_circle_packing_hc_vertex_candidate(
    map_data::FKMap,
    tri_faces,
    triangle_edge_ids,
    pack_edges_old,
    pack_green_pairs,
    candidate::Integer,
)
    n = num_vertices(map_data)
    order = _vertex_link_order_from_triangle_edge_ids(n, tri_faces, triangle_edge_ids, candidate)
    length(order) >= 3 || throw(ArgumentError("chosen FK outer vertex $candidate does not induce a large enough boundary cycle"))

    outer = Int32(candidate)
    keep_tri = BitVector(undef, size(tri_faces, 1))
    for i in 1:size(tri_faces, 1)
        keep_tri[i] = all(tri_faces[i, j] != outer for j in 1:3)
    end
    pack_triangles_old = tri_faces[keep_tri, :]
    candidate_triangle_edge_ids = triangle_edge_ids[keep_tri, :]
    candidate_keep_vertices = Int32[v for v in 0:(n - 1) if v != outer]
    old_to_new = Dict{Int,Int}(Int(v) => i - 1 for (i, v) in enumerate(candidate_keep_vertices))

    keep_edge = BitVector(undef, size(pack_edges_old, 1))
    for i in 1:size(pack_edges_old, 1)
        keep_edge[i] = (pack_edges_old[i, 1] != outer) && (pack_edges_old[i, 2] != outer)
    end
    candidate_edges = _drop_self_loops(_reindex_edges(pack_edges_old[keep_edge, :], old_to_new; deduplicate=false))
    candidate_triangles = _reindex_triangles(pack_triangles_old, old_to_new; deduplicate=false)
    candidate_boundary = Int32[Int32(old_to_new[Int(v)]) for v in order]

    return (
        outer=Int32(candidate),
        boundary_old=Int32.(order),
        keep_vertices=Int32.(candidate_keep_vertices),
        pack_edges=Int32.(candidate_edges),
        pack_boundary=Int32.(candidate_boundary),
        pack_triangles=Int32.(candidate_triangles),
        pack_triangle_edge_ids=Int32.(candidate_triangle_edge_ids),
        pack_green_pairs=Int32.(pack_green_pairs),
        boundary_length=length(order),
    )
end

function _sphere_circle_packing_hc_outer_vertex_candidate(map_data::FKMap, edge_groups; outer_vertex=nothing)
    n = num_vertices(map_data)
    tri_faces, triangle_edge_ids = _surface_triangle_data(map_data)
    pack_edges_old = active_edge_pairs(map_data; collapse=false, drop_loops=true)
    pack_green_pairs = Int32.(map_data.triangle_green_edges)
    candidates = outer_vertex === nothing ? collect(0:(n - 1)) : [Int(outer_vertex)]
    last_error = nothing

    for candidate in candidates
        try
            raw_candidate = _sphere_circle_packing_hc_vertex_candidate(
                map_data,
                tri_faces,
                triangle_edge_ids,
                pack_edges_old,
                pack_green_pairs,
                candidate,
            )
            return _finalize_hc_sphere_circle_packing_candidate(map_data, edge_groups, raw_candidate)
        catch err
            last_error = err
        end
    end

    if outer_vertex !== nothing
        msg = last_error === nothing ? "failed to find a valid FK outer vertex for sphere circle packing" : sprint(showerror, last_error)
        throw(ArgumentError("invalid FK layout.options.outer_vertex=$(outer_vertex): $msg. This must be an original removable FK vertex, not a rendered synthetic outer vertex index."))
    end
    throw(last_error === nothing ? ArgumentError("failed to find a valid FK outer vertex for sphere circle packing") : last_error)
end

function _finalize_hc_sphere_circle_packing_candidate(map_data::FKMap, edge_groups, candidate_data)
    keep_vertices = Int32.(candidate_data.keep_vertices)
    keep_set = Set(keep_vertices)
    removed_vertices = Set(Int32[v for v in 0:(num_vertices(map_data)-1) if !(Int32(v) in keep_set)])
    old_to_new = Dict(Int(v) => i - 1 for (i, v) in enumerate(keep_vertices))
    packing_edge_groups = _fk_cleanup_reindex_edge_groups(edge_groups, removed_vertices, old_to_new)

    packing_num_vertices, packing_edges, packing_edge_groups, packing_boundary_vertices, packing_triangles, packing_triangle_edge_ids, cleanup_extra =
        _cleanup_fk_gasket_layout_graph(
            length(keep_vertices),
            candidate_data.pack_edges,
            packing_edge_groups,
            candidate_data.pack_boundary,
            candidate_data.pack_triangles,
            candidate_data.pack_triangle_edge_ids,
        )

    metadata = Dict{String,Any}(
        "circle_packing_removed_outer_vertex" => Int(candidate_data.outer),
        "circle_packing_outer_boundary_vertices" => Int.(candidate_data.boundary_old),
        "circle_packing_packing_to_render_vertices" => Int.(keep_vertices),
        "circle_packing_render_vertex_count" => num_vertices(map_data),
        "circle_packing_outer_boundary_length" => Int(candidate_data.boundary_length),
        "exploration_green_edge_pairs" => Int32.(candidate_data.pack_green_pairs),
    )
    merge!(metadata, cleanup_extra)
    metadata["circle_packing_packing_boundary_vertices"] = Int.(packing_boundary_vertices)

    local_keep = Int32.(collect(get(cleanup_extra, "gasket_kept_vertex_indices", Int[])))
    if !isempty(local_keep)
        metadata["circle_packing_packing_to_render_vertices"] = Int.(keep_vertices[Int.(local_keep) .+ 1])
        metadata["circle_packing_render_vertex_count"] = num_vertices(map_data)
    end

    length(packing_boundary_vertices) >= 3 || throw(ArgumentError("chosen FK outer vertex $(candidate_data.outer) does not induce a large enough boundary cycle"))
    length(unique(packing_boundary_vertices)) == length(packing_boundary_vertices) || throw(ArgumentError("chosen FK outer vertex $(candidate_data.outer) does not induce a simple boundary cycle"))
    _validate_disk_triangulation(
        packing_num_vertices,
        packing_triangles,
        packing_boundary_vertices;
        triangle_edge_ids=packing_triangle_edge_ids,
    )

    return (
        packing_num_vertices=packing_num_vertices,
        packing_edges=packing_edges,
        packing_edge_groups=packing_edge_groups,
        packing_boundary_vertices=packing_boundary_vertices,
        packing_triangles=packing_triangles,
        packing_triangle_edge_ids=packing_triangle_edge_ids,
        metadata=metadata,
    )
end

function _prepare_schnyder_problem(map_data::SchnyderMap; dimension::Int, boundary_scale::Float64)
    edges = layout_edges(map_data; drop_loops=true)
    edge_groups = _edge_groups_for_map(map_data)
    faces = _faces_as_lists(map_data.faces)

    if dimension == 3
        outer_vertices = _schnyder_outer_vertices(edge_groups)
        outer_set = Set(outer_vertices)
        interior_faces = [face for face in faces if Set(face) != outer_set]
        return LayoutProblem(
            num_vertices=map_data.nverts,
            edges=edges,
            edge_groups=edge_groups,
            faces=faces,
            packing_boundary_vertices=outer_vertices,
            packing_faces=interior_faces,
            metadata=Dict(
                "model" => "schnyder",
                "dimension" => 3,
                "circle_packing_topology" => "sphere",
                "circle_packing_outer_boundary_vertices" => Int.(outer_vertices),
            ),
            render_map_data=map_data,
        )
    end

    outer_vertices = _schnyder_outer_vertices(edge_groups)
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

function _mated_crt_raw_edges(map_data::MatedCRTMap)
    return sanitize_edge_array(hcat(map_data.edge_u, map_data.edge_v))
end

function _mated_crt_reconstruct_keep_vertices(map_data::MatedCRTMap)
    r = Int(map_data.refinement)
    raw_n = fld(length(map_data.L) - 1, r)
    raw_n >= num_vertices(map_data) ||
        throw(ArgumentError("stored mated_crt path is inconsistent with the quotient vertex count"))

    edge_u, edge_v, edge_kind = _build_mated_crt_edges(map_data.L, map_data.R, raw_n, r)
    raw_faces, raw_face_edge_ids = _mated_crt_face_cycles(raw_n, edge_u, edge_v, edge_kind, Int32[])

    faces = Vector{Vector{Int32}}()
    face_edge_ids = Vector{Vector{Int32}}()
    for (face, ids) in zip(raw_faces, raw_face_edge_ids)
        split_faces, split_ids = _mated_crt_split_pinched_face_walk(face, ids)
        append!(faces, split_faces)
        append!(face_edge_ids, split_ids)
    end

    edge_u, edge_v, edge_kind, _, _ =
        _collapse_mated_crt_parallel_edges(edge_u, edge_v, edge_kind, faces, face_edge_ids)

    intervals = _mated_crt_outermost_intervals(map_data.topology, raw_n, edge_u, edge_v, edge_kind)
    removed = falses(raw_n)
    for (a, b) in intervals
        for v in (Int(a) + 1):(Int(b) - 1)
            removed[v + 1] = true
        end
    end

    keep_vertices = Int32[v for v in 0:(raw_n - 1) if !removed[v + 1]]
    length(keep_vertices) == num_vertices(map_data) ||
        throw(ArgumentError("reconstructed mated_crt quotient vertices do not match the stored sphere map"))
    return raw_n, keep_vertices
end

function _mated_crt_rotate_loop(values::Vector{Float64}, refinement::Integer, shift_vertices::Integer)
    r = Int(refinement)
    n = fld(length(values) - 1, r)
    n >= 1 || throw(ArgumentError("mated_crt loop rotation needs at least one cell"))
    total_steps = n * r
    shift = mod(Int(shift_vertices), n)
    shift == 0 && return copy(values)

    out = similar(values)
    for j in 0:total_steps
        src = mod(j + shift * r, total_steps)
        out[j + 1] = values[src + 1]
    end
    out .-= out[1]
    return out
end

function _mated_crt_forward_mid_min_table(values::Vector{Float64}, refinement::Integer)
    r = Int(refinement)
    n = fld(length(values) - 1, r)
    vertex_vals = Float64[values[p * r + 1] for p in 0:(n - 1)]
    cell_min = Float64[minimum(@view values[p * r + 1:(p + 1) * r + 1]) for p in 0:(n - 1)]

    table = fill(Inf, n, n)
    for a in 0:(n - 1)
        current = vertex_vals[a + 1]
        table[a + 1, a + 1] = current
        for step in 1:(n - 1)
            if step >= 2
                current = min(current, cell_min[mod(a + step - 1, n) + 1])
            end
            b = mod(a + step, n)
            table[a + 1, b + 1] = current
        end
    end
    return table, cell_min
end

function _mated_crt_choose_cyclic_cut(L::Vector{Float64}, R::Vector{Float64}, refinement::Integer)
    r = Int(refinement)
    n = fld(length(L) - 1, r)
    n <= 1 && return 0

    lower_mid, cell_min_L = _mated_crt_forward_mid_min_table(L, r)
    upper_mid, cell_min_R = _mated_crt_forward_mid_min_table(R, r)
    covered = falses(n)

    for a in 0:(n - 1)
        for step in 2:(n - 1)
            b = mod(a + step, n)
            lower_ok = max(cell_min_L[a + 1], cell_min_L[b + 1]) <= lower_mid[a + 1, b + 1] + 1.0e-10
            upper_ok = max(cell_min_R[a + 1], cell_min_R[b + 1]) <= upper_mid[a + 1, b + 1] + 1.0e-10
            if lower_ok && upper_ok
                for k in 1:(step - 1)
                    covered[mod(a + k, n) + 1] = true
                end
            end
        end
    end

    for shift in 1:(n - 1)
        !covered[shift + 1] && return shift
    end
    return 1
end

function _mated_crt_graph_only_quotient_edges(
    L::Vector{Float64},
    R::Vector{Float64},
    refinement::Integer,
)
    r = Int(refinement)
    raw_n = fld(length(L) - 1, r)
    edge_u, edge_v, edge_kind = _build_mated_crt_edges(L, R, raw_n, r)
    edge_u, edge_v, edge_kind, _, _ =
        _collapse_mated_crt_parallel_edges(edge_u, edge_v, edge_kind, Vector{Vector{Int32}}(), Vector{Vector{Int32}}())

    intervals = _mated_crt_outermost_intervals(:sphere, raw_n, edge_u, edge_v, edge_kind)
    removed = falses(raw_n)
    for (a, b) in intervals
        for v in (Int(a) + 1):(Int(b) - 1)
            removed[v + 1] = true
        end
    end

    keep_vertices = Int32[v for v in 0:(raw_n - 1) if !removed[v + 1]]
    old_to_new = Dict(Int(v) => i - 1 for (i, v) in enumerate(keep_vertices))

    remapped = NTuple{2,Int32}[]
    for i in eachindex(edge_u)
        removed[Int(edge_u[i]) + 1] && continue
        removed[Int(edge_v[i]) + 1] && continue
        u = Int32(old_to_new[Int(edge_u[i])])
        v = Int32(old_to_new[Int(edge_v[i])])
        u == v && continue
        push!(remapped, (u, v))
    end

    edge_array = Matrix{Int32}(undef, length(remapped), 2)
    for (i, (u, v)) in enumerate(remapped)
        edge_array[i, 1] = u
        edge_array[i, 2] = v
    end
    return keep_vertices, first(collapse_undirected_edges(edge_array; drop_loops=true))
end

function _mated_crt_sphere_sfdp_edges(map_data::MatedCRTMap)
    map_data.topology == :sphere ||
        throw(ArgumentError("sphere-topology sfdp layout edges are defined only for sphere topology"))

    raw_n, keep_old = _mated_crt_reconstruct_keep_vertices(map_data)
    shift = _mated_crt_choose_cyclic_cut(map_data.L, map_data.R, map_data.refinement)
    L_rot = _mated_crt_rotate_loop(map_data.L, map_data.refinement, shift)
    R_rot = _mated_crt_rotate_loop(map_data.R, map_data.refinement, shift)
    keep_rot, edges_rot = _mated_crt_graph_only_quotient_edges(L_rot, R_rot, map_data.refinement)

    keep_raw = Int32[mod(Int(v) + shift, raw_n) for v in keep_rot]
    sort(copy(keep_raw)) == sort(copy(keep_old)) ||
        throw(ArgumentError("cyclic mated_crt sphere quotient changed the retained raw vertices"))

    raw_to_old = Dict(Int(v) => Int32(i - 1) for (i, v) in enumerate(keep_old))
    remapped = Matrix{Int32}(undef, size(edges_rot, 1), 2)
    for i in 1:size(edges_rot, 1)
        remapped[i, 1] = raw_to_old[Int(keep_raw[Int(edges_rot[i, 1]) + 1])]
        remapped[i, 2] = raw_to_old[Int(keep_raw[Int(edges_rot[i, 2]) + 1])]
    end
    return sanitize_edge_array(remapped), shift
end

function _mated_crt_edge_groups(map_data::MatedCRTMap)
    generic = _mated_crt_raw_edges(map_data)
    upper_pairs = NTuple{2,Int32}[]
    lower_pairs = NTuple{2,Int32}[]

    for i in eachindex(map_data.edge_u)
        u = map_data.edge_u[i]
        v = map_data.edge_v[i]
        if map_data.edge_kind[i] == :upper
            push!(upper_pairs, (u, v))
        elseif map_data.edge_kind[i] == :lower
            push!(lower_pairs, (u, v))
        end
    end

    out = Dict{String,Matrix{Int32}}("generic" => generic)
    if !isempty(upper_pairs)
        arr = Matrix{Int32}(undef, length(upper_pairs), 2)
        for (i, (u, v)) in enumerate(upper_pairs)
            arr[i, 1] = u
            arr[i, 2] = v
        end
        out["red"] = arr
    end
    if !isempty(lower_pairs)
        arr = Matrix{Int32}(undef, length(lower_pairs), 2)
        for (i, (u, v)) in enumerate(lower_pairs)
            arr[i, 1] = u
            arr[i, 2] = v
        end
        out["blue"] = arr
    end
    return out
end

function _mated_crt_face_lists_and_edge_ids(map_data::MatedCRTMap; drop_outer::Bool=false)
    outer = Int(map_data.outer_face_index)
    faces = Vector{Vector{Int32}}()
    edge_ids = Vector{Vector{Int32}}()
    for i in eachindex(map_data.faces)
        drop_outer && outer >= 0 && i == outer + 1 && continue
        push!(faces, Int32.(map_data.faces[i]))
        push!(edge_ids, Int32.(map_data.face_edge_ids[i]))
    end
    return faces, edge_ids
end

function _mated_crt_triangle_matrices(map_data::MatedCRTMap; drop_outer::Bool=false)
    face_lists, edge_lists = _mated_crt_face_lists_and_edge_ids(map_data; drop_outer=drop_outer)
    triangles = Vector{Vector{Int32}}()
    triangle_edge_ids = Vector{Vector{Int32}}()
    for (face, ids) in zip(face_lists, edge_lists)
        length(face) == 3 || throw(ArgumentError("mated_crt internal faces must be triangles"))
        push!(triangles, face)
        push!(triangle_edge_ids, ids)
    end

    tri = Matrix{Int32}(undef, length(triangles), 3)
    tri_ids = Matrix{Int32}(undef, length(triangle_edge_ids), 3)
    for i in eachindex(triangles)
        tri[i, :] .= triangles[i]
        tri_ids[i, :] .= triangle_edge_ids[i]
    end
    return tri, tri_ids, face_lists
end

function _mated_crt_sphere_face_data(map_data::MatedCRTMap)
    tri_faces, triangle_edge_ids, face_lists = _mated_crt_triangle_matrices(map_data; drop_outer=false)
    choice = _first_simple_face(face_lists, 3)
    choice === nothing && throw(ArgumentError("failed to find a simple mated_crt triangular face for sphere layouts"))
    idx, boundary = choice
    length(boundary) == 3 || throw(ArgumentError("failed to find a simple mated_crt outer triangle"))

    keep = trues(length(face_lists))
    keep[Int(idx) + 1] = false
    return (
        boundary=Int32.(boundary[1:3]),
        faces=face_lists[keep],
        triangles=_drop_triangle_row(tri_faces, idx),
        triangle_edge_ids=_drop_triangle_row(triangle_edge_ids, idx),
        metadata=Dict{String,Any}(
            "circle_packing_outer_triangle_index" => idx,
            "circle_packing_outer_boundary_vertices" => Int.(boundary[1:3]),
        ),
    )
end

function _mated_crt_metadata(map_data::MatedCRTMap, dimension::Int)
    metadata = Dict{String,Any}(
        "model" => "mated_crt",
        "dimension" => dimension,
        "topology" => String(map_data.topology),
        "vertices" => num_vertices(map_data),
        "gamma" => map_data.gamma,
        "gamma_prime" => map_data.gamma_prime,
        "kappa" => map_data.kappa,
        "kappa_prime" => map_data.kappa_prime,
        "correlation" => map_data.brownian_correlation,
        "refinement" => map_data.refinement,
        "sampler" => map_data.sampler,
    )
    if map_data.topology == :sphere && map_data.sphere_layout_map !== nothing
        metadata["sphere_construction"] = occursin("approx", lowercase(map_data.sampler)) ?
            "approx_hybrid_cone_walk_spanning_tree" :
            "exact_cone_walk_spanning_tree"
        metadata["linearized_word_length"] = length(map_data.sphere_layout_map.word)
        metadata["linearized_num_productions"] = length(map_data.sphere_layout_map.production_steps)
    end
    return metadata
end

function _prepare_mated_crt_problem(
    map_data::MatedCRTMap;
    dimension::Int,
    boundary_scale::Float64,
    engine=nothing,
)
    topo = map_data.topology
    layout_engine = lowercase(strip(string(something(engine, topo == :disk ? "tutte" : (dimension == 3 ? "circle_packing" : "tutte")))))
    edges = _mated_crt_raw_edges(map_data)
    edge_groups = _mated_crt_edge_groups(map_data)

    if topo == :disk
        dimension == 2 || throw(ArgumentError("disk-topology mated_crt maps currently support only `dimension: 2`"))
        tri, tri_ids, face_lists = _mated_crt_triangle_matrices(map_data; drop_outer=true)
        boundary = copy(map_data.boundary_vertices)
        length(boundary) >= 3 || throw(ArgumentError("mated_crt disk topology needs at least three boundary vertices"))
        metadata = _mated_crt_metadata(map_data, 2)
        metadata["boundary_positions_mode"] = "harmonic_measure"
        metadata["outer_face_removed"] = true
        metadata["layout_num_edges"] = size(edges, 1)
        metadata["layout_num_collapsed_edges"] = size(first(collapse_undirected_edges(edges; drop_loops=true)), 1)
        return LayoutProblem(
            num_vertices=num_vertices(map_data),
            edges=edges,
            edge_groups=edge_groups,
            boundary_vertices=boundary,
            boundary_positions=nothing,
            faces=face_lists,
            surface_triangles=tri,
            surface_triangle_edge_ids=tri_ids,
            packing_boundary_vertices=layout_engine == "circle_packing" ? boundary : nothing,
            packing_faces=layout_engine == "circle_packing" ? face_lists : nothing,
            packing_triangles=layout_engine == "circle_packing" ? tri : nothing,
            packing_triangle_edge_ids=layout_engine == "circle_packing" ? tri_ids : nothing,
            metadata=metadata,
            render_map_data=map_data,
        )
    end

    topo == :sphere || throw(ArgumentError("unsupported mated_crt topology $(repr(topo))"))
    metadata = _mated_crt_metadata(map_data, dimension)
    base_layout_edges = size(edges, 1)
    base_collapsed_edges = size(first(collapse_undirected_edges(edges; drop_loops=true)), 1)
    metadata["layout_num_edges"] = base_layout_edges
    metadata["layout_num_collapsed_edges"] = base_collapsed_edges

    if map_data.sphere_layout_map !== nothing
        tri_all, tri_ids_all, face_lists_all = _mated_crt_triangle_matrices(map_data; drop_outer=false)
        if dimension == 3
            if layout_engine == "circle_packing"
                try
                    chosen = _sphere_circle_packing_hc_outer_vertex_candidate(map_data.sphere_layout_map, edge_groups)
                    merge!(metadata, chosen.metadata)
                    metadata["circle_packing_topology"] = "sphere"
                    metadata["circle_packing_outer_vertex_selection"] = "first_valid"
                    metadata["layout_graph_mode"] = "spanning_tree_raw"
                    return LayoutProblem(
                        num_vertices=chosen.packing_num_vertices,
                        edges=chosen.packing_edges,
                        edge_groups=chosen.packing_edge_groups,
                        surface_triangles=chosen.packing_triangles,
                        surface_triangle_edge_ids=chosen.packing_triangle_edge_ids,
                        packing_boundary_vertices=chosen.packing_boundary_vertices,
                        packing_triangles=chosen.packing_triangles,
                        packing_triangle_edge_ids=chosen.packing_triangle_edge_ids,
                        metadata=metadata,
                        render_map_data=map_data,
                    )
                catch err
                    if !(err isa ArgumentError)
                        rethrow()
                    end

                    sphere_face = _mated_crt_sphere_face_data(map_data)
                    merge!(metadata, sphere_face.metadata)
                    metadata["circle_packing_topology"] = "sphere"
                    metadata["circle_packing_outer_vertex_selection"] = "mated_crt_face_fallback"
                    metadata["layout_graph_mode"] = "spanning_tree_raw"
                    return LayoutProblem(
                        num_vertices=num_vertices(map_data),
                        edges=edges,
                        edge_groups=edge_groups,
                        faces=face_lists_all,
                        surface_triangles=tri_all,
                        surface_triangle_edge_ids=tri_ids_all,
                        packing_boundary_vertices=sphere_face.boundary,
                        packing_faces=sphere_face.faces,
                        packing_triangles=sphere_face.triangles,
                        packing_triangle_edge_ids=sphere_face.triangle_edge_ids,
                        metadata=metadata,
                        render_map_data=map_data,
                    )
                end
            end

            layout_engine == "sfdp" || throw(ArgumentError("sphere-topology mated_crt maps support `layout.engine: sfdp` or `circle_packing` in 3D"))
            metadata["layout_graph_mode"] = "spanning_tree_raw"
            return LayoutProblem(
                num_vertices=num_vertices(map_data),
                edges=edges,
                edge_groups=edge_groups,
                faces=face_lists_all,
                surface_triangles=tri_all,
                surface_triangle_edge_ids=tri_ids_all,
                metadata=metadata,
                render_map_data=map_data,
            )
        elseif dimension == 2
            sphere_face = _mated_crt_sphere_face_data(map_data)
            merge!(metadata, sphere_face.metadata)
            metadata["outer_face_removed"] = true
            metadata["boundary_positions_mode"] = "explicit_triangle"
            return LayoutProblem(
                num_vertices=num_vertices(map_data),
                edges=edges,
                edge_groups=edge_groups,
                boundary_vertices=sphere_face.boundary,
                boundary_positions=equilateral_triangle(boundary_scale),
                faces=sphere_face.faces,
                surface_triangles=sphere_face.triangles,
                surface_triangle_edge_ids=sphere_face.triangle_edge_ids,
                packing_boundary_vertices=layout_engine == "circle_packing" ? sphere_face.boundary : nothing,
                packing_faces=layout_engine == "circle_packing" ? sphere_face.faces : nothing,
                packing_triangles=layout_engine == "circle_packing" ? sphere_face.triangles : nothing,
                packing_triangle_edge_ids=layout_engine == "circle_packing" ? sphere_face.triangle_edge_ids : nothing,
                metadata=metadata,
                render_map_data=map_data,
            )
        end

        throw(ArgumentError("dimension must be 2 or 3"))
    end

    tri_all, tri_ids_all, face_lists_all = _mated_crt_triangle_matrices(map_data; drop_outer=false)
    sphere_face = _mated_crt_sphere_face_data(map_data)

    merge!(metadata, sphere_face.metadata)
    metadata["circle_packing_preserves_render_vertices"] = true

    if dimension == 3
        if layout_engine == "circle_packing"
            metadata["circle_packing_topology"] = "sphere"
            return LayoutProblem(
                num_vertices=num_vertices(map_data),
                edges=edges,
                edge_groups=edge_groups,
                faces=face_lists_all,
                surface_triangles=tri_all,
                surface_triangle_edge_ids=tri_ids_all,
                packing_boundary_vertices=sphere_face.boundary,
                packing_faces=sphere_face.faces,
                packing_triangles=sphere_face.triangles,
                packing_triangle_edge_ids=sphere_face.triangle_edge_ids,
                metadata=metadata,
                render_map_data=map_data,
            )
        end

        layout_engine == "sfdp" || throw(ArgumentError("sphere-topology mated_crt maps support `layout.engine: sfdp` or `circle_packing` in 3D"))
        sfdp_edges, cut_shift = _mated_crt_sphere_sfdp_edges(map_data)
        metadata["layout_graph_mode"] = "sphere_cyclic_quotient"
        metadata["layout_num_base_edges"] = base_layout_edges
        metadata["layout_num_base_collapsed_edges"] = base_collapsed_edges
        metadata["layout_graph_cut_shift"] = cut_shift
        metadata["layout_num_edges"] = size(sfdp_edges, 1)
        metadata["layout_num_collapsed_edges"] = size(first(collapse_undirected_edges(sfdp_edges; drop_loops=true)), 1)
        return LayoutProblem(
            num_vertices=num_vertices(map_data),
            edges=sfdp_edges,
            edge_groups=edge_groups,
            faces=face_lists_all,
            surface_triangles=tri_all,
            surface_triangle_edge_ids=tri_ids_all,
            metadata=metadata,
            render_map_data=map_data,
        )
    elseif dimension == 2
        metadata["outer_face_removed"] = true
        metadata["boundary_positions_mode"] = "explicit_triangle"
        return LayoutProblem(
            num_vertices=num_vertices(map_data),
            edges=edges,
            edge_groups=edge_groups,
            boundary_vertices=sphere_face.boundary,
            boundary_positions=equilateral_triangle(boundary_scale),
            faces=sphere_face.faces,
            surface_triangles=sphere_face.triangles,
            surface_triangle_edge_ids=sphere_face.triangle_edge_ids,
            packing_boundary_vertices=layout_engine == "circle_packing" ? sphere_face.boundary : nothing,
            packing_faces=layout_engine == "circle_packing" ? sphere_face.faces : nothing,
            packing_triangles=layout_engine == "circle_packing" ? sphere_face.triangles : nothing,
            packing_triangle_edge_ids=layout_engine == "circle_packing" ? sphere_face.triangle_edge_ids : nothing,
            metadata=metadata,
            render_map_data=map_data,
        )
    end

    throw(ArgumentError("dimension must be 2 or 3"))
end

function _prepare_uniform_problem(map_data::UniformMap; dimension::Int, boundary_scale::Float64, outer_face_index=nothing, engine=nothing)
    layout_engine = lowercase(strip(string(something(engine, dimension == 3 ? "sfdp" : "tutte"))))

    if layout_engine == "circle_packing"
        throw(ArgumentError("uniform maps do not currently support circle_packing; use `tutte` in 2D or `sfdp` in 3D"))
    end

    edges = layout_edges(map_data; drop_loops=true)
    edge_groups = _edge_groups_for_map(map_data)
    faces = _faces_as_lists(map_data.faces)

    if dimension == 3
        raw_edges = edge_pairs(map_data; collapse=false, drop_loops=true)
        return LayoutProblem(
            num_vertices=num_vertices(map_data),
            edges=raw_edges,
            edge_groups=Dict("generic" => raw_edges),
            faces=faces,
            metadata=Dict(
                "model" => "uniform",
                "dimension" => 3,
                "layout_num_edges" => size(raw_edges, 1),
                "layout_num_collapsed_edges" => size(first(collapse_undirected_edges(raw_edges; drop_loops=true)), 1),
            ),
            render_map_data=map_data,
        )
    end

    boundary_vertices4, interior_faces, interior_triangles, interior_triangle_edge_ids, disk_meta =
        _uniform_disk_face_data(map_data; outer_face_index=outer_face_index)

    return LayoutProblem(
        num_vertices=num_vertices(map_data),
        edges=edges,
        edge_groups=edge_groups,
        boundary_vertices=boundary_vertices4,
        boundary_positions=square(boundary_scale),
        faces=interior_faces,
        surface_triangles=interior_triangles,
        surface_triangle_edge_ids=interior_triangle_edge_ids,
        metadata=merge(Dict(
            "model" => "uniform",
            "dimension" => 2,
            "outer_face_removed" => true,
            "outer_face_index" => disk_meta["outer_face_index"],
        ), disk_meta),
        render_map_data=map_data,
    )
end

function _prepare_half_plane_meandric_problem(
    map_data::HalfPlaneMeandricSystemMap;
    dimension::Int,
    engine=nothing,
)
    dimension == 2 || throw(ArgumentError("half-plane meandric systems currently support only `dimension: 2` with a Tutte embedding"))

    layout_engine = lowercase(strip(string(something(engine, "tutte"))))
    layout_engine == "tutte" || throw(ArgumentError("half-plane meandric systems currently support only `layout.engine: tutte`"))

    return LayoutProblem(
        num_vertices=num_vertices(map_data),
        edges=layout_edges(map_data; drop_loops=true),
        edge_groups=grouped_edges(map_data),
        boundary_vertices=copy(map_data.boundary_vertices),
        metadata=Dict(
            "model" => "half_plane_meandric",
            "dimension" => 2,
            "boundary_half_length" => Int(map_data.boundary_half_length),
            "boundary_position_source" => "harmonic_measure",
            "raw_component_sizes" => Int.(map_data.raw_component_sizes),
            "linked_component_sizes" => Int.(map_data.linked_component_sizes),
            "cle_component_sizes" => Int.(map_data.cle_component_sizes),
            "linked_interface" => map_data.linked_interface === nothing ? nothing : Int[map_data.linked_interface[1], map_data.linked_interface[2]],
        ),
        render_map_data=map_data,
    )
end

function _prepare_closed_meandric_problem(
    map_data::UniformMeandricSystemMap;
    dimension::Int,
    engine=nothing,
)
    dimension == 3 || throw(ArgumentError("closed meandric systems currently support only `dimension: 3`; use `layout.engine: sfdp`"))

    layout_engine = lowercase(strip(string(something(engine, "sfdp"))))
    layout_engine == "sfdp" || throw(ArgumentError("closed meandric systems currently support only `layout.engine: sfdp`; circle packing needs a sphere triangulation that is not wired in yet"))

    faces = [Int32.(collect(face)) for face in map_data.faces]
    triangles, triangle_edge_ids = _fan_triangulate_faces_with_face_edge_ids(
        faces,
        map_data.face_edge_ids;
        next_synthetic=length(map_data.face_edge_group),
        drop_degenerate=true,
    )
    layout_edges_augmented, layout_meta = _augment_closed_meandric_sfdp_graph(map_data, triangles)
    metadata = Dict{String,Any}(
        "model" => _meandric_model_name(map_data),
        "dimension" => 3,
        "num_loops" => num_loops(map_data),
        "component_sizes" => Int.(map_data.component_sizes),
        "num_faces" => length(faces),
        "upper_face_count" => Int(map_data.order) + 1,
        "lower_face_count" => Int(map_data.order) + 1,
    )
    merge!(metadata, layout_meta)
    merge!(metadata, map_data.sampler_metadata)

    return LayoutProblem(
        num_vertices=Int(map_data.auxiliary_num_vertices),
        edges=layout_edges_augmented,
        edge_groups=grouped_edges(map_data),
        faces=faces,
        surface_triangles=triangles,
        surface_triangle_edge_ids=triangle_edge_ids,
        metadata=metadata,
        render_map_data=map_data,
        render_vertex_indices=copy(map_data.render_vertex_indices),
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

function _hc_gasket_subgraph_raw(map_data::FKMap, gasket::AbstractString; metadata_override=nothing)
    kind, prod_idx, start_step, end_step, _, boundary_color, fresh_color = _hc_gasket_selection(map_data, gasket)
    boundary_vertices_old, meta = _hc_gasket_boundary(map_data, gasket)
    metadata_override === nothing || merge!(meta, metadata_override)

    root_tree = _build_current_gasket_tree(map_data, kind[1], start_step, end_step)
    get!(meta, "gasket_tree_word", _current_gasket_word(root_tree))

    keep_steps = [s for s in start_step:end_step if s != Int(map_data.order_face[prod_idx+1])]
    tri_faces_all = Int32.(map_data.triangulation_faces)
    tri_green_all = Int32.(map_data.triangle_green_edges)
    tri_diag_all = Int32.(map_data.triangle_diag_edge)

    selected_triangles_old = tri_faces_all[keep_steps.+1, :]
    selected_green = tri_green_all[keep_steps.+1, :]
    selected_diag = tri_diag_all[keep_steps.+1]

    valid_triangle_rows, drop_edge_ids = _gasket_cleanup_mask(
        map_data,
        selected_triangles_old,
        selected_green,
        selected_diag,
    )
    raw_triangle_count = size(selected_triangles_old, 1)
    raw_edge_ids = unique(Int32[v for v in vcat(vec(selected_green), selected_diag) if v >= 0])

    selected_triangles_old = selected_triangles_old[valid_triangle_rows, :]
    selected_green = selected_green[valid_triangle_rows, :]
    selected_diag = selected_diag[valid_triangle_rows]

    edge_ids = Int32[e for e in unique(Int32[v for v in vcat(vec(selected_green), selected_diag) if v >= 0]) if !(e in drop_edge_ids)]
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

    triangles = _reindex_triangles(selected_triangles_old, old_to_new; deduplicate=false)
    boundary_vertices = Int32[old_to_new[Int(v)] for v in boundary_vertices_old]

    edge_groups = Dict{String,Matrix{Int32}}()
    color_names = Dict(0 => "green", 1 => "red", 2 => "blue", 3 => "purple", 4 => "orange")
    for (color_id, name) in color_names
        mask = edge_color[edge_ids.+1] .== color_id
        any(mask) || continue
        pairs_old = hcat(edge_u[edge_ids[mask].+1], edge_v[edge_ids[mask].+1])
        pairs_new = _drop_self_loops(_reindex_edges(pairs_old, old_to_new; deduplicate=false))
        if size(pairs_new, 1) > 0
            edge_groups[name] = pairs_new
        end
    end

    all_edges = if isempty(edge_groups)
        Matrix{Int32}(undef, 0, 2)
    else
        sanitize_edge_array(reduce(vcat, collect(values(edge_groups))))
    end

    meta["boundary_mode"] = "$(kind)_gasket"
    meta["gasket_removed_exterior_triangle_step"] = Int(map_data.order_face[prod_idx+1])
    meta["gasket_original_vertex_ids"] = Int.(original_vertices)
    meta["gasket_num_vertices"] = length(original_vertices)
    meta["gasket_num_raw_edges"] = length(raw_edge_ids)
    meta["gasket_num_edges"] = size(all_edges, 1)
    meta["gasket_removed_degenerate_triangles"] = raw_triangle_count - size(selected_triangles_old, 1)
    meta["gasket_removed_degenerate_edges"] = length(raw_edge_ids) - length(edge_ids)
    meta["gasket_num_collapsed_edges"] = size(_deduplicate_undirected_edges(all_edges), 1)
    meta["gasket_num_raw_triangles"] = raw_triangle_count
    meta["gasket_num_triangles"] = size(triangles, 1)
    meta["gasket_num_collapsed_triangles"] = size(sanitize_triangles(triangles; drop_degenerate=true, deduplicate=true), 1)
    get!(meta, "gasket_clean_tree_word", meta["gasket_tree_word"])
    meta["gasket_boundary_vertex_color"] = boundary_color
    meta["gasket_boundary_edge_colors"] = [boundary_color, fresh_color]
    meta["exploration_green_edge_pairs"] = copy(selected_green)

    triangle_edge_ids = _triangle_side_edge_ids(selected_triangles_old, selected_green, selected_diag, edge_u, edge_v)

    return length(original_vertices), all_edges, edge_groups, boundary_vertices, triangles, triangle_edge_ids, meta, map_data
end

function _hc_gasket_subgraph(map_data::FKMap, gasket::AbstractString)
    kind, _, start_step, end_step, _, _, _ = _hc_gasket_selection(map_data, gasket)
    root_tree = _build_current_gasket_tree(map_data, kind[1], start_step, end_step)
    small_size_one, small_size_two = _count_small_current_gasket_holes(root_tree)
    raw_subword = map_data.word[start_step+1:end_step+1]
    cleaned_subword = _clean_current_gasket_fk_subword(map_data, root_tree)

    if small_size_one + small_size_two > 0 || cleaned_subword != raw_subword
        cleaned_tree = _clean_current_gasket_root(root_tree)
        decorated_tree = _decorate_current_gasket_tree(map_data, cleaned_tree)
        metadata_override = Dict{String,Any}(
            "gasket_tree_word" => _current_gasket_word(root_tree),
            "gasket_removed_size1_holes" => small_size_one,
            "gasket_removed_size2_holes" => small_size_two,
            "gasket_clean_tree_word" => _current_gasket_word(cleaned_tree),
            "gasket_clean_decorated_tree_word" => _current_gasket_word(decorated_tree),
            "gasket_clean_fk_word" => cleaned_subword,
        )
        try
            cleaned_map, reduced, balanced_word = _build_cleaned_fk_gasket_map(cleaned_subword)
            metadata_override["gasket_clean_reduced_word"] = reduced
            metadata_override["gasket_clean_balanced_word"] = balanced_word
            metadata_override["gasket_cleanup_mode"] = "fk_word_small_holes"
            metadata_override["gasket_boundary_follows_clean_word"] = true
            return _hc_gasket_subgraph_raw(cleaned_map, gasket; metadata_override=metadata_override)
        catch err
            metadata_override["gasket_clean_word_error"] = sprint(showerror, err)
            metadata_override["gasket_clean_word_fallback"] = true
            try
                nverts, gasket_edges, gasket_edge_groups, boundary_vertices, _, triangles, extra = _cleaned_current_gasket_layout(
                    kind[1],
                    decorated_tree,
                    kind == "h" ? "red" : "blue",
                    kind == "h" ? "orange" : "purple",
                )
                if length(boundary_vertices) >= 3 && size(triangles, 1) > 0
                    merge!(metadata_override, extra)
                    metadata_override["boundary_mode"] = "$(kind)_gasket"
                    metadata_override["gasket_boundary_vertex_color"] = kind == "h" ? "red" : "blue"
                    metadata_override["gasket_boundary_edge_colors"] = kind == "h" ? ["red", "orange"] : ["blue", "purple"]
                    metadata_override["gasket_boundary_follows_clean_tree"] = true
                    metadata_override["gasket_original_vertex_ids"] = Int[]
                    metadata_override["gasket_num_collapsed_edges"] = size(_deduplicate_undirected_edges(gasket_edges), 1)
                    metadata_override["gasket_num_collapsed_triangles"] = size(sanitize_triangles(triangles; drop_degenerate=true, deduplicate=true), 1)
                    get!(metadata_override, "gasket_tree_word", _current_gasket_word(root_tree))
                    return nverts, gasket_edges, gasket_edge_groups, boundary_vertices, triangles, nothing, metadata_override, nothing
                end
            catch tree_err
                metadata_override["gasket_clean_tree_error"] = sprint(showerror, tree_err)
            end
            metadata_override["gasket_clean_tree_fallback"] = true
            return _hc_gasket_subgraph_raw(map_data, gasket; metadata_override=metadata_override)
        end
    end

    metadata_override = Dict{String,Any}(
        "gasket_removed_size1_holes" => 0,
        "gasket_removed_size2_holes" => 0,
    )
    return _hc_gasket_subgraph_raw(map_data, gasket; metadata_override=metadata_override)
end

function _prepare_hc_problem(
    map_data::FKMap;
    dimension::Int,
    boundary_scale::Float64,
    hc_boundary_mode=nothing,
    engine=nothing,
    outer_vertex=nothing,
)
    edges = active_edge_pairs(map_data; collapse=false, drop_loops=true)
    edge_groups = _edge_groups_for_map(map_data)
    map_variant = _hc_variant(map_data)
    metadata = Dict{String,Any}("model" => "fk", "dimension" => dimension, "variant" => map_variant)
    metadata["layout_num_edges"] = size(edges, 1)
    metadata["layout_num_collapsed_edges"] = size(first(collapse_undirected_edges(edges; drop_loops=true)), 1)

    if dimension == 3
        packing_boundary_vertices = nothing
        packing_triangles = nothing
        packing_triangle_edge_ids = nothing
        packing_num_vertices = num_vertices(map_data)
        packing_edges = edges
        packing_edge_groups = edge_groups
        layout_engine = lowercase(strip(string(something(engine, "sfdp"))))
        if layout_engine == "circle_packing"
            chosen = _sphere_circle_packing_hc_outer_vertex_candidate(map_data, edge_groups; outer_vertex=outer_vertex)
            packing_num_vertices = chosen.packing_num_vertices
            packing_edges = chosen.packing_edges
            packing_edge_groups = chosen.packing_edge_groups
            packing_boundary_vertices = chosen.packing_boundary_vertices
            packing_triangles = chosen.packing_triangles
            packing_triangle_edge_ids = chosen.packing_triangle_edge_ids
            merge!(metadata, chosen.metadata)
            metadata["circle_packing_topology"] = "sphere"
            metadata["circle_packing_outer_vertex_selection"] = outer_vertex === nothing ? "first_valid" : "explicit"
        end
        surface_triangles, surface_triangle_edge_ids = layout_engine == "circle_packing" ? (packing_triangles, packing_triangle_edge_ids) : _surface_triangle_data(map_data)
        return LayoutProblem(
            num_vertices=packing_num_vertices,
            edges=packing_edges,
            edge_groups=packing_edge_groups,
            surface_triangles=surface_triangles,
            surface_triangle_edge_ids=surface_triangle_edge_ids,
            packing_boundary_vertices=packing_boundary_vertices,
            packing_triangles=packing_triangles,
            packing_triangle_edge_ids=packing_triangle_edge_ids,
            metadata=metadata,
            render_map_data=map_data,
        )
    end

    map_variant == "spanning_tree" &&
        throw(ArgumentError("spanning_tree maps do not currently support 2D embeddings; use `dimension: 3` with `sfdp` or `circle_packing`"))

    boundary_mode = lowercase(strip(string(something(hc_boundary_mode, "h_gasket"))))

    if boundary_mode in ("h", "c", "h_gasket", "c_gasket")
        nverts, gasket_edges, gasket_edge_groups, boundary_vertices, triangles, triangle_edge_ids, extra, gasket_render_map = _hc_gasket_subgraph(map_data, boundary_mode)
        nverts, gasket_edges, gasket_edge_groups, boundary_vertices, triangles, triangle_edge_ids, cleanup_extra =
            _cleanup_fk_gasket_layout_graph(nverts, gasket_edges, gasket_edge_groups, boundary_vertices, triangles, triangle_edge_ids)
        merge!(metadata, extra)
        merge!(metadata, cleanup_extra)
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
            surface_triangle_edge_ids=triangle_edge_ids,
            metadata=metadata,
            render_map_data=gasket_render_map,
        )
    end

    throw(ArgumentError("FK 2D layouts only support `hc_boundary_mode: h_gasket` or `hc_boundary_mode: c_gasket`"))
end

function prepare_layout_problem(map_data; dimension::Integer, boundary_scale=nothing, options=Dict{String,Any}())
    dim = Int(dimension)
    dim in (2, 3) || throw(ArgumentError("dimension must be 2 or 3"))
    scale = boundary_scale === nothing ? max(sqrt(max(num_vertices(map_data), 1)), 2.0) : float(boundary_scale)
    opts = Dict{String,Any}(string(k) => v for (k, v) in pairs(options))
    obsolete_circle_packing_options = ("outer_vertex_selection", "outer_vertex_preview_candidates", "outer_vertex_preview_maxiter", "outer_vertex_preview_tol")
    if get(opts, "engine", nothing) == "circle_packing"
        for key in obsolete_circle_packing_options
            haskey(opts, key) || continue
            throw(ArgumentError("layout.options.$key is no longer supported"))
        end
    end
    if map_data isa UniformMap && haskey(opts, "outer_face_triangle")
        throw(ArgumentError("uniform layouts do not support layout.options.outer_face_triangle"))
    end

    if map_data isa SchnyderMap
        return _prepare_schnyder_problem(map_data; dimension=dim, boundary_scale=scale)
    elseif map_data isa FKMap
        return _prepare_hc_problem(
            map_data;
            dimension=dim,
            boundary_scale=scale,
            hc_boundary_mode=get(opts, "hc_boundary_mode", nothing),
            engine=get(opts, "engine", nothing),
            outer_vertex=get(opts, "outer_vertex", nothing),
        )
    elseif map_data isa MatedCRTMap
        return _prepare_mated_crt_problem(
            map_data;
            dimension=dim,
            boundary_scale=scale,
            engine=get(opts, "engine", nothing),
        )
    elseif map_data isa UniformMap
        return _prepare_uniform_problem(
            map_data;
            dimension=dim,
            boundary_scale=scale,
            outer_face_index=get(opts, "outer_face_index", nothing),
            engine=get(opts, "engine", nothing),
        )
    elseif map_data isa HalfPlaneMeandricSystemMap
        return _prepare_half_plane_meandric_problem(
            map_data;
            dimension=dim,
            engine=get(opts, "engine", nothing),
        )
    elseif map_data isa UniformMeandricSystemMap
        return _prepare_closed_meandric_problem(
            map_data;
            dimension=dim,
            engine=get(opts, "engine", nothing),
        )
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

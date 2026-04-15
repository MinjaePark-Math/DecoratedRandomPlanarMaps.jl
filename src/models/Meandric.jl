# Half-plane and closed meandric-system models.

struct HalfPlaneMeandricSystemMap <: AbstractRandomPlanarMap
    order::Int32
    boundary_half_length::Int32
    upper_steps::Vector{Int8}
    lower_steps::Vector{Int8}
    upper_adj::Vector{Int32}            # 0-based partner, -1 if unmatched
    lower_adj::Vector{Int32}            # 0-based partner, -1 if unmatched
    boundary_vertices::Vector{Int32}    # unmatched lower openings, in boundary order
    edge_u::Vector{Int32}               # render graph: upper/lower arcs only
    edge_v::Vector{Int32}
    edge_group::Vector{String}
    layout_edge_u::Vector{Int32}        # layout graph: contour cycle + arcs
    layout_edge_v::Vector{Int32}
    raw_component_sizes::Vector{Int32}
    linked_component_sizes::Vector{Int32}
    cle_component_sizes::Vector{Int32}
    linked_interface::Union{Nothing,NTuple{2,Int32}}
end

const HalfPlaneMeandricMap = HalfPlaneMeandricSystemMap

struct UniformMeandricSystemMap <: AbstractRandomPlanarMap
    variant::Symbol
    order::Int32
    upper_steps::Vector{Int8}
    lower_steps::Vector{Int8}
    upper_adj::Vector{Int32}
    lower_adj::Vector{Int32}
    edge_u::Vector{Int32}               # render graph on the original strand vertices
    edge_v::Vector{Int32}
    edge_group::Vector{String}
    faces::Vector{Vector{Int32}}
    face_edge_ids::Vector{Vector{Int32}}
    face_edge_group::Vector{String}
    component_sizes::Vector{Int32}
    auxiliary_num_vertices::Int32       # glued-tree graph used for SFDP
    auxiliary_edge_u::Vector{Int32}
    auxiliary_edge_v::Vector{Int32}
    auxiliary_edge_group::Vector{String}
    render_vertex_indices::Vector{Int32}
    sampler_metadata::Dict{String,Any}
end

const UniformMeandricMap = UniformMeandricSystemMap
const UniformMeanderMap = UniformMeandricSystemMap

mutable struct _DyckSwapIndexPool
    indices::Vector{Int32}
    slot_of::Vector{Int32}
end

num_vertices(m::HalfPlaneMeandricSystemMap) = length(m.upper_adj)
num_faces(m::HalfPlaneMeandricSystemMap) = Int(m.order)

function layout_edges(m::HalfPlaneMeandricSystemMap; drop_loops::Bool=true)
    pairs = hcat(m.layout_edge_u, m.layout_edge_v)
    return first(collapse_undirected_edges(pairs; drop_loops=drop_loops))
end

num_vertices(m::UniformMeandricSystemMap) = length(m.upper_adj)
num_faces(m::UniformMeandricSystemMap) = length(m.faces)
num_edges(m::UniformMeandricSystemMap) = length(m.edge_u)
num_loops(m::UniformMeandricSystemMap) = length(m.component_sizes)

function layout_edges(m::UniformMeandricSystemMap; drop_loops::Bool=true)
    pairs = hcat(m.edge_u, m.edge_v)
    return first(collapse_undirected_edges(pairs; drop_loops=drop_loops))
end

function auxiliary_layout_edges(m::UniformMeandricSystemMap; drop_loops::Bool=true)
    pairs = hcat(m.auxiliary_edge_u, m.auxiliary_edge_v)
    return first(collapse_undirected_edges(pairs; drop_loops=drop_loops))
end

is_meander(m::UniformMeandricSystemMap) = m.variant == :uniform_meander

function _meandric_model_name(m::UniformMeandricSystemMap)
    return is_meander(m) ? "uniform_meander" : "uniform_meandric"
end

function sample_uniform_dyck_path_with_adjacency(order::Integer, rng::AbstractRNG=Random.default_rng())
    n = Int(order)
    n >= 1 || throw(ArgumentError("order must be >= 1"))

    steps = Vector{Int8}(undef, 2n)
    adj = fill(Int32(-1), 2n)
    open_stack = Int32[]
    sizehint!(open_stack, n)

    opens = 0
    closes = 0
    for i in 1:(2n)
        total = (2n - opens - closes) * (opens - closes + 1)
        up_threshold = (n - opens) * (opens - closes + 2)
        if rand(rng, 0:(total - 1)) < up_threshold
            steps[i] = Int8(1)
            push!(open_stack, Int32(i - 1))
            opens += 1
        else
            steps[i] = Int8(-1)
            isempty(open_stack) && throw(ArgumentError("uniform Dyck sampler encountered an invalid close step"))
            j = pop!(open_stack)
            adj[i] = j
            adj[Int(j) + 1] = Int32(i - 1)
            closes += 1
        end
    end

    isempty(open_stack) || throw(ArgumentError("uniform Dyck sampler left unmatched openings"))
    return steps, adj
end

function sample_half_plane_bracket_path_with_adjacency(order::Integer, boundary_half_length::Integer, rng::AbstractRNG=Random.default_rng())
    n = Int(order)
    h = Int(boundary_half_length)
    n >= 1 || throw(ArgumentError("order must be >= 1"))
    0 <= h <= n || throw(ArgumentError("boundary_half_length must lie in [0, order]"))

    steps = Vector{Int8}(undef, 2n)
    adj = fill(Int32(-1), 2n)
    open_stack = Int32[]
    sizehint!(open_stack, n)

    opens = 0
    closes = 0
    prod1 = BigInt(n + h + 1)
    prod2 = BigInt(n - h)

    for i in 1:(2n)
        total = BigInt(2n - opens - closes) * (prod1 - prod2)
        up_threshold = BigInt(n + h - opens) * prod1 - BigInt(n - h - opens - 1) * prod2
        total > 0 || throw(ArgumentError("failed to count half-plane bracket completions"))

        if rand(rng, big(0):(total - 1)) < up_threshold
            steps[i] = Int8(1)
            push!(open_stack, Int32(i - 1))
            opens += 1
            prod1 *= BigInt(n + h - opens + 1)
            prod2 *= BigInt(n - h - opens)
        else
            steps[i] = Int8(-1)
            isempty(open_stack) && throw(ArgumentError("half-plane bracket sampler encountered an invalid close step"))
            j = pop!(open_stack)
            adj[i] = j
            adj[Int(j) + 1] = Int32(i - 1)
            prod1 = div(prod1, BigInt(n + h - closes + 1))
            denom2 = n - h - closes
            prod2 = denom2 == 0 ? big(0) : div(prod2, BigInt(denom2))
            closes += 1
        end
    end

    return steps, adj, copy(open_stack)
end

function sample_half_plane_bracket_path(order::Integer, boundary_half_length::Integer, rng::AbstractRNG=Random.default_rng())
    steps, _, _ = sample_half_plane_bracket_path_with_adjacency(order, boundary_half_length, rng)
    return steps
end

function _pair_bracket_steps(steps::AbstractVector{<:Integer}; allow_unmatched::Bool=false)
    seq = Int8.(collect(steps))
    adj = fill(Int32(-1), length(seq))
    stack = Int[]

    for i in eachindex(seq)
        if seq[i] == 1
            push!(stack, i - 1)
        elseif seq[i] == -1
            isempty(stack) && throw(ArgumentError("invalid bracket path: too many down-steps"))
            j = pop!(stack)
            adj[i] = Int32(j)
            adj[j + 1] = Int32(i - 1)
        else
            throw(ArgumentError("bracket steps must be +/-1"))
        end
    end

    if allow_unmatched
        return adj, Int32.(stack)
    end

    isempty(stack) || throw(ArgumentError("invalid balanced bracket path: unmatched up-steps remain"))
    return adj
end

function _matching_edges(adj::AbstractVector{<:Integer}, group::AbstractString)
    u = Int32[]
    v = Int32[]
    g = String[]
    for i in eachindex(adj)
        j = Int(adj[i])
        j > i - 1 || continue
        push!(u, Int32(i - 1))
        push!(v, Int32(j))
        push!(g, String(group))
    end
    return u, v, g
end

function _cycle_edge_vectors(num_vertices::Integer)
    n = Int(num_vertices)
    u = Vector{Int32}(undef, n)
    v = Vector{Int32}(undef, n)
    for i in 0:(n - 1)
        u[i + 1] = Int32(i)
        v[i + 1] = Int32(mod(i + 1, n))
    end
    return u, v
end

function _path_edge_vectors(num_vertices::Integer)
    n = Int(num_vertices)
    n <= 1 && return Int32[], Int32[]
    u = Vector{Int32}(undef, n - 1)
    v = Vector{Int32}(undef, n - 1)
    for i in 0:(n - 2)
        u[i + 1] = Int32(i)
        v[i + 1] = Int32(i + 1)
    end
    return u, v
end

function _union_find_components(num_vertices::Integer, edges)
    n = Int(num_vertices)
    parent = collect(1:n)
    rank = zeros(Int, n)

    function find_root(x::Int)
        y = x
        while parent[y] != y
            parent[y] = parent[parent[y]]
            y = parent[y]
        end
        return y
    end

    function union!(a::Int, b::Int)
        ra = find_root(a)
        rb = find_root(b)
        ra == rb && return
        if rank[ra] < rank[rb]
            parent[ra] = rb
        elseif rank[ra] > rank[rb]
            parent[rb] = ra
        else
            parent[rb] = ra
            rank[ra] += 1
        end
        return
    end

    arr = sanitize_edge_array(edges)
    for i in 1:size(arr, 1)
        union!(Int(arr[i, 1]) + 1, Int(arr[i, 2]) + 1)
    end

    labels = Vector{Int32}(undef, n)
    counts = Dict{Int,Int}()
    for i in 1:n
        root = find_root(i)
        labels[i] = Int32(root - 1)
        counts[root] = get(counts, root, 0) + 1
    end

    sizes = sort!(Int32.(collect(values(counts))); rev=true)
    return labels, sizes
end

function _meandric_component_labels_and_sizes(
    upper_adj::AbstractVector{<:Integer},
    lower_adj::AbstractVector{<:Integer};
    extra_edges=nothing,
)
    n = length(upper_adj)
    length(lower_adj) == n || throw(ArgumentError("upper and lower matchings must have the same number of vertices"))

    extra_partner = fill(Int32(-1), n)
    if extra_edges !== nothing
        arr = sanitize_edge_array(extra_edges)
        for i in 1:size(arr, 1)
            u = Int(arr[i, 1]) + 1
            v = Int(arr[i, 2]) + 1
            extra_partner[u] = Int32(v - 1)
            extra_partner[v] = Int32(u - 1)
        end
    end

    labels = fill(Int32(-1), n)
    counts = Int32[]
    stack = Int32[]
    sizehint!(counts, max(1, fld(n, 2)))
    sizehint!(stack, max(1, fld(n, 2)))

    for start in 0:(n - 1)
        labels[start + 1] >= 0 && continue
        label = Int32(length(counts))
        push!(counts, Int32(0))
        push!(stack, Int32(start))

        while !isempty(stack)
            v = pop!(stack)
            idx = Int(v) + 1
            labels[idx] >= 0 && continue
            labels[idx] = label
            counts[Int(label) + 1] += 1

            upper_nb = Int(upper_adj[idx])
            if upper_nb >= 0 && labels[upper_nb + 1] < 0
                push!(stack, Int32(upper_nb))
            end

            lower_nb = Int(lower_adj[idx])
            if lower_nb >= 0 && lower_nb != upper_nb && labels[lower_nb + 1] < 0
                push!(stack, Int32(lower_nb))
            end

            extra_nb = Int(extra_partner[idx])
            if extra_nb >= 0 && extra_nb != upper_nb && extra_nb != lower_nb && labels[extra_nb + 1] < 0
                push!(stack, Int32(extra_nb))
            end
        end
    end

    return labels, counts
end

function _sorted_meandric_component_labels(
    component_sizes::AbstractVector{<:Integer};
    excluded_labels::AbstractSet{Int32}=Set{Int32}(),
)
    labels = Int32.(collect(0:(length(component_sizes) - 1)))
    isempty(excluded_labels) || filter!(label -> !(label in excluded_labels), labels)
    sort!(labels; by=label -> (-Int(component_sizes[Int(label) + 1]), Int(label)))
    return labels
end

function _sorted_meandric_component_sizes(
    component_sizes::AbstractVector{<:Integer};
    excluded_labels::AbstractSet{Int32}=Set{Int32}(),
)
    labels = _sorted_meandric_component_labels(component_sizes; excluded_labels=excluded_labels)
    return Int32[Int(component_sizes[Int(label) + 1]) for label in labels]
end

function _component_labels_touching_vertices(component_labels::AbstractVector{Int32}, vertices)
    touched = Set{Int32}()
    for v in vertices
        push!(touched, component_labels[Int(v) + 1])
    end
    return touched
end

function _boundary_pair_edges(boundary_vertices, target_index::Union{Nothing,Int}=nothing)
    boundary = Int32.(collect(boundary_vertices))
    p = length(boundary)
    p == 0 && return Matrix{Int32}(undef, 0, 2), nothing
    iseven(p) || throw(ArgumentError("half-plane meandric boundary length must be even"))

    if target_index === nothing
        edges = Matrix{Int32}(undef, fld(p, 2), 2)
        row = 1
        for i in 1:2:p
            edges[row, 1] = boundary[i]
            edges[row, 2] = boundary[i + 1]
            row += 1
        end
        return edges, nothing
    end

    target = Int(target_index)
    1 < target <= p || throw(ArgumentError("target_index must lie in 2:length(boundary_vertices)"))

    rows = fld(p, 2) - 1
    edges = Matrix{Int32}(undef, max(rows, 0), 2)
    row = 1

    for i in 2:2:(target - 2)
        edges[row, 1] = boundary[i]
        edges[row, 2] = boundary[i + 1]
        row += 1
    end
    for i in (target + 1):2:(p - 1)
        edges[row, 1] = boundary[i]
        edges[row, 2] = boundary[i + 1]
        row += 1
    end

    return edges[1:(row - 1), :], (boundary[1], boundary[target])
end

function _select_link_target(boundary_vertices, raw_labels)
    boundary = Int32.(collect(boundary_vertices))
    p = length(boundary)
    p >= 2 || return nothing

    start_label = raw_labels[Int(boundary[1]) + 1]
    candidate = max(2, 2 * round(Int, p / 4))
    candidate = min(candidate, p)
    isodd(candidate) && (candidate -= 1)
    candidate < 2 && (candidate = 2)

    j = candidate
    while j <= p && raw_labels[Int(boundary[j]) + 1] == start_label
        j += 2
    end
    if j > p
        j = 2
        while j <= p && raw_labels[Int(boundary[j]) + 1] == start_label
            j += 2
        end
    end
    return j <= p ? j : nothing
end

function _tree_from_dyck_steps(steps)
    seq = Int8.(collect(steps))
    parent = Int32[-1]
    edge_u = Int32[]
    edge_v = Int32[]
    current = Int32(0)
    contour_vertices = Int32[0]

    for step in seq
        if step == 1
            child = Int32(length(parent))
            push!(parent, current)
            push!(edge_u, current)
            push!(edge_v, child)
            current = child
        elseif step == -1
            current > 0 || throw(ArgumentError("invalid Dyck path for tree contour"))
            current = parent[Int(current) + 1]
        else
            throw(ArgumentError("Dyck path steps must be +/-1"))
        end
        push!(contour_vertices, current)
    end

    current == 0 || throw(ArgumentError("Dyck path contour must return to the root"))
    return parent, edge_u, edge_v, contour_vertices
end

function _build_glued_tree_graph(upper_steps, lower_steps)
    _, upper_u, upper_v, upper_contour = _tree_from_dyck_steps(upper_steps)
    _, lower_u, lower_v, lower_contour = _tree_from_dyck_steps(lower_steps)

    n_upper = isempty(upper_u) ? 1 : Int(maximum(vcat(upper_u, upper_v))) + 1
    n_lower = isempty(lower_u) ? 1 : Int(maximum(vcat(lower_u, lower_v))) + 1
    num_points = length(upper_steps)
    length(lower_steps) == num_points || throw(ArgumentError("upper and lower arc diagrams must have the same size"))

    offset_upper = 0
    offset_lower = n_upper
    offset_interface = n_upper + n_lower

    edge_u = Int32[]
    edge_v = Int32[]
    groups = String[]

    function add_edge!(u::Int32, v::Int32, group::AbstractString)
        push!(edge_u, u)
        push!(edge_v, v)
        push!(groups, String(group))
        return nothing
    end

    for i in eachindex(upper_u)
        add_edge!(Int32(offset_upper) + upper_u[i], Int32(offset_upper) + upper_v[i], "tree_upper")
    end
    for i in eachindex(lower_u)
        add_edge!(Int32(offset_lower) + lower_u[i], Int32(offset_lower) + lower_v[i], "tree_lower")
    end

    render_vertex_indices = Vector{Int32}(undef, num_points)
    for i in 1:num_points
        interface = Int32(offset_interface + i - 1)
        render_vertex_indices[i] = interface
        add_edge!(interface, Int32(offset_upper) + upper_contour[i], "glue_upper")
        add_edge!(interface, Int32(offset_upper) + upper_contour[i + 1], "glue_upper")
        add_edge!(interface, Int32(offset_lower) + lower_contour[i], "glue_lower")
        add_edge!(interface, Int32(offset_lower) + lower_contour[i + 1], "glue_lower")
    end

    return Int32(offset_interface + num_points), edge_u, edge_v, groups, render_vertex_indices
end

function _matching_edges_with_opening_ids(
    adj::AbstractVector{<:Integer},
    group::AbstractString;
    edge_id_offset::Integer=0,
)
    u = Int32[]
    v = Int32[]
    g = String[]
    opening_edge_ids = fill(Int32(-1), length(adj))
    next_edge_id = Int32(edge_id_offset)
    for i in eachindex(adj)
        j = Int(adj[i])
        j > i - 1 || continue
        push!(u, Int32(i - 1))
        push!(v, Int32(j))
        push!(g, String(group))
        opening_edge_ids[i] = next_edge_id
        next_edge_id += 1
    end
    return u, v, g, opening_edge_ids
end

function _matching_children_from_steps(steps::AbstractVector{<:Integer})
    seq = Int8.(collect(steps))
    children = [Int32[] for _ in eachindex(seq)]
    open_stack = Int32[]
    top_level_openings = Int32[]
    sizehint!(open_stack, fld(length(seq), 2))

    for i in eachindex(seq)
        step = Int(seq[i])
        if step == 1
            isempty(open_stack) && push!(top_level_openings, Int32(i - 1))
            push!(open_stack, Int32(i - 1))
        elseif step == -1
            isempty(open_stack) && throw(ArgumentError("invalid balanced bracket path while building meandric faces"))
            opening = pop!(open_stack)
            if !isempty(open_stack)
                push!(children[Int(open_stack[end]) + 1], opening)
            end
        else
            throw(ArgumentError("bracket steps must be +/-1"))
        end
    end

    isempty(open_stack) || throw(ArgumentError("invalid balanced bracket path while building meandric faces"))
    return children, top_level_openings
end

function _matching_face_cycles(
    steps::AbstractVector{<:Integer},
    adj::AbstractVector{<:Integer},
    opening_edge_ids::AbstractVector{<:Integer},
    root_closure_edge_id::Integer,
)
    n = length(adj)
    length(steps) == n || throw(ArgumentError("steps and adjacency must have the same length"))
    length(opening_edge_ids) == n || throw(ArgumentError("opening_edge_ids must align with adjacency"))

    children, top_level_openings = _matching_children_from_steps(steps)
    faces = Vector{Vector{Int32}}()
    face_edge_ids = Vector{Vector{Int32}}()
    sizehint!(faces, fld(n, 2) + 1)
    sizehint!(face_edge_ids, fld(n, 2) + 1)

    for opening in 0:(n - 1)
        closing = Int(adj[opening + 1])
        closing > opening || continue

        vertices = Int32[Int32(opening)]
        edge_cycle = Int32[]
        current = opening

        for child_opening32 in children[opening + 1]
            child_opening = Int(child_opening32)
            child_opening == current + 1 ||
                throw(ArgumentError("meandric face construction expected alternating generic/loop edges"))

            push!(edge_cycle, Int32(current))
            current += 1
            push!(vertices, Int32(current))

            child_edge_id = Int32(opening_edge_ids[child_opening + 1])
            child_edge_id >= 0 || throw(ArgumentError("missing child arc edge id in meandric face construction"))
            push!(edge_cycle, child_edge_id)

            current = Int(adj[child_opening + 1])
            push!(vertices, Int32(current))
        end

        closing == current + 1 ||
            throw(ArgumentError("meandric face construction expected alternating generic/loop edges"))

        push!(edge_cycle, Int32(current))
        current += 1
        push!(vertices, Int32(current))

        outer_edge_id = Int32(opening_edge_ids[opening + 1])
        outer_edge_id >= 0 || throw(ArgumentError("missing outer arc edge id in meandric face construction"))
        push!(edge_cycle, outer_edge_id)

        length(vertices) == length(edge_cycle) ||
            throw(ArgumentError("meandric face construction produced inconsistent face data"))

        push!(faces, vertices)
        push!(face_edge_ids, edge_cycle)
    end

    isempty(top_level_openings) && return faces, face_edge_ids

    root_vertices = Int32[]
    root_edge_cycle = Int32[]
    current = Int(top_level_openings[1])
    current == 0 || throw(ArgumentError("balanced meandric matching must start at the leftmost vertex"))
    push!(root_vertices, Int32(current))

    for (idx, opening32) in enumerate(top_level_openings)
        opening = Int(opening32)
        idx == 1 || begin
            opening == current + 1 ||
                throw(ArgumentError("meandric root face construction expected alternating generic/loop edges"))
            push!(root_edge_cycle, Int32(current))
            current += 1
            push!(root_vertices, Int32(current))
        end

        opening == current || throw(ArgumentError("meandric root face construction expected consecutive top-level arcs"))
        arc_edge_id = Int32(opening_edge_ids[opening + 1])
        arc_edge_id >= 0 || throw(ArgumentError("missing top-level arc edge id in meandric root face construction"))
        push!(root_edge_cycle, arc_edge_id)

        current = Int(adj[opening + 1])
        push!(root_vertices, Int32(current))
    end

    current == n - 1 || throw(ArgumentError("balanced meandric matching must end at the rightmost vertex"))
    push!(root_edge_cycle, Int32(root_closure_edge_id))

    length(root_vertices) == length(root_edge_cycle) ||
        throw(ArgumentError("meandric root face construction produced inconsistent face data"))

    push!(faces, root_vertices)
    push!(face_edge_ids, root_edge_cycle)

    return faces, face_edge_ids
end

function _closed_meandric_faces(
    upper_steps::AbstractVector{<:Integer},
    upper_adj::AbstractVector{<:Integer},
    upper_opening_edge_ids::AbstractVector{<:Integer},
    lower_steps::AbstractVector{<:Integer},
    lower_adj::AbstractVector{<:Integer},
    lower_opening_edge_ids::AbstractVector{<:Integer},
)
    root_closure_edge_id = length(upper_adj) - 1 + fld(length(upper_adj), 2) + fld(length(lower_adj), 2)
    upper_faces, upper_face_edge_ids = _matching_face_cycles(upper_steps, upper_adj, upper_opening_edge_ids, root_closure_edge_id)
    lower_faces, lower_face_edge_ids = _matching_face_cycles(lower_steps, lower_adj, lower_opening_edge_ids, root_closure_edge_id)
    face_edge_group = vcat(
        fill("generic", length(upper_adj) - 1),
        fill("upper", fld(length(upper_adj), 2)),
        fill("lower", fld(length(lower_adj), 2)),
        ["generic"],
    )
    return vcat(upper_faces, lower_faces), vcat(upper_face_edge_ids, lower_face_edge_ids), face_edge_group
end

function _closed_meandric_component_sizes(upper_adj, lower_adj)
    _, component_sizes = _meandric_component_labels_and_sizes(upper_adj, lower_adj)
    return _sorted_meandric_component_sizes(component_sizes)
end

function _build_closed_meandric_map(
    variant::Symbol,
    upper_steps,
    lower_steps;
    sampler_metadata=Dict{String,Any}(),
)
    variant in (:uniform_meandric, :uniform_meander) ||
        throw(ArgumentError("unsupported meandric variant $(repr(variant))"))

    upper_seq = Int8.(collect(upper_steps))
    lower_seq = Int8.(collect(lower_steps))
    length(upper_seq) == length(lower_seq) || throw(ArgumentError("upper and lower step sequences must have the same length"))
    !isempty(upper_seq) || throw(ArgumentError("meandric systems require non-empty step sequences"))
    iseven(length(upper_seq)) || throw(ArgumentError("meandric step sequences must have even length"))

    order = length(upper_seq) ÷ 2
    upper_adj = _pair_bracket_steps(upper_seq; allow_unmatched=false)
    lower_adj = _pair_bracket_steps(lower_seq; allow_unmatched=false)

    strand_u, strand_v = _path_edge_vectors(length(upper_seq))
    upper_u, upper_v, upper_g, upper_opening_edge_ids = _matching_edges_with_opening_ids(
        upper_adj,
        "upper";
        edge_id_offset=length(strand_u),
    )
    lower_u, lower_v, lower_g, lower_opening_edge_ids = _matching_edges_with_opening_ids(
        lower_adj,
        "lower";
        edge_id_offset=length(strand_u) + length(upper_u),
    )
    render_u = vcat(strand_u, upper_u, lower_u)
    render_v = vcat(strand_v, upper_v, lower_v)
    render_groups = vcat(fill("generic", length(strand_u)), upper_g, lower_g)
    faces, face_edge_ids, face_edge_group = _closed_meandric_faces(
        upper_seq,
        upper_adj,
        upper_opening_edge_ids,
        lower_seq,
        lower_adj,
        lower_opening_edge_ids,
    )

    component_sizes = _closed_meandric_component_sizes(upper_adj, lower_adj)
    aux_n, aux_u, aux_v, aux_groups, render_vertex_indices = _build_glued_tree_graph(upper_seq, lower_seq)

    sampler_info = Dict{String,Any}(string(k) => v for (k, v) in pairs(sampler_metadata))

    return UniformMeandricSystemMap(
        variant,
        Int32(order),
        upper_seq,
        lower_seq,
        upper_adj,
        lower_adj,
        render_u,
        render_v,
        render_groups,
        faces,
        face_edge_ids,
        face_edge_group,
        component_sizes,
        aux_n,
        aux_u,
        aux_v,
        aux_groups,
        render_vertex_indices,
        sampler_info,
    )
end

function _find_peak_swap_parent_opening(adj::AbstractVector{<:Integer}, opening_pos::Int32)
    probe = Int(opening_pos) - 1
    leaf_close = Int(opening_pos) + 1

    @inbounds while probe >= 0
        mate = Int(adj[probe + 1])
        if mate > probe
            mate > leaf_close && return Int32(probe)
            probe -= 1
        else
            probe = mate - 1
        end
    end

    return Int32(-1)
end

function _build_dyck_swap_index_pool(steps::AbstractVector{<:Integer})
    num_boundaries = max(length(steps) - 1, 0)
    indices = Int32[]
    slot_of = fill(Int32(0), num_boundaries)

    for idx in 1:num_boundaries
        if steps[idx] != steps[idx + 1]
            push!(indices, Int32(idx))
            slot_of[idx] = Int32(length(indices))
        end
    end

    return _DyckSwapIndexPool(indices, slot_of)
end

function _refresh_dyck_swap_index!(pool::_DyckSwapIndexPool, steps::AbstractVector{<:Integer}, idx::Integer)
    1 <= idx < length(steps) || return nothing

    should_be_active = steps[idx] != steps[idx + 1]
    slot = Int(pool.slot_of[idx])

    if should_be_active
        if slot == 0
            push!(pool.indices, Int32(idx))
            pool.slot_of[idx] = Int32(length(pool.indices))
        end
    elseif slot != 0
        moved = Int(pop!(pool.indices))
        if slot <= length(pool.indices)
            pool.indices[slot] = Int32(moved)
            pool.slot_of[moved] = Int32(slot)
        end
        pool.slot_of[idx] = Int32(0)
    end

    return nothing
end

function _refresh_dyck_swap_index_window!(pool::_DyckSwapIndexPool, steps::AbstractVector{<:Integer}, idx::Integer)
    _refresh_dyck_swap_index!(pool, steps, idx - 1)
    _refresh_dyck_swap_index!(pool, steps, idx)
    _refresh_dyck_swap_index!(pool, steps, idx + 1)
    return nothing
end

function _dyck_adjacent_swap_effect(steps::AbstractVector{<:Integer}, adj::AbstractVector{<:Integer}, idx::Integer)
    n = length(steps)
    1 <= idx < n || return nothing

    left = Int(steps[idx])
    right = Int(steps[idx + 1])
    left == right && return nothing

    if left == -1 && right == 1
        close_pos = Int32(idx - 1)
        open_pos = Int32(idx)
        outer_open = Int32(adj[idx])
        outer_close = Int32(adj[idx + 1])
        outer_open >= 0 || return nothing
        outer_close >= 0 || return nothing
        return (
            test_u=outer_open,
            test_v=open_pos,
            new1_u=outer_open,
            new1_v=outer_close,
            new2_u=close_pos,
            new2_v=open_pos,
        )
    elseif left == 1 && right == -1
        open_pos = Int32(idx - 1)
        close_pos = Int32(idx)
        Int32(adj[idx]) == close_pos || return nothing
        parent_open = _find_peak_swap_parent_opening(adj, open_pos)
        parent_open >= 0 || return nothing
        parent_close = Int32(adj[Int(parent_open) + 1])
        return (
            test_u=parent_open,
            test_v=open_pos,
            new1_u=parent_open,
            new1_v=open_pos,
            new2_u=close_pos,
            new2_v=parent_close,
        )
    end

    return nothing
end

function _meandric_vertices_share_component(
    upper_adj::AbstractVector{<:Integer},
    lower_adj::AbstractVector{<:Integer},
    start_vertex::Int32,
    target_vertex::Int32,
)
    start_vertex == target_vertex && return true
    curr = start_vertex

    @inbounds while true
        curr = Int32(upper_adj[Int(curr) + 1])
        curr == target_vertex && return true
        curr = Int32(lower_adj[Int(curr) + 1])
        curr == target_vertex && return true
        curr == start_vertex && return false
    end
end

@inline function _set_matching_pair!(adj::Vector{Int32}, u::Int32, v::Int32)
    @inbounds begin
        adj[Int(u) + 1] = v
        adj[Int(v) + 1] = u
    end
    return nothing
end

function _apply_dyck_adjacent_swap!(
    steps::Vector{Int8},
    adj::Vector{Int32},
    pool::_DyckSwapIndexPool,
    idx::Integer,
    effect,
)
    @inbounds steps[idx], steps[idx + 1] = steps[idx + 1], steps[idx]
    _set_matching_pair!(adj, effect.new1_u, effect.new1_v)
    _set_matching_pair!(adj, effect.new2_u, effect.new2_v)
    _refresh_dyck_swap_index_window!(pool, steps, idx)
    return nothing
end

function _accept_loop_penalty_delta(delta::Integer, q::Real, rng::AbstractRNG)
    delta <= 0 && return true

    qf = float(q)
    qf <= 0.0 && return false
    return rand(rng) < qf^delta
end

function _meandric_mcmc_step(
    upper_steps::Vector{Int8},
    lower_steps::Vector{Int8},
    upper_adj::Vector{Int32},
    lower_adj::Vector{Int32},
    upper_pool::_DyckSwapIndexPool,
    lower_pool::_DyckSwapIndexPool,
    current_loops::Int,
    q::Real,
    rng::AbstractRNG,
)
    update_upper = rand(rng, Bool)
    steps = update_upper ? upper_steps : lower_steps
    adj = update_upper ? upper_adj : lower_adj
    pool = update_upper ? upper_pool : lower_pool
    isempty(pool.indices) && return upper_steps, lower_steps, upper_adj, lower_adj, current_loops, false, false
    idx = Int(rand(rng, pool.indices))
    effect = _dyck_adjacent_swap_effect(steps, adj, idx)

    if effect === nothing
        return upper_steps, lower_steps, upper_adj, lower_adj, current_loops, false, false
    end

    same_component = _meandric_vertices_share_component(upper_adj, lower_adj, effect.test_u, effect.test_v)
    delta_loops = same_component ? 1 : -1

    if _accept_loop_penalty_delta(delta_loops, q, rng)
        _apply_dyck_adjacent_swap!(steps, adj, pool, idx, effect)
        new_loop_count = current_loops + delta_loops
        return upper_steps, lower_steps, upper_adj, lower_adj, new_loop_count, true, true
    end

    return upper_steps, lower_steps, upper_adj, lower_adj, current_loops, true, false
end

function _normalize_tempering_schedule(schedule)
    values = if schedule === nothing
        Float64[1.0, 0.5, 0.25, 0.125, 0.0625]
    else
        [float(q) for q in schedule]
    end

    isempty(values) && throw(ArgumentError("temper_schedule may not be empty"))
    for q in values
        0.0 <= q <= 1.0 || throw(ArgumentError("temper_schedule values must lie in [0, 1]"))
    end

    values = sort!(unique(values); rev=true)
    values[1] < 1.0 && pushfirst!(values, 1.0)
    positive = [q for q in values if q > 0.0]
    isempty(positive) && (positive = Float64[0.0625])
    return values, positive[end]
end

function _sample_uniform_meander_steps(
    order::Integer,
    rng::AbstractRNG;
    temper_schedule=nothing,
    sweeps_per_temperature=nothing,
    search_sweeps=nothing,
    descent_sweeps=nothing,
    mixing_sweeps=nothing,
    restarts::Integer=8,
)
    n = Int(order)
    n >= 1 || throw(ArgumentError("order must be >= 1"))

    schedule, search_q = _normalize_tempering_schedule(temper_schedule)
    anneal_schedule = [q for q in schedule if q > 0.0]
    isempty(anneal_schedule) && (anneal_schedule = [search_q])

    sweeps_temp = sweeps_per_temperature === nothing ? max(8, 2 * ceil(Int, sqrt(n))) : Int(sweeps_per_temperature)
    search_sweeps_int = search_sweeps === nothing ? max(24, 4 * sweeps_temp) : Int(search_sweeps)
    descent_sweeps_int = descent_sweeps === nothing ? max(8, sweeps_temp) : Int(descent_sweeps)
    mixing_sweeps_int = mixing_sweeps === nothing ? max(12, 2 * sweeps_temp) : Int(mixing_sweeps)
    restart_count = Int(restarts)

    sweeps_temp >= 1 || throw(ArgumentError("sweeps_per_temperature must be >= 1"))
    search_sweeps_int >= 1 || throw(ArgumentError("search_sweeps must be >= 1"))
    descent_sweeps_int >= 0 || throw(ArgumentError("descent_sweeps must be >= 0"))
    mixing_sweeps_int >= 0 || throw(ArgumentError("mixing_sweeps must be >= 0"))
    restart_count >= 1 || throw(ArgumentError("restarts must be >= 1"))

    steps_per_sweep = max(1, 2n)
    total_proposals = 0
    valid_proposals = 0
    accepted_proposals = 0
    best_loops_seen = typemax(Int)
    attempts_used = restart_count

    for attempt in 1:restart_count
        upper_steps, upper_adj = sample_uniform_dyck_path_with_adjacency(n, rng)
        lower_steps, lower_adj = sample_uniform_dyck_path_with_adjacency(n, rng)
        upper_pool = _build_dyck_swap_index_pool(upper_steps)
        lower_pool = _build_dyck_swap_index_pool(lower_steps)
        current_loops = length(_closed_meandric_component_sizes(upper_adj, lower_adj))
        best_loops_seen = min(best_loops_seen, current_loops)

        for q in anneal_schedule
            for _ in 1:(sweeps_temp * steps_per_sweep)
                upper_steps, lower_steps, upper_adj, lower_adj, current_loops, valid_move, accepted_move =
                    _meandric_mcmc_step(
                        upper_steps,
                        lower_steps,
                        upper_adj,
                        lower_adj,
                        upper_pool,
                        lower_pool,
                        current_loops,
                        q,
                        rng,
                    )
                total_proposals += 1
                valid_move && (valid_proposals += 1)
                accepted_move && (accepted_proposals += 1)
                best_loops_seen = min(best_loops_seen, current_loops)
            end
        end

        for _ in 1:(search_sweeps_int * steps_per_sweep)
            current_loops == 1 && break
            upper_steps, lower_steps, upper_adj, lower_adj, current_loops, valid_move, accepted_move =
                _meandric_mcmc_step(
                    upper_steps,
                    lower_steps,
                    upper_adj,
                    lower_adj,
                    upper_pool,
                    lower_pool,
                    current_loops,
                    search_q,
                    rng,
                )
            total_proposals += 1
            valid_move && (valid_proposals += 1)
            accepted_move && (accepted_proposals += 1)
            best_loops_seen = min(best_loops_seen, current_loops)
        end

        for _ in 1:(descent_sweeps_int * steps_per_sweep)
            current_loops == 1 && break
            upper_steps, lower_steps, upper_adj, lower_adj, current_loops, valid_move, accepted_move =
                _meandric_mcmc_step(
                    upper_steps,
                    lower_steps,
                    upper_adj,
                    lower_adj,
                    upper_pool,
                    lower_pool,
                    current_loops,
                    0.0,
                    rng,
                )
            total_proposals += 1
            valid_move && (valid_proposals += 1)
            accepted_move && (accepted_proposals += 1)
            best_loops_seen = min(best_loops_seen, current_loops)
        end

        if current_loops == 1
            for _ in 1:(mixing_sweeps_int * steps_per_sweep)
                upper_steps, lower_steps, upper_adj, lower_adj, current_loops, valid_move, accepted_move =
                    _meandric_mcmc_step(
                        upper_steps,
                        lower_steps,
                        upper_adj,
                        lower_adj,
                        upper_pool,
                        lower_pool,
                        current_loops,
                        0.0,
                        rng,
                    )
                total_proposals += 1
                valid_move && (valid_proposals += 1)
                accepted_move && (accepted_proposals += 1)
            end

            current_loops == 1 || throw(ErrorException("zero-temperature meander mixing escaped the single-loop state"))

            attempts_used = attempt
            metadata = Dict{String,Any}(
                "sampler" => "tempered_mcmc",
                "temper_schedule" => schedule,
                "search_q" => search_q,
                "sweeps_per_temperature" => sweeps_temp,
                "search_sweeps" => search_sweeps_int,
                "descent_sweeps" => descent_sweeps_int,
                "mixing_sweeps" => mixing_sweeps_int,
                "restarts" => restart_count,
                "attempts_used" => attempts_used,
                "proposals" => total_proposals,
                "valid_proposals" => valid_proposals,
                "accepted_proposals" => accepted_proposals,
                "acceptance_rate" => total_proposals == 0 ? 0.0 : accepted_proposals / total_proposals,
                "valid_acceptance_rate" => valid_proposals == 0 ? 0.0 : accepted_proposals / valid_proposals,
                "final_loops" => 1,
                "best_loops_seen" => best_loops_seen,
            )
            return upper_steps, lower_steps, metadata
        end
    end

    throw(ErrorException(
        "failed to sample a single-loop meander after $(restart_count) tempered-MCMC restart(s); " *
        "best_loops_seen=$(best_loops_seen). Try increasing `restarts`, `search_sweeps`, `descent_sweeps`, or `sweeps_per_temperature`.",
    ))
end

function build_half_plane_meandric_system(; order::Integer, boundary_half_length=nothing, seed::Integer=1)
    n = Int(order)
    n >= 1 || throw(ArgumentError("order must be >= 1"))
    bh = boundary_half_length === nothing ? max(1, round(Int, sqrt(float(n)))) : Int(boundary_half_length)
    1 <= bh <= n || throw(ArgumentError("boundary_half_length must lie in [1, order]"))

    rng = MersenneTwister(Int(seed))
    upper_steps, upper_adj = sample_uniform_dyck_path_with_adjacency(n, rng)
    lower_steps, lower_adj, boundary = sample_half_plane_bracket_path_with_adjacency(n, bh, rng)

    upper_u, upper_v, upper_g = _matching_edges(upper_adj, "upper")
    lower_u, lower_v, lower_g = _matching_edges(lower_adj, "lower")
    contour_u, contour_v = _cycle_edge_vectors(2n)
    strand_u, strand_v = _path_edge_vectors(2n)

    render_u = vcat(strand_u, upper_u, lower_u)
    render_v = vcat(strand_v, upper_v, lower_v)
    render_groups = vcat(fill("generic", length(strand_u)), upper_g, lower_g)

    layout_u = vcat(contour_u, upper_u, lower_u)
    layout_v = vcat(contour_v, upper_v, lower_v)

    raw_labels, raw_counts = _meandric_component_labels_and_sizes(upper_adj, lower_adj)
    raw_sizes = _sorted_meandric_component_sizes(raw_counts)
    target_index = _select_link_target(boundary, raw_labels)
    linked_edges, linked_interface = target_index === nothing ?
        (Matrix{Int32}(undef, 0, 2), nothing) :
        _boundary_pair_edges(boundary, target_index)
    cle_edges, _ = _boundary_pair_edges(boundary, nothing)

    _, linked_counts = _meandric_component_labels_and_sizes(upper_adj, lower_adj; extra_edges=linked_edges)
    _, cle_counts = _meandric_component_labels_and_sizes(upper_adj, lower_adj; extra_edges=cle_edges)
    linked_sizes = _sorted_meandric_component_sizes(linked_counts)
    cle_sizes = _sorted_meandric_component_sizes(cle_counts)

    return HalfPlaneMeandricSystemMap(
        Int32(n),
        Int32(bh),
        upper_steps,
        lower_steps,
        upper_adj,
        lower_adj,
        boundary,
        render_u,
        render_v,
        render_groups,
        layout_u,
        layout_v,
        raw_sizes,
        linked_sizes,
        cle_sizes,
        linked_interface,
    )
end

generate_half_plane_meandric_system(; kwargs...) = build_half_plane_meandric_system(; kwargs...)

function build_uniform_meandric_system(; order::Integer, seed::Integer=1)
    n = Int(order)
    n >= 1 || throw(ArgumentError("order must be >= 1"))
    rng = MersenneTwister(Int(seed))

    upper_steps, _ = sample_uniform_dyck_path_with_adjacency(n, rng)
    lower_steps, _ = sample_uniform_dyck_path_with_adjacency(n, rng)

    return _build_closed_meandric_map(
        :uniform_meandric,
        upper_steps,
        lower_steps;
        sampler_metadata=Dict{String,Any}("sampler" => "independent_uniform_dyck_paths"),
    )
end

generate_uniform_meandric_system(; kwargs...) = build_uniform_meandric_system(; kwargs...)

function build_uniform_meander(
    ;
    order::Integer,
    seed::Integer=1,
    temper_schedule=nothing,
    sweeps_per_temperature=nothing,
    search_sweeps=nothing,
    descent_sweeps=nothing,
    mixing_sweeps=nothing,
    restarts::Integer=8,
)
    rng = MersenneTwister(Int(seed))
    upper_steps, lower_steps, sampler_metadata = _sample_uniform_meander_steps(
        order,
        rng;
        temper_schedule=temper_schedule,
        sweeps_per_temperature=sweeps_per_temperature,
        search_sweeps=search_sweeps,
        descent_sweeps=descent_sweeps,
        mixing_sweeps=mixing_sweeps,
        restarts=restarts,
    )

    map_data = _build_closed_meandric_map(
        :uniform_meander,
        upper_steps,
        lower_steps;
        sampler_metadata=sampler_metadata,
    )
    length(map_data.component_sizes) == 1 || throw(ErrorException("uniform meander sampler returned a multi-loop meandric system"))
    return map_data
end

generate_uniform_meander(; kwargs...) = build_uniform_meander(; kwargs...)

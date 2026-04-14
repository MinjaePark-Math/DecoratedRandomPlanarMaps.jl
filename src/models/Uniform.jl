# Uniform quadrangulation model and Schaeffer-style construction.

const ROOT_COLOR = Int8(0)
const EDGE_COLOR = Int8(1)

const UNIFORM_COLOR_NAME = Dict(
    Int(ROOT_COLOR) => "root",
    Int(EDGE_COLOR) => "edge",
)

struct UniformMap <: AbstractRandomPlanarMap
    dyck_path::Vector{Int8}
    matching::Vector{Int32}          # values are 0-based step indices
    edge_labels::Vector{Int8}
    corner_labels::Vector{Int32}
    vertex_labels::Vector{Int32}
    faces::Vector{Vector{Int32}}     # vertex ids are 0-based
    face_edge_ids::Vector{Vector{Int32}}
    edge_u::Vector{Int32}
    edge_v::Vector{Int32}
    edge_color::Vector{Int8}
    root_edge::Int32                 # 0-based undirected edge id
    construction_step::Vector{Int32}
    construction_match::Vector{Int32}
    construction_edge_label::Vector{Int8}
    construction_cur_label::Vector{Int32}
    construction_next_label::Vector{Int32}
    construction_stack_height::Vector{Int32}
    construction_case::Vector{String}
end

const UniformQuadrangulationMap = UniformMap
const UniformQuadrangulation = UniformMap

num_vertices(m::UniformMap) = length(m.vertex_labels)
num_faces(m::UniformMap) = length(m.faces)
num_edges(m::UniformMap) = length(m.edge_u)
min_label(m::UniformMap) = minimum(m.corner_labels)
max_label(m::UniformMap) = maximum(m.corner_labels)

function edge_pairs(
    m::UniformMap;
    colors=nothing,
    collapse::Bool=false,
    drop_loops::Bool=false,
)::Matrix{Int32}
    mask = trues(length(m.edge_u))
    if colors !== nothing
        vals = Int8[Int(c) for c in colors]
        if isempty(vals)
            return Matrix{Int32}(undef, 0, 2)
        elseif length(vals) == 1
            mask .&= (m.edge_color .== vals[1])
        else
            allowed = Set(vals)
            mask .&= [c in allowed for c in m.edge_color]
        end
    end

    pairs = Matrix{Int32}(undef, count(mask), 2)
    k = 0
    for i in eachindex(mask)
        mask[i] || continue
        u = m.edge_u[i]
        v = m.edge_v[i]
        if drop_loops && u == v
            continue
        end
        k += 1
        pairs[k, 1] = u
        pairs[k, 2] = v
    end
    pairs = pairs[1:k, :]
    if collapse
        return first(collapse_undirected_edges(pairs; drop_loops=drop_loops))
    end
    return pairs
end

collapsed_layout_edges(m::UniformMap; drop_loops::Bool=true) = edge_pairs(m; collapse=true, drop_loops=drop_loops)
layout_edges(m::UniformMap; drop_loops::Bool=true) = collapsed_layout_edges(m; drop_loops=drop_loops)

function matching_from_dyck_path(path)::Vector{Int32}
    steps = Int8.(collect(path))
    stack = Int[]
    matching = fill(Int32(-1), length(steps))
    for i in eachindex(steps)
        step = Int(steps[i])
        if step == 1
            push!(stack, i - 1)
        elseif step == -1
            isempty(stack) && throw(ArgumentError("invalid Dyck path: too many down-steps"))
            j = pop!(stack)
            matching[i] = Int32(j)
            matching[j + 1] = Int32(i - 1)
        else
            throw(ArgumentError("Dyck path steps must be ±1"))
        end
    end
    isempty(stack) || throw(ArgumentError("invalid Dyck path: unmatched up-steps"))
    return matching
end

function sample_schaeffer_data(n::Integer; seed=nothing)
    n_int = Int(n)
    n_int >= 1 || throw(ArgumentError("n must be >= 1"))
    rng = seed === nothing ? Random.default_rng() : MersenneTwister(Int(seed))
    dyck = sample_uniform_dyck_path(n_int, rng)
    matching = matching_from_dyck_path(dyck)

    edge_labels = Vector{Int8}(undef, 2n_int)
    for i in 0:(2n_int - 1)
        j = Int(matching[i + 1])
        if j > i
            edge_labels[i + 1] = Int8(rand(rng, -1:1))
        else
            edge_labels[i + 1] = Int8(-Int(edge_labels[j + 1]))
        end
    end
    return dyck, matching, edge_labels
end

function extract_undirected_edges(m::PlanarMap, root_halfedge::Integer)
    edge_u = Int32[]
    edge_v = Int32[]
    edge_color = Int8[]
    halfedge_to_undirected = fill(Int32(-1), length(m.edges))
    root_pair = sort((Int(root_halfedge), Int(m.edges[Int(root_halfedge)].adj)))
    root_edge_id = -1

    for eid in eachindex(m.edges)
        adj = Int(m.edges[eid].adj)
        eid > adj && continue

        if m.edges[eid].vertex < 0 || m.edges[adj].vertex < 0
            throw(ArgumentError("vertex ids must be assigned before extracting edges"))
        end

        push!(edge_u, m.edges[eid].vertex)
        push!(edge_v, m.edges[adj].vertex)
        undirected_id = Int32(length(edge_u) - 1)
        halfedge_to_undirected[eid] = undirected_id
        halfedge_to_undirected[adj] = undirected_id

        pair = sort((eid, adj))
        color = pair == root_pair ? ROOT_COLOR : EDGE_COLOR
        push!(edge_color, color)
        if color == ROOT_COLOR
            root_edge_id = length(edge_u) - 1
        end
    end

    root_edge_id >= 0 || error("failed to identify the root edge")
    return edge_u, edge_v, edge_color, Int32(root_edge_id), halfedge_to_undirected
end

function extract_faces_with_edge_ids(m::PlanarMap, halfedge_to_undirected)::Tuple{Vector{Vector{Int32}},Vector{Vector{Int32}}}
    seen = Set{Int}()
    faces = Vector{Vector{Int32}}()
    face_edge_ids = Vector{Vector{Int32}}()
    for eid in eachindex(m.edges)
        eid in seen && continue
        cur = eid
        face = Int32[]
        edge_ids = Int32[]
        while true
            push!(seen, cur)
            push!(face, m.edges[cur].vertex)
            undirected_id = halfedge_to_undirected[cur]
            undirected_id >= 0 || throw(ArgumentError("missing undirected edge id for halfedge $cur"))
            push!(edge_ids, undirected_id)
            cur = Int(m.edges[cur].next)
            cur == eid && break
        end
        push!(faces, face)
        push!(face_edge_ids, edge_ids)
    end
    return faces, face_edge_ids
end

function build_uniform_quadrangulation_map(dyck_path, matching, edge_labels)::UniformMap
    dyck = Int8.(collect(dyck_path))
    match = Int32.(collect(matching))
    labels = Int8.(collect(edge_labels))

    ndims(dyck) == 1 || throw(ArgumentError("dyck_path must be a vector"))
    !isempty(dyck) || throw(ArgumentError("dyck_path must be non-empty"))
    iseven(length(dyck)) || throw(ArgumentError("dyck_path length must be even"))
    length(match) == length(dyck) || throw(ArgumentError("matching must have the same length as dyck_path"))
    length(labels) == length(dyck) || throw(ArgumentError("edge_labels must have the same length as dyck_path"))

    n = length(dyck) ÷ 2

    m = PlanarMap()
    root_halfedge = new_double_edge!(m; color=Int(ROOT_COLOR))

    stack = [(0, Int(root_halfedge))]               # (label, halfedge index)
    corner_edges = Int32[root_halfedge]
    corner_labels = Vector{Int32}(undef, 2n)
    corner_labels[1] = Int32(0)

    construction_step = Int32.(collect(0:(2n - 2)))
    construction_match = Vector{Int32}(undef, 2n - 1)
    construction_edge_label = Vector{Int8}(undef, 2n - 1)
    construction_cur_label = Vector{Int32}(undef, 2n - 1)
    construction_next_label = Vector{Int32}(undef, 2n - 1)
    construction_stack_height = Vector{Int32}(undef, 2n - 1)
    construction_case = String[]

    cur_label = 0
    for step in 0:(2n - 2)
        edge_label = Int(labels[step + 1])
        next_label = cur_label + edge_label

        construction_match[step + 1] = match[step + 1]
        construction_edge_label[step + 1] = Int8(edge_label)
        construction_cur_label[step + 1] = Int32(cur_label)
        construction_next_label[step + 1] = Int32(next_label)

        if next_label >= cur_label
            if Int(match[step + 1]) >= step
                if !isempty(stack) && stack[end][1] == next_label
                    e = insert_edge!(m, get_next(m, stack[end][2], 1); color=Int(EDGE_COLOR))
                    stack[end] = (stack[end][1], Int(m.edges[Int(e)].adj))
                    push!(construction_case, "nondecreasing-opening-reuse-level")
                else
                    push!(stack, (next_label, Int(new_double_edge!(m; color=Int(EDGE_COLOR)))))
                    push!(construction_case, "nondecreasing-opening-new-level")
                end
            else
                matched_corner = Int(corner_edges[Int(match[step + 1]) + 1])
                if !isempty(stack) && stack[end][1] == next_label
                    e = insert_edge!(m, matched_corner, get_next(m, stack[end][2], 1); color=Int(EDGE_COLOR))
                    stack[end] = (stack[end][1], Int(e))
                    push!(construction_case, "nondecreasing-closing-reuse-level")
                else
                    e = insert_edge!(m, matched_corner; color=Int(EDGE_COLOR))
                    push!(stack, (next_label, Int(e)))
                    push!(construction_case, "nondecreasing-closing-new-level")
                end
            end
        else
            attach_edge = Int(get_next(m, stack[end][2], 1))
            pop!(stack)

            if !isempty(stack) && stack[end][1] == next_label
                e = insert_edge!(m, attach_edge, get_next(m, stack[end][2], 1); color=Int(EDGE_COLOR))
                stack[end] = (stack[end][1], Int(e))
                push!(construction_case, "descending-reuse-level")
            else
                e = insert_edge!(m, attach_edge; color=Int(EDGE_COLOR))
                push!(stack, (next_label, Int(e)))
                push!(construction_case, "descending-new-level")
            end

            if Int(match[step + 1]) < step
                contract_vertices!(m, stack[end][2], Int(corner_edges[Int(match[step + 1]) + 1]))
                construction_case[end] *= "+contract"
            end
        end

        push!(corner_edges, Int32(stack[end][2]))
        corner_labels[step + 2] = Int32(next_label)
        construction_stack_height[step + 1] = Int32(length(stack))
        cur_label = next_label
    end

    contract_edge = Int(get_previous(m, m.edges[Int(corner_edges[1])].adj, 1 - stack[end][1]))
    while length(stack) > 1
        contract_vertices!(m, get_next(m, stack[end][2], 1), get_next(m, contract_edge, 1))
        contract_edge = Int(get_previous(m, contract_edge, 1))
        pop!(stack)
    end

    nverts = assign_vertex_ids!(m)
    vertex_labels = fill(typemin(Int32), nverts)
    for (label, edge) in zip(corner_labels, corner_edges)
        vid = Int(m.edges[Int(edge)].vertex) + 1
        if vertex_labels[vid] == typemin(Int32)
            vertex_labels[vid] = label
        end
    end

    edge_u, edge_v, edge_color, root_edge, halfedge_to_undirected = extract_undirected_edges(m, root_halfedge)
    faces, face_edge_ids = extract_faces_with_edge_ids(m, halfedge_to_undirected)

    return UniformMap(
        dyck,
        match,
        labels,
        corner_labels,
        vertex_labels,
        faces,
        face_edge_ids,
        edge_u,
        edge_v,
        edge_color,
        root_edge,
        construction_step,
        construction_match,
        construction_edge_label,
        construction_cur_label,
        construction_next_label,
        construction_stack_height,
        construction_case,
    )
end

function generate_uniform_map(; faces::Integer, seed::Integer=1)
    dyck, matching, edge_labels = sample_schaeffer_data(faces; seed=seed)
    return build_uniform_quadrangulation_map(dyck, matching, edge_labels)
end

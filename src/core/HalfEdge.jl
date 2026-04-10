# Minimal half-edge structure used by the map constructors.

mutable struct HalfEdge
    next::Int32
    prev::Int32
    adj::Int32
    color::Any
    vertex::Int32
end

HalfEdge() = HalfEdge(0, 0, 0, nothing, Int32(-1))

mutable struct PlanarMap
    edges::Vector{HalfEdge}
    root::Int32
end

PlanarMap() = PlanarMap(HalfEdge[], Int32(0))

function new_edge!(m::PlanarMap)::Int32
    push!(m.edges, HalfEdge())
    return Int32(length(m.edges))
end

function new_double_edge!(m::PlanarMap; color=nothing)::Int32
    a = Int(new_edge!(m))
    b = Int(new_edge!(m))
    m.edges[a].adj = Int32(b)
    m.edges[b].adj = Int32(a)

    m.edges[a].next = Int32(b)
    m.edges[a].prev = Int32(b)
    m.edges[b].next = Int32(a)
    m.edges[b].prev = Int32(a)

    m.edges[a].color = color
    m.edges[b].color = color
    return Int32(a)
end

function get_next(m::PlanarMap, edge::Integer, n::Integer=1)::Int32
    cur = Int(edge)
    if n >= 0
        for _ in 1:n
            cur = Int(m.edges[cur].next)
        end
    else
        for _ in 1:(-n)
            cur = Int(m.edges[cur].prev)
        end
    end
    return Int32(cur)
end

function get_previous(m::PlanarMap, edge::Integer, n::Integer=1)::Int32
    cur = Int(edge)
    if n >= 0
        for _ in 1:n
            cur = Int(m.edges[cur].prev)
        end
    else
        for _ in 1:(-n)
            cur = Int(m.edges[cur].next)
        end
    end
    return Int32(cur)
end

rotate_cw(m::PlanarMap, edge::Integer) = m.edges[Int(m.edges[Int(edge)].adj)].next

function make_polygon!(m::PlanarMap, n::Integer; color="outer")
    n > 0 || throw(ArgumentError("n must be positive"))
    empty!(m.edges)
    m.root = Int32(0)

    for _ in 1:(2*Int(n))
        push!(m.edges, HalfEdge())
    end

    n_int = Int(n)
    for i in 1:n_int
        out = i
        inn = i + n_int

        m.edges[out].adj = Int32(inn)
        m.edges[out].next = Int32(mod1(i + 1, n_int))
        m.edges[out].prev = Int32(mod1(i - 1, n_int))
        m.edges[out].color = color

        prev_in = n_int + mod1(i - 1, n_int)
        next_in = n_int + mod1(i + 1, n_int)
        m.edges[inn].adj = Int32(out)
        m.edges[inn].next = Int32(prev_in)
        m.edges[inn].prev = Int32(next_in)
        m.edges[inn].color = color
    end

    m.root = Int32(1)
    return m
end

function contract_vertices!(m::PlanarMap, edge1::Integer, edge2::Integer)
    e1 = Int(edge1)
    e2 = Int(edge2)
    e1 == e2 && return
    e1_prev = Int(m.edges[e1].prev)
    e2_prev = Int(m.edges[e2].prev)

    m.edges[e1_prev].next = Int32(e2)
    m.edges[e2_prev].next = Int32(e1)
    m.edges[e1].prev = Int32(e2_prev)
    m.edges[e2].prev = Int32(e1_prev)
    return
end

function insert_edge!(m::PlanarMap, edge1::Integer, edge2::Union{Nothing,Integer}=nothing; color=nothing)::Int32
    new_edge_idx = Int(new_double_edge!(m; color=color))
    new_adj_idx = Int(m.edges[new_edge_idx].adj)

    m.edges[new_adj_idx].next = Int32(edge1)
    m.edges[new_edge_idx].prev = m.edges[Int(edge1)].prev
    m.edges[Int(m.edges[Int(edge1)].prev)].next = Int32(new_edge_idx)
    m.edges[Int(edge1)].prev = Int32(new_adj_idx)

    if edge2 !== nothing
        contract_vertices!(m, m.edges[new_edge_idx].next, Int(edge2))
    end
    return Int32(new_edge_idx)
end

function assign_vertex_ids!(m::PlanarMap)::Int
    vid = 0
    for eid in eachindex(m.edges)
        if m.edges[eid].vertex != -1
            continue
        end
        cur = eid
        while true
            m.edges[cur].vertex = Int32(vid)
            cur = Int(rotate_cw(m, cur))
            cur == eid && break
        end
        vid += 1
    end
    return vid
end

function extract_faces(m::PlanarMap)::Vector{Vector{Int32}}
    seen = Set{Int}()
    faces = Vector{Vector{Int32}}()
    for eid in eachindex(m.edges)
        eid in seen && continue
        cur = eid
        face = Int32[]
        while true
            push!(seen, cur)
            push!(face, m.edges[cur].vertex)
            cur = Int(m.edges[cur].next)
            cur == eid && break
        end
        push!(faces, face)
    end
    return faces
end

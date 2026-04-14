# Circle packing layout for triangulated disk maps.

function _circle_packing_boundary(boundary_vertices, num_vertices::Integer)
    boundary = Int.(collect(boundary_vertices))
    !isempty(boundary) || throw(ArgumentError("boundary_vertices may not be empty"))
    minimum(boundary) >= 0 || throw(ArgumentError("boundary vertices must be >= 0"))
    maximum(boundary) < Int(num_vertices) || throw(ArgumentError("boundary vertices must lie in [0, num_vertices)"))
    length(unique(boundary)) == length(boundary) || throw(ArgumentError("boundary vertices must be distinct"))
    length(boundary) >= 3 || throw(ArgumentError("circle packing requires at least three boundary vertices"))
    return boundary
end

function _circle_packing_triangles(num_vertices::Integer; faces=nothing, triangles=nothing)
    tri = if triangles !== nothing
        sanitize_triangles(triangles; drop_degenerate=true, deduplicate=false)
    elseif faces !== nothing
        fan_triangulate_faces(faces; drop_degenerate=true, deduplicate=false)
    else
        throw(ArgumentError("circle packing requires `faces` or `triangles` to describe a disk triangulation"))
    end

    size(tri, 1) > 0 || throw(ArgumentError("circle packing requires at least one interior triangle"))
    minimum(tri) >= 0 || throw(ArgumentError("triangle vertices must be >= 0"))
    maximum(tri) < Int(num_vertices) || throw(ArgumentError("triangle vertices must lie in [0, num_vertices)"))
    return tri, (triangles !== nothing ? "provided" : "faces_fan")
end

function _boundary_cycle_edges(boundary_vertices)
    boundary = Int32.(collect(boundary_vertices))
    m = length(boundary)
    edges = Matrix{Int32}(undef, m, 2)
    for i in 1:m
        j = i == m ? 1 : i + 1
        edges[i, 1] = boundary[i]
        edges[i, 2] = boundary[j]
    end
    return edges
end

function _triangle_edge_array(triangles)
    tri = sanitize_triangles(triangles; drop_degenerate=true, deduplicate=false)
    size(tri, 1) == 0 && return _empty_edges()

    raw = Matrix{Int32}(undef, 3 * size(tri, 1), 2)
    for i in 1:size(tri, 1)
        a = tri[i, 1]
        b = tri[i, 2]
        c = tri[i, 3]
        row = 3i - 2
        raw[row, 1] = a
        raw[row, 2] = b
        raw[row+1, 1] = b
        raw[row+1, 2] = c
        raw[row+2, 1] = c
        raw[row+2, 2] = a
    end
    return first(collapse_undirected_edges(raw; drop_loops=true))
end

function _validate_disk_triangulation(num_vertices::Integer, triangles, boundary_vertices; triangle_edge_ids=nothing)
    tri = sanitize_triangles(triangles; drop_degenerate=true, deduplicate=false)
    ordered_neighbors = triangle_edge_ids === nothing ?
                        _ordered_neighbors_from_triangles(num_vertices, tri, boundary_vertices) :
                        _ordered_neighbors_from_triangle_edge_ids(num_vertices, tri, triangle_edge_ids, boundary_vertices)
    source = triangle_edge_ids === nothing ? "triangle_vertices" : "triangle_edge_ids"
    return tri, ordered_neighbors, source
end

@inline _clamped_acos(x::Real) = acos(clamp(float(x), -1.0, 1.0))
@inline _clamped_asin(x::Real) = asin(clamp(float(x), -1.0, 1.0))

function _triangle_angle(x::Real, y::Real, z::Real)
    x > 0 || return π
    z > 0 || return π
    return _clamped_acos((x * x + z * z - y * y) / (2 * x * z))
end

function _boundary_angle(r_current::Real, r_prev::Real, outer_radius::Real)
    outer = float(outer_radius)
    rc = float(r_current)
    rp = float(r_prev)
    outer <= rc + rp && return Inf
    denom_current = outer / rc - 1
    denom_prev = outer / rp - 1
    (denom_current > 0 && denom_prev > 0) || return Inf
    return _clamped_acos(1 - 2 / (denom_current * denom_prev))
end

function _boundary_angle_sum_error(outer_radius::Real, boundary_vertices, radii::AbstractVector{<:Real})
    boundary = Int.(collect(boundary_vertices))
    total = 0.0
    for i in eachindex(boundary)
        j = i == length(boundary) ? 1 : i + 1
        total += _boundary_angle(radii[boundary[i]+1], radii[boundary[j]+1], outer_radius)
    end
    return total - 2π
end

function _find_outer_radius(boundary_vertices, radii::AbstractVector{<:Real}; tol::Float64=1.0e-12, maxiter::Int=96)
    boundary = Int.(collect(boundary_vertices))
    lower = 0.0
    for i in eachindex(boundary)
        j = i == length(boundary) ? 1 : i + 1
        lower = max(lower, float(radii[boundary[i]+1]) + float(radii[boundary[j]+1]))
    end
    lower += max(tol, eps(Float64))

    upper = max(2lower, 1.0)
    err_upper = _boundary_angle_sum_error(upper, boundary, radii)
    expand_steps = 0
    while (!isfinite(err_upper) || err_upper > 0.0) && expand_steps < 128
        upper *= 2
        err_upper = _boundary_angle_sum_error(upper, boundary, radii)
        expand_steps += 1
    end
    (isfinite(err_upper) && err_upper <= 0.0) || throw(ArgumentError("failed to bracket the outer circle radius for circle packing"))

    lo = lower
    hi = upper
    for _ in 1:maxiter
        mid = 0.5 * (lo + hi)
        err_mid = _boundary_angle_sum_error(mid, boundary, radii)
        if !isfinite(err_mid) || err_mid > 0.0
            lo = mid
        else
            hi = mid
        end
        hi - lo <= tol * max(1.0, hi) && break
    end
    return hi
end

function _normalize_circle_radii!(radii::AbstractVector{<:Real}, boundary_vertices; min_radius::Float64=1.0e-8)
    @inbounds for i in eachindex(radii)
        radii[i] = max(float(radii[i]), min_radius)
    end
    outer = _find_outer_radius(boundary_vertices, radii)
    radii ./= outer
    return radii
end

function _boundary_cycle_feasible(radii::AbstractVector{<:Real}, boundary_vertices; margin::Float64=1.0e-12)
    boundary = Int.(collect(boundary_vertices))
    for i in eachindex(boundary)
        j = i == length(boundary) ? 1 : i + 1
        if !(float(radii[boundary[i]+1]) + float(radii[boundary[j]+1]) < 1.0 - margin)
            return false
        end
    end
    return true
end

function _ordered_neighbors_from_triangles(num_vertices::Integer, triangles, boundary_vertices)
    n = Int(num_vertices)
    tri = sanitize_triangles(triangles; drop_degenerate=true, deduplicate=false)
    boundary = Int.(collect(boundary_vertices))
    boundary_mask = falses(n)
    boundary_mask[boundary.+1] .= true
    boundary_pos = Dict(v => i for (i, v) in enumerate(boundary))

    wedges = [Vector{NTuple{2,Int32}}() for _ in 1:n]

    for i in 1:size(tri, 1)
        a = tri[i, 1]
        b = tri[i, 2]
        c = tri[i, 3]
        push!(wedges[Int(a)+1], (b, c))
        push!(wedges[Int(b)+1], (c, a))
        push!(wedges[Int(c)+1], (a, b))
    end

    orders = Vector{Vector{Int32}}(undef, n)
    for vertex in 0:(n-1)
        local_wedges = wedges[vertex+1]
        if isempty(local_wedges)
            orders[vertex+1] = Int32[]
            continue
        end

        local_adj = Dict{Int32,Vector{Tuple{Int32,Int}}}()
        local_deg = Dict{Int32,Int}()
        for (edge_id, (u, w)) in enumerate(local_wedges)
            push!(get!(local_adj, u, Tuple{Int32,Int}[]), (w, edge_id))
            push!(get!(local_adj, w, Tuple{Int32,Int}[]), (u, edge_id))
            local_deg[u] = get(local_deg, u, 0) + 1
            local_deg[w] = get(local_deg, w, 0) + 1
        end
        for values in values(local_adj)
            sort!(values; by=x -> (x[1], x[2]))
        end

        odd_vertices = sort!(Int32[v for (v, d) in local_deg if isodd(d)])

        component_count = 0
        visited_neighbors = Set{Int32}()
        for start in keys(local_adj)
            start in visited_neighbors && continue
            component_count += 1
            queue = Int32[start]
            push!(visited_neighbors, start)
            while !isempty(queue)
                current = popfirst!(queue)
                for (next_vertex, _) in get(local_adj, current, Tuple{Int32,Int}[])
                    if !(next_vertex in visited_neighbors)
                        push!(visited_neighbors, next_vertex)
                        push!(queue, next_vertex)
                    end
                end
            end
        end
        component_count == 1 || throw(ArgumentError(
            "circle packing currently requires a simply connected disk triangulation with one cyclic local sector per vertex; vertex $vertex induces $component_count disconnected sectors, which usually means the input has holes or needs edge-identity-aware triangulation data",
        ))

        function eulerian_trail(start::Int32)
            used = falses(length(local_wedges))
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
                    next_vertex, edge_id = pop!(adjlist)
                    used[edge_id] = true
                    push!(stack, next_vertex)
                end
            end

            all(used) || throw(ArgumentError("failed to recover a connected neighbor ordering for circle packing around vertex $vertex"))
            reverse!(path)
            return path
        end

        if boundary_mask[vertex+1]
            idx = boundary_pos[vertex]
            prev_idx = idx == 1 ? length(boundary) : idx - 1
            next_idx = idx == length(boundary) ? 1 : idx + 1
            prev_boundary = Int32(boundary[prev_idx])
            next_boundary = Int32(boundary[next_idx])
            haskey(local_adj, prev_boundary) || throw(ArgumentError("boundary cycle is incompatible with the circle-packing triangulation"))
            haskey(local_adj, next_boundary) || throw(ArgumentError("boundary cycle is incompatible with the circle-packing triangulation"))
            odd_vertices == sort(Int32[prev_boundary, next_boundary]) ||
                throw(ArgumentError(
                    "circle packing expects each boundary vertex to lie on a single disk boundary arc; boundary vertex $vertex has local odd endpoints $(Tuple(Int.(odd_vertices))) but expected $(Tuple(sort(Int32[prev_boundary, next_boundary])))",
                ))

            order = eulerian_trail(prev_boundary)
            order[end] == next_boundary || throw(ArgumentError("boundary neighbor order endpoints do not match the boundary cycle"))
            length(order) == length(local_wedges) + 1 || throw(ArgumentError("boundary neighbor order has the wrong length"))
            orders[vertex+1] = order
        else
            isempty(odd_vertices) || throw(ArgumentError("circle packing expects an interior disk vertex to induce an Eulerian cycle, but vertex $vertex has odd endpoints $(Tuple(Int.(odd_vertices)))"))
            start = minimum(collect(keys(local_adj)))
            cycle = eulerian_trail(start)
            length(cycle) >= 2 && cycle[1] == cycle[end] || throw(ArgumentError("failed to close an interior neighbor cycle"))
            order = cycle[1:end-1]
            length(order) == length(local_wedges) || throw(ArgumentError("interior neighbor order has the wrong length"))
            orders[vertex+1] = order
        end
    end

    return orders
end

function _ordered_neighbors_from_triangle_edge_ids(num_vertices::Integer, triangles, triangle_edge_ids, boundary_vertices)
    n = Int(num_vertices)
    tri = sanitize_triangles(triangles; drop_degenerate=true, deduplicate=false)
    tri_edge_ids = Int32.(triangle_edge_ids)
    size(tri_edge_ids, 1) == size(tri, 1) || throw(ArgumentError("triangle_edge_ids must align with the triangle rows"))
    size(tri_edge_ids, 2) == 3 || throw(ArgumentError("triangle_edge_ids must have three columns ordered as (ab, bc, ca)"))

    boundary = Int.(collect(boundary_vertices))
    boundary_mask = falses(n)
    boundary_mask[boundary.+1] .= true
    boundary_pos = Dict(v => i for (i, v) in enumerate(boundary))

    local_adj = [Dict{Int32,Vector{Tuple{Int32,Int}}}() for _ in 1:n]
    local_deg = [Dict{Int32,Int}() for _ in 1:n]
    local_neighbor = [Dict{Int32,Int32}() for _ in 1:n]
    local_wedge_count = zeros(Int, n)

    function register_neighbor!(vertex_idx::Int, edge_id::Int32, neighbor::Int32)
        existing = get(local_neighbor[vertex_idx], edge_id, Int32(-1))
        if existing != Int32(-1) && existing != neighbor
            throw(ArgumentError("triangle_edge_ids are inconsistent with the triangulation around vertex $(vertex_idx - 1)"))
        end
        local_neighbor[vertex_idx][edge_id] = neighbor
    end

    function add_local_wedge!(vertex::Int32, left_edge::Int32, left_neighbor::Int32, right_edge::Int32, right_neighbor::Int32)
        vidx = Int(vertex) + 1
        register_neighbor!(vidx, left_edge, left_neighbor)
        register_neighbor!(vidx, right_edge, right_neighbor)
        local_wedge_count[vidx] += 1
        wedge_id = local_wedge_count[vidx]
        push!(get!(local_adj[vidx], left_edge, Tuple{Int32,Int}[]), (right_edge, wedge_id))
        push!(get!(local_adj[vidx], right_edge, Tuple{Int32,Int}[]), (left_edge, wedge_id))
        local_deg[vidx][left_edge] = get(local_deg[vidx], left_edge, 0) + 1
        local_deg[vidx][right_edge] = get(local_deg[vidx], right_edge, 0) + 1
    end

    for i in 1:size(tri, 1)
        a = tri[i, 1]
        b = tri[i, 2]
        c = tri[i, 3]
        eab = tri_edge_ids[i, 1]
        ebc = tri_edge_ids[i, 2]
        eca = tri_edge_ids[i, 3]
        add_local_wedge!(a, eca, c, eab, b)
        add_local_wedge!(b, eab, a, ebc, c)
        add_local_wedge!(c, ebc, b, eca, a)
    end

    orders = Vector{Vector{Int32}}(undef, n)
    for vertex in 0:(n-1)
        local_adj_v = local_adj[vertex+1]
        if isempty(local_adj_v)
            orders[vertex+1] = Int32[]
            continue
        end

        for values in values(local_adj_v)
            sort!(values; by=x -> (x[1], x[2]))
        end

        odd_edges = sort!(Int32[eid for (eid, d) in local_deg[vertex+1] if isodd(d)])

        component_count = 0
        components = Vector{Vector{Int32}}()
        visited_edges = Set{Int32}()
        for start in keys(local_adj_v)
            start in visited_edges && continue
            component_count += 1
            component = Int32[]
            queue = Int32[start]
            push!(visited_edges, start)
            while !isempty(queue)
                current = popfirst!(queue)
                push!(component, current)
                for (next_edge, _) in get(local_adj_v, current, Tuple{Int32,Int}[])
                    if !(next_edge in visited_edges)
                        push!(visited_edges, next_edge)
                        push!(queue, next_edge)
                    end
                end
            end
            push!(components, component)
        end

        function eulerian_trail(local_adj_component::Dict{Int32,Vector{Tuple{Int32,Int}}}, start::Int32)
            wedge_ids = sort!(unique(Int[wedge_id for values in values(local_adj_component) for (_, wedge_id) in values]))
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

            all(used[wedge_ids]) || throw(ArgumentError("failed to recover a connected neighbor ordering for circle packing around vertex $vertex"))
            reverse!(path)
            return path
        end

        neighbor_of = local_neighbor[vertex+1]
        function stitch_component_orders(component_orders::Vector{Vector{Int32}})
            assembled = copy(component_orders[1])
            pending = [copy(order) for order in component_orders[2:end]]

            while !isempty(pending)
                attached = false
                for i in eachindex(pending)
                    order = pending[i]
                    if assembled[end] == order[1]
                        assembled = vcat(assembled, order)
                    elseif assembled[end] == order[end]
                        assembled = vcat(assembled, reverse(order))
                    elseif assembled[1] == order[end]
                        assembled = vcat(order, assembled)
                    elseif assembled[1] == order[1]
                        assembled = vcat(reverse(order), assembled)
                    else
                        continue
                    end
                    deleteat!(pending, i)
                    attached = true
                    break
                end
                attached || return nothing
            end

            return assembled
        end

        if boundary_mask[vertex+1]
            idx = boundary_pos[vertex]
            prev_idx = idx == 1 ? length(boundary) : idx - 1
            next_idx = idx == length(boundary) ? 1 : idx + 1
            prev_boundary = Int32(boundary[prev_idx])
            next_boundary = Int32(boundary[next_idx])
            if component_count == 1
                length(odd_edges) == 2 || throw(ArgumentError("boundary cycle is incompatible with the edge-identity-aware circle-packing triangulation"))

                odd_neighbors = sort!(Int32[neighbor_of[eid] for eid in odd_edges])
                odd_neighbors == sort(Int32[prev_boundary, next_boundary]) ||
                    throw(ArgumentError(
                        "circle packing expects each boundary vertex to lie on a single disk boundary arc; boundary vertex $vertex has local odd endpoints $(Tuple(Int.(odd_neighbors))) but expected $(Tuple(sort(Int32[prev_boundary, next_boundary])))",
                    ))

                start_edge = if neighbor_of[odd_edges[1]] == prev_boundary
                    odd_edges[1]
                else
                    odd_edges[2]
                end

                edge_path = eulerian_trail(local_adj_v, start_edge)
                neighbor_of[edge_path[end]] == next_boundary || throw(ArgumentError("boundary neighbor order endpoints do not match the boundary cycle"))
                length(edge_path) == local_wedge_count[vertex+1] + 1 || throw(ArgumentError("boundary neighbor order has the wrong length"))
                orders[vertex+1] = Int32[neighbor_of[eid] for eid in edge_path]
            else
                component_orders = Vector{Vector{Int32}}()
                for component in components
                    comp_adj = Dict{Int32,Vector{Tuple{Int32,Int}}}()
                    comp_set = Set(component)
                    comp_odd = Int32[]
                    for edge_id in component
                        values = [(next_edge, wedge_id) for (next_edge, wedge_id) in get(local_adj_v, edge_id, Tuple{Int32,Int}[]) if next_edge in comp_set]
                        comp_adj[edge_id] = values
                        isodd(length(values)) && push!(comp_odd, edge_id)
                    end
                    length(comp_odd) == 2 || throw(ArgumentError(
                        "circle packing currently requires a simply connected disk triangulation with one cyclic local sector per vertex; vertex $vertex induces $component_count disconnected sectors, which usually means the input has holes or needs edge-identity-aware triangulation data",
                    ))
                    path = eulerian_trail(comp_adj, comp_odd[1])
                    order = Int32[neighbor_of[eid] for eid in path]
                    if order[1] != neighbor_of[comp_odd[1]]
                        reverse!(order)
                    end
                    push!(component_orders, order)
                end

                start_component = findfirst(order -> order[1] == prev_boundary || order[end] == prev_boundary, component_orders)
                start_component === nothing && throw(ArgumentError(
                    "circle packing currently requires a simply connected disk triangulation with one cyclic local sector per vertex; vertex $vertex induces $component_count disconnected sectors, which usually means the input has holes or needs edge-identity-aware triangulation data",
                ))
                if component_orders[start_component][end] == prev_boundary
                    reverse!(component_orders[start_component])
                end
                if start_component != 1
                    component_orders[1], component_orders[start_component] = component_orders[start_component], component_orders[1]
                end

                stitched = stitch_component_orders(component_orders)
                stitched === nothing && throw(ArgumentError(
                    "circle packing currently requires a simply connected disk triangulation with one cyclic local sector per vertex; vertex $vertex induces $component_count disconnected sectors, which usually means the input has holes or needs edge-identity-aware triangulation data",
                ))
                stitched[1] == prev_boundary || reverse!(stitched)
                stitched[1] == prev_boundary && stitched[end] == next_boundary || throw(ArgumentError("boundary neighbor order endpoints do not match the boundary cycle"))
                length(stitched) == length(local_deg[vertex+1]) || throw(ArgumentError("boundary neighbor order has the wrong length"))
                orders[vertex+1] = stitched
            end
        else
            if component_count == 1
                if isempty(odd_edges)
                    start_edge = minimum(collect(keys(local_adj_v)))
                    edge_cycle = eulerian_trail(local_adj_v, start_edge)
                    length(edge_cycle) >= 2 && edge_cycle[1] == edge_cycle[end] || throw(ArgumentError("failed to close an interior neighbor cycle"))
                    order_edges = edge_cycle[1:end-1]
                    length(order_edges) == local_wedge_count[vertex+1] || throw(ArgumentError("interior neighbor order has the wrong length"))
                    orders[vertex+1] = Int32[neighbor_of[eid] for eid in order_edges]
                elseif length(odd_edges) == 2 && neighbor_of[odd_edges[1]] == neighbor_of[odd_edges[2]]
                    edge_path = eulerian_trail(local_adj_v, odd_edges[1])
                    edge_path[end] == odd_edges[2] || throw(ArgumentError("failed to recover an interior repeated-neighbor cycle"))
                    length(edge_path) == local_wedge_count[vertex+1] + 1 || throw(ArgumentError("interior neighbor order has the wrong length"))
                    orders[vertex+1] = Int32[neighbor_of[eid] for eid in edge_path]
                else
                    throw(ArgumentError("circle packing expects an interior disk vertex to induce an Eulerian cycle, but vertex $vertex has odd endpoints $(Tuple(Int.(odd_edges)))"))
                end
            else
                component_orders = Vector{Vector{Int32}}()
                for component in components
                    comp_adj = Dict{Int32,Vector{Tuple{Int32,Int}}}()
                    comp_set = Set(component)
                    comp_odd = Int32[]
                    for edge_id in component
                        values = [(next_edge, wedge_id) for (next_edge, wedge_id) in get(local_adj_v, edge_id, Tuple{Int32,Int}[]) if next_edge in comp_set]
                        comp_adj[edge_id] = values
                        isodd(length(values)) && push!(comp_odd, edge_id)
                    end
                    length(comp_odd) == 2 || throw(ArgumentError(
                        "circle packing currently requires a simply connected disk triangulation with one cyclic local sector per vertex; vertex $vertex induces $component_count disconnected sectors, which usually means the input has holes or needs edge-identity-aware triangulation data",
                    ))
                    path = eulerian_trail(comp_adj, comp_odd[1])
                    push!(component_orders, Int32[neighbor_of[eid] for eid in path])
                end

                stitched = stitch_component_orders(component_orders)
                stitched === nothing && throw(ArgumentError(
                    "circle packing currently requires a simply connected disk triangulation with one cyclic local sector per vertex; vertex $vertex induces $component_count disconnected sectors, which usually means the input has holes or needs edge-identity-aware triangulation data",
                ))
                length(stitched) == length(local_deg[vertex+1]) || throw(ArgumentError("interior neighbor order has the wrong length"))
                orders[vertex+1] = stitched
            end
        end
    end

    return orders
end

function _update_boundary_centers!(positions::AbstractMatrix{<:Real}, radii::AbstractVector{<:Real}, boundary_vertices)
    boundary = Int.(collect(boundary_vertices))
    fill!(positions, 0.0)

    first_vertex = boundary[1] + 1
    positions[first_vertex, 1] = 1.0 - radii[first_vertex]
    positions[first_vertex, 2] = 0.0

    angle = 0.0
    for i in 2:length(boundary)
        vertex = boundary[i] + 1
        prev_vertex = boundary[i-1] + 1
        angle += _boundary_angle(radii[vertex], radii[prev_vertex], 1.0)
        radius_to_center = 1.0 - radii[vertex]
        positions[vertex, 1] = radius_to_center * cos(angle)
        positions[vertex, 2] = radius_to_center * sin(angle)
    end
    return positions
end

function _boundary_mask(boundary_vertices, n::Integer)
    mask = falses(Int(n))
    mask[Int.(collect(boundary_vertices)).+1] .= true
    return mask
end

function _boundary_neighbor_indices(boundary_vertices, n::Integer)
    prev_idx = zeros(Int, Int(n))
    next_idx = zeros(Int, Int(n))
    boundary = Int.(collect(boundary_vertices))
    m = length(boundary)
    for i in 1:m
        v = boundary[i] + 1
        prev_idx[v] = boundary[i == 1 ? m : i - 1] + 1
        next_idx[v] = boundary[i == m ? 1 : i + 1] + 1
    end
    return prev_idx, next_idx
end

function _neighbor_index_lists(ordered_neighbors)
    out = Vector{Vector{Int}}(undef, length(ordered_neighbors))
    for i in eachindex(ordered_neighbors)
        out[i] = Int.(ordered_neighbors[i]) .+ 1
    end
    return out
end

function _initialize_circle_packing_radii(
    ordered_neighbors_idx,
    boundary;
    initial_radius::Float64=0.5,
    min_radius::Float64=1.0e-8,
)
    n = length(ordered_neighbors_idx)
    boundary_mask = _boundary_mask(boundary, n)
    radii = fill(float(initial_radius), n)
    @inbounds for i in eachindex(radii)
        boundary_mask[i] && continue
        deg = max(length(ordered_neighbors_idx[i]), 3)
        s = sin(π / deg)
        guess = s <= 0.0 ? initial_radius : (1.0 / s - 1.0)
        radii[i] = max(float(guess), min_radius)
    end
    _normalize_circle_radii!(radii, boundary; min_radius=min_radius)
    return radii, boundary_mask
end

function _circle_packing_residual_at_vertex(
    i::Integer,
    radii::AbstractVector{<:Real},
    ordered_neighbors_idx,
    boundary_mask::AbstractVector{Bool},
    boundary_prev_idx::AbstractVector{Int},
    boundary_next_idx::AbstractVector{Int},
)
    nbd = ordered_neighbors_idx[i]
    isempty(nbd) && return 0.0

    if boundary_mask[i]
        length(nbd) >= 2 || return 0.0
        angle_sum = 0.0
        for j in 1:(length(nbd)-1)
            u = nbd[j]
            v = nbd[j+1]
            angle_sum += _triangle_angle(
                radii[i] + radii[u],
                radii[u] + radii[v],
                radii[i] + radii[v],
            )
        end

        prev_idx = boundary_prev_idx[i]
        next_idx = boundary_next_idx[i]
        prev_idx > 0 || return Inf
        next_idx > 0 || return Inf
        gap_i = 1.0 - float(radii[i])
        gap_prev = 1.0 - float(radii[prev_idx])
        gap_next = 1.0 - float(radii[next_idx])
        (gap_i > 0.0 && gap_prev > 0.0 && gap_next > 0.0) || return Inf

        beta_prev = _triangle_angle(gap_i, gap_prev, radii[i] + radii[prev_idx])
        beta_next = _triangle_angle(gap_i, gap_next, radii[i] + radii[next_idx])
        return angle_sum - beta_prev - beta_next
    end

    angle_sum = 0.0
    deg = length(nbd)
    for j in 1:deg
        u = nbd[j]
        v = nbd[j == deg ? 1 : j + 1]
        angle_sum += _triangle_angle(
            radii[i] + radii[u],
            radii[u] + radii[v],
            radii[i] + radii[v],
        )
    end
    return angle_sum - 2π
end

function _circle_packing_residuals!(
    residual::AbstractVector{<:Real},
    radii::AbstractVector{<:Real},
    ordered_neighbors_idx,
    boundary_mask::AbstractVector{Bool},
    boundary_prev_idx::AbstractVector{Int},
    boundary_next_idx::AbstractVector{Int},
)
    @inbounds for i in eachindex(residual)
        residual[i] = _circle_packing_residual_at_vertex(
            i,
            radii,
            ordered_neighbors_idx,
            boundary_mask,
            boundary_prev_idx,
            boundary_next_idx,
        )
    end
    return residual
end

function _circle_packing_jacobian_workspace(ordered_neighbors_idx)
    n = length(ordered_neighbors_idx)
    rows = Int[]
    cols = Int[]
    vals = Float64[]
    affected_rows = Vector{Vector{Int}}(undef, n)

    for j in 1:n
        local_rows = sort!(unique(vcat([j], ordered_neighbors_idx[j])))
        affected_rows[j] = local_rows
        for row in local_rows
            push!(rows, row)
            push!(cols, j)
            push!(vals, 0.0)
            push!(rows, j)
            push!(cols, row)
            push!(vals, 0.0)
        end
    end

    J = sparse(rows, cols, vals, n, n)
    pair_to_nz = Dict{Tuple{Int,Int},Int}()
    for col in 1:n
        for ptr in nzrange(J, col)
            pair_to_nz[(J.rowval[ptr], col)] = ptr
        end
    end

    jacobian_slots = Vector{Vector{Int}}(undef, n)
    diag_slots = zeros(Int, n)
    for j in 1:n
        jacobian_slots[j] = [pair_to_nz[(row, j)] for row in affected_rows[j]]
        diag_slots[j] = pair_to_nz[(j, j)]
    end

    symmetric_slots = zeros(Int, nnz(J))
    for col in 1:n
        for ptr in nzrange(J, col)
            symmetric_slots[ptr] = get(pair_to_nz, (col, J.rowval[ptr]), ptr)
        end
    end

    return J, affected_rows, jacobian_slots, diag_slots, symmetric_slots
end

function _symmetrize_sparse_values!(A, symmetric_slots::AbstractVector{Int})
    @inbounds for ptr in eachindex(A.nzval)
        partner = symmetric_slots[ptr]
        partner < ptr && continue
        partner == ptr && continue
        value = 0.5 * (A.nzval[ptr] + A.nzval[partner])
        A.nzval[ptr] = value
        A.nzval[partner] = value
    end
    return A
end

function _fill_circle_packing_jacobian!(
    J,
    radii::AbstractVector{<:Real},
    residual::AbstractVector{<:Real},
    ordered_neighbors_idx,
    boundary_mask::AbstractVector{Bool},
    boundary_prev_idx::AbstractVector{Int},
    boundary_next_idx::AbstractVector{Int},
    affected_rows,
    jacobian_slots;
    log_eps::Float64=1.0e-7,
)
    fill!(J.nzval, 0.0)
    scale = exp(log_eps)

    for j in eachindex(radii)
        r_old = float(radii[j])
        radii[j] = r_old * scale
        rows = affected_rows[j]
        slots = jacobian_slots[j]
        for k in eachindex(rows)
            row = rows[k]
            perturbed = _circle_packing_residual_at_vertex(
                row,
                radii,
                ordered_neighbors_idx,
                boundary_mask,
                boundary_prev_idx,
                boundary_next_idx,
            )
            J.nzval[slots[k]] = (perturbed - residual[row]) / log_eps
        end
        radii[j] = r_old
    end

    return J
end

@inline function _row_distance(positions::AbstractMatrix{<:Real}, i::Integer, j::Integer)
    acc = 0.0
    @inbounds for k in axes(positions, 2)
        δ = float(positions[i, k]) - float(positions[j, k])
        acc += δ * δ
    end
    return sqrt(acc)
end

function _edge_conductance(r_center::Real, r_current::Real, r_prev::Real, r_next::Real)
    value = (
        sqrt(max(0.0, float(r_center) * float(r_prev) * float(r_current) / (float(r_center) + float(r_prev) + float(r_current)))) +
        sqrt(max(0.0, float(r_center) * float(r_next) * float(r_current) / (float(r_center) + float(r_next) + float(r_current))))
    ) / (float(r_center) + float(r_current))
    return isfinite(value) ? value : 0.0
end

@inline function _circle_packing_cg_statevars(n::Integer)
    m = Int(n)
    return IterativeSolvers.CGStateVariables(
        zeros(Float64, m),
        zeros(Float64, m),
        zeros(Float64, m),
    )
end

@inline _circle_packing_cg_maxiter(n::Integer; cap::Int=512) = min(max(64, 4 * ceil(Int, sqrt(max(Int(n), 1)))), cap)

mutable struct _CirclePackingNewtonCGWorkspace
    preconditioner::Any
    statevars::Any
    step::Vector{Float64}
end

mutable struct _CirclePackingPositionCGWorkspace
    preconditioner::Any
    x_statevars::Any
    y_statevars::Any
    x::Vector{Float64}
    y::Vector{Float64}
end

function _circle_packing_newton_cg_workspace(n::Integer)
    return _CirclePackingNewtonCGWorkspace(
        nothing,
        _circle_packing_cg_statevars(n),
        zeros(Float64, Int(n)),
    )
end

function _circle_packing_position_cg_workspace(n::Integer)
    return _CirclePackingPositionCGWorkspace(
        nothing,
        _circle_packing_cg_statevars(n),
        _circle_packing_cg_statevars(n),
        zeros(Float64, Int(n)),
        zeros(Float64, Int(n)),
    )
end

function _refresh_ilu0_preconditioner(preconditioner, A::SparseMatrixCSC{Float64,Int})
    try
        return preconditioner === nothing ? ILUZero.ilu0(A) : ILUZero.ilu0!(preconditioner, A)
    catch
        return nothing
    end
end

function _circle_packing_cg_solve!(
    x::Vector{Float64},
    A::SparseMatrixCSC{Float64,Int},
    b::AbstractVector{Float64},
    statevars;
    preconditioner=nothing,
    reltol::Float64=1.0e-6,
    abstol::Float64=0.0,
    maxiter::Int=size(A, 1),
)
    fill!(x, 0.0)
    if preconditioner !== nothing
        try
            _, history = IterativeSolvers.cg!(
                x,
                A,
                b;
                Pl=preconditioner,
                reltol=reltol,
                abstol=abstol,
                maxiter=maxiter,
                log=true,
                statevars=statevars,
                initially_zero=true,
            )
            all(isfinite, x) && return x, history, true
        catch
        end
    end

    fill!(x, 0.0)
    try
        _, history = IterativeSolvers.cg!(
            x,
            A,
            b;
            reltol=reltol,
            abstol=abstol,
            maxiter=maxiter,
            log=true,
            statevars=statevars,
            initially_zero=true,
        )
        return x, history, false
    catch
        return x, nothing, false
    end
end

function _conductance_workspace(ordered_neighbors_idx, boundary_mask::AbstractVector{Bool})
    n = length(ordered_neighbors_idx)
    interior_index = zeros(Int, n)
    interior_vertices = Int[]
    for i in 1:n
        if !boundary_mask[i] && !isempty(ordered_neighbors_idx[i])
            push!(interior_vertices, i)
            interior_index[i] = length(interior_vertices)
        end
    end

    nin = length(interior_vertices)
    rows = Int[]
    cols = Int[]
    vals = Float64[]
    insertion_slots = [Int[] for _ in 1:n]
    diag_insertions = zeros(Int, n)

    sizehint!(rows, max(4nin, 16))
    sizehint!(cols, max(4nin, 16))
    sizehint!(vals, max(4nin, 16))

    for i in interior_vertices
        row = interior_index[i]
        nbd = ordered_neighbors_idx[i]
        local_insertions = zeros(Int, length(nbd))
        for j in eachindex(nbd)
            col = interior_index[nbd[j]]
            col == 0 && continue
            push!(rows, row)
            push!(cols, col)
            push!(vals, 0.0)
            local_insertions[j] = length(rows)
        end
        insertion_slots[i] = local_insertions

        push!(rows, row)
        push!(cols, row)
        push!(vals, 1.0)
        diag_insertions[i] = length(rows)
    end

    A = sparse(rows, cols, vals, nin, nin)
    pair_to_nz = Dict{Tuple{Int,Int},Int}()
    for col in 1:nin
        for ptr in nzrange(A, col)
            pair_to_nz[(A.rowval[ptr], col)] = ptr
        end
    end

    neighbor_slots = [Int[] for _ in 1:n]
    diag_slots = zeros(Int, n)
    max_degree = 0
    for i in interior_vertices
        row = interior_index[i]
        local_insertions = insertion_slots[i]
        if isempty(local_insertions)
            diag_slots[i] = pair_to_nz[(row, row)]
            continue
        end

        slots = similar(local_insertions)
        for j in eachindex(local_insertions)
            insert_idx = local_insertions[j]
            slots[j] = insert_idx == 0 ? 0 : pair_to_nz[(rows[insert_idx], cols[insert_idx])]
        end
        neighbor_slots[i] = slots
        diag_slots[i] = pair_to_nz[(row, row)]
        max_degree = max(max_degree, length(slots))
    end

    return A, interior_vertices, interior_index, neighbor_slots, diag_slots, max_degree
end

function _update_conductance_matrix!(
    A::SparseMatrixCSC{Float64,Int},
    rhs::AbstractMatrix{Float64},
    positions::AbstractMatrix{Float64},
    radii::AbstractVector{<:Real},
    ordered_neighbors_idx,
    boundary_mask::AbstractVector{Bool},
    interior_vertices,
    interior_index::AbstractVector{Int},
    neighbor_slots,
    diag_slots,
    scratch::AbstractVector{<:Real},
)
    fill!(A.nzval, 0.0)
    fill!(rhs, 0.0)

    for i in interior_vertices
        nbd = ordered_neighbors_idx[i]
        slots = neighbor_slots[i]
        deg = length(nbd)
        csum = 0.0
        @inbounds for j in 1:deg
            prev = j == 1 ? deg : j - 1
            next = j == deg ? 1 : j + 1
            c = _edge_conductance(
                radii[i],
                radii[nbd[j]],
                radii[nbd[prev]],
                radii[nbd[next]],
            )
            scratch[j] = c
            csum += c
        end

        row = interior_index[i]
        if csum > 0.0
            @inbounds for j in 1:deg
                c = scratch[j]
                nb = nbd[j]
                slot = slots[j]
                if slot == 0
                    rhs[row, 1] += c * float(positions[nb, 1])
                    rhs[row, 2] += c * float(positions[nb, 2])
                else
                    A.nzval[slot] += -c
                end
            end
            A.nzval[diag_slots[i]] = csum
        else
            A.nzval[diag_slots[i]] = 1.0
        end
    end

    return A
end

function _conductance_matrix(radii::AbstractVector{<:Real}, ordered_neighbors, boundary_vertices)
    n = length(radii)
    boundary_mask = falses(n)
    boundary_mask[Int.(collect(boundary_vertices)).+1] .= true

    rows = Int[]
    cols = Int[]
    vals = Float64[]
    sizehint!(rows, max(4n, 16))
    sizehint!(cols, max(4n, 16))
    sizehint!(vals, max(4n, 16))

    for i in 1:n
        if boundary_mask[i] || isempty(ordered_neighbors[i])
            push!(rows, i)
            push!(cols, i)
            push!(vals, 1.0)
            continue
        end

        nbd = ordered_neighbors[i]
        deg = length(nbd)
        cvals = zeros(Float64, deg)
        csum = 0.0
        for j in 1:deg
            prev = j == 1 ? deg : j - 1
            next = j == deg ? 1 : j + 1
            c = _edge_conductance(
                radii[i],
                radii[Int(nbd[j])+1],
                radii[Int(nbd[prev])+1],
                radii[Int(nbd[next])+1],
            )
            cvals[j] = c
            csum += c
        end

        if csum <= 0.0
            cvals .= 1.0 / deg
        else
            cvals ./= csum
        end

        for j in 1:deg
            push!(rows, i)
            push!(cols, Int(nbd[j]) + 1)
            push!(vals, -cvals[j])
        end
        push!(rows, i)
        push!(cols, i)
        push!(vals, 1.0)
    end

    return sparse(rows, cols, vals, n, n)
end

function _update_circle_radii!(radii::AbstractVector{<:Real}, ordered_neighbors, boundary_vertices, positions::AbstractMatrix{<:Real})
    n = length(radii)
    boundary_mask = falses(n)
    boundary_mask[Int.(collect(boundary_vertices)).+1] .= true

    for i in 1:n
        nbd = ordered_neighbors[i]
        isempty(nbd) && continue

        if boundary_mask[i]
            length(nbd) >= 2 || continue
            phi_total = 0.0
            radius_acc = 0.0
            for j in 1:(length(nbd)-1)
                u = Int(nbd[j]) + 1
                v = Int(nbd[j+1]) + 1
                x = norm(view(positions, i, :) .- view(positions, u, :))
                y = norm(view(positions, u, :) .- view(positions, v, :))
                z = norm(view(positions, i, :) .- view(positions, v, :))
                phi = _triangle_angle(x, y, z)
                radius_eff = max(0.0, 0.5 * (x + z - y))
                phi_total += phi
                radius_acc += phi * radius_eff * radius_eff
            end
        else
            phi_total = 2π
            radius_acc = 0.0
            deg = length(nbd)
            for j in 1:deg
                k = j == deg ? 1 : j + 1
                u = Int(nbd[j]) + 1
                v = Int(nbd[k]) + 1
                x = norm(view(positions, i, :) .- view(positions, u, :))
                y = norm(view(positions, u, :) .- view(positions, v, :))
                z = norm(view(positions, i, :) .- view(positions, v, :))
                phi = _triangle_angle(x, y, z)
                radius_eff = max(0.0, 0.5 * (x + z - y))
                radius_acc += phi * radius_eff * radius_eff
            end
        end

        if phi_total > 0.0 && radius_acc > 0.0
            candidate = sqrt(radius_acc / phi_total)
            isfinite(candidate) && candidate > 0.0 && (radii[i] = candidate)
        end
    end
    return radii
end

function _update_circle_radii!(radii::AbstractVector{<:Real}, ordered_neighbors_idx, boundary_mask::AbstractVector{Bool}, positions::AbstractMatrix{<:Real})
    for i in eachindex(radii)
        nbd = ordered_neighbors_idx[i]
        isempty(nbd) && continue

        if boundary_mask[i]
            length(nbd) >= 2 || continue
            phi_total = 0.0
            radius_acc = 0.0
            for j in 1:(length(nbd)-1)
                u = nbd[j]
                v = nbd[j+1]
                x = _row_distance(positions, i, u)
                y = _row_distance(positions, u, v)
                z = _row_distance(positions, i, v)
                phi = _triangle_angle(x, y, z)
                radius_eff = max(0.0, 0.5 * (x + z - y))
                phi_total += phi
                radius_acc += phi * radius_eff * radius_eff
            end
        else
            phi_total = 2π
            radius_acc = 0.0
            deg = length(nbd)
            for j in 1:deg
                k = j == deg ? 1 : j + 1
                u = nbd[j]
                v = nbd[k]
                x = _row_distance(positions, i, u)
                y = _row_distance(positions, u, v)
                z = _row_distance(positions, i, v)
                phi = _triangle_angle(x, y, z)
                radius_eff = max(0.0, 0.5 * (x + z - y))
                radius_acc += phi * radius_eff * radius_eff
            end
        end

        if phi_total > 0.0 && radius_acc > 0.0
            candidate = sqrt(radius_acc / phi_total)
            isfinite(candidate) && candidate > 0.0 && (radii[i] = candidate)
        end
    end
    return radii
end

function _circle_packing_interior_residual(radii::AbstractVector{<:Real}, ordered_neighbors, boundary_vertices)
    boundary_mask = falses(length(radii))
    boundary_mask[Int.(collect(boundary_vertices)).+1] .= true
    residual = 0.0
    for i in eachindex(ordered_neighbors)
        boundary_mask[i] && continue
        nbd = ordered_neighbors[i]
        isempty(nbd) && continue
        angle_sum = 0.0
        deg = length(nbd)
        for j in 1:deg
            prev = j == 1 ? deg : j - 1
            angle_sum += 2 * _clamped_asin(
                sqrt(max(
                    0.0,
                    1 / ((1 + radii[i] / radii[Int(nbd[prev])+1]) * (1 + radii[i] / radii[Int(nbd[j])+1])),
                )),
            )
        end
        residual = max(residual, abs(angle_sum - 2π))
    end
    return residual
end

function _circle_packing_interior_residual(radii::AbstractVector{<:Real}, ordered_neighbors_idx, boundary_mask::AbstractVector{Bool})
    residual = 0.0
    for i in eachindex(ordered_neighbors_idx)
        boundary_mask[i] && continue
        nbd = ordered_neighbors_idx[i]
        isempty(nbd) && continue
        angle_sum = 0.0
        deg = length(nbd)
        for j in 1:deg
            prev = j == 1 ? deg : j - 1
            angle_sum += 2 * _clamped_asin(
                sqrt(max(
                    0.0,
                    1 / ((1 + radii[i] / radii[nbd[prev]]) * (1 + radii[i] / radii[nbd[j]])),
                )),
            )
        end
        residual = max(residual, abs(angle_sum - 2π))
    end
    return residual
end

function _circle_packing_edge_residual(positions::AbstractMatrix{<:Real}, radii::AbstractVector{<:Real}, edges)
    arr = sanitize_edge_array(edges)
    size(arr, 1) == 0 && return 0.0
    residual = 0.0
    for i in 1:size(arr, 1)
        u = Int(arr[i, 1]) + 1
        v = Int(arr[i, 2]) + 1
        dist = norm(view(positions, u, :) .- view(positions, v, :))
        residual = max(residual, abs(dist - radii[u] - radii[v]))
    end
    return residual
end

function _circle_packing_newton_solve!(
    radii::AbstractVector{Float64},
    boundary,
    ordered_neighbors_idx;
    maxiter::Int=200,
    tol::Float64=1.0e-8,
    min_radius::Float64=1.0e-8,
)
    n = length(radii)
    boundary_mask = _boundary_mask(boundary, n)
    boundary_prev_idx, boundary_next_idx = _boundary_neighbor_indices(boundary, n)
    residual = zeros(Float64, n)
    trial_residual = similar(residual)
    u = log.(max.(radii, min_radius))
    trial_u = similar(u)
    trial_radii = similar(radii)
    rhs = similar(residual)
    J, affected_rows, jacobian_slots, diag_slots, symmetric_slots = _circle_packing_jacobian_workspace(ordered_neighbors_idx)
    J_solver = copy(J)
    cg_workspace = _circle_packing_newton_cg_workspace(n)
    converged = false
    iter_count = 0
    max_radius_delta = Inf
    damping = 1.0e-4
    linear_iterations = 0
    cg_maxiter = _circle_packing_cg_maxiter(n; cap=400)

    for iter in 1:maxiter
        _circle_packing_residuals!(residual, radii, ordered_neighbors_idx, boundary_mask, boundary_prev_idx, boundary_next_idx)
        residual_norm = maximum(abs.(residual))
        iter_count = iter
        if !(isfinite(residual_norm))
            break
        elseif residual_norm <= tol
            converged = true
            max_radius_delta = 0.0
            break
        end

        _fill_circle_packing_jacobian!(
            J,
            radii,
            residual,
            ordered_neighbors_idx,
            boundary_mask,
            boundary_prev_idx,
            boundary_next_idx,
            affected_rows,
            jacobian_slots,
        )

        rhs .= -residual
        accepted = false
        λ = damping
        while λ <= 1.0e8
            copyto!(J_solver.nzval, J.nzval)
            _symmetrize_sparse_values!(J_solver, symmetric_slots)
            @inbounds for slot in diag_slots
                J_solver.nzval[slot] += λ
            end
            cg_workspace.preconditioner = _refresh_ilu0_preconditioner(cg_workspace.preconditioner, J_solver)
            step, history, _ = _circle_packing_cg_solve!(
                cg_workspace.step,
                J_solver,
                rhs,
                cg_workspace.statevars;
                preconditioner=cg_workspace.preconditioner,
                reltol=1.0e-3,
                abstol=max(1.0e-10, 0.1 * tol),
                maxiter=cg_maxiter,
            )
            linear_iterations += history === nothing ? 0 : Int(getproperty(history, :iters))
            if history === nothing || !all(isfinite, step)
                λ *= 4.0
                continue
            end

            trial_u .= u .+ step
            @inbounds for i in eachindex(trial_radii)
                trial_radii[i] = max(exp(trial_u[i]), min_radius)
            end
            valid = _boundary_cycle_feasible(trial_radii, boundary)
            @inbounds for i in eachindex(trial_radii)
                if valid && boundary_mask[i] && !(trial_radii[i] < 1.0 - 1.0e-12)
                    valid = false
                    break
                end
            end
            if valid
                _circle_packing_residuals!(trial_residual, trial_radii, ordered_neighbors_idx, boundary_mask, boundary_prev_idx, boundary_next_idx)
                trial_norm = maximum(abs.(trial_residual))
                if isfinite(trial_norm) && trial_norm < residual_norm
                    max_radius_delta = maximum(abs.(trial_radii .- radii) ./ max.(radii, min_radius))
                    radii .= trial_radii
                    u .= trial_u
                    damping = max(1.0e-8, 0.25 * λ)
                    accepted = true
                    break
                end
            end
            λ *= 4.0
        end

        accepted || break
    end

    _circle_packing_residuals!(residual, radii, ordered_neighbors_idx, boundary_mask, boundary_prev_idx, boundary_next_idx)
    return Dict{String,Any}(
        "converged" => converged,
        "iterations" => iter_count,
        "max_radius_delta" => max_radius_delta,
        "residual_inf_norm" => maximum(abs.(residual)),
        "boundary_mask" => boundary_mask,
        "linear_solver" => "cg",
        "linear_iterations" => linear_iterations,
    )
end

function _circle_packing_fixed_point_solve!(
    radii::AbstractVector{Float64},
    boundary,
    ordered_neighbors_idx,
    boundary_mask::AbstractVector{Bool};
    maxiter::Int=200,
    tol::Float64=1.0e-8,
    residual_tol::Float64=max(1.0e-6, sqrt(tol)),
    relaxation::Float64=1.0,
    min_radius::Float64=1.0e-8,
)
    current = copy(radii)
    candidate = similar(current)
    positions = zeros(Float64, length(current), 2)
    A, interior_vertices, interior_index, neighbor_slots, diag_slots, max_degree = _conductance_workspace(ordered_neighbors_idx, boundary_mask)
    rhs = zeros(Float64, size(A, 1), 2)
    conductance_scratch = zeros(Float64, max(max_degree, 1))
    solver_state = nothing
    converged = false
    iter_count = 0
    max_radius_delta = Inf
    residual_inf_norm = Inf

    for iter in 1:maxiter
        _update_boundary_centers!(positions, current, boundary)
        solver_state = _solve_circle_packing_positions!(
            positions,
            rhs,
            solver_state,
            A,
            current,
            ordered_neighbors_idx,
            boundary_mask,
            interior_vertices,
            interior_index,
            neighbor_slots,
            diag_slots,
            conductance_scratch;
            first_solve=(iter == 1),
        )

        copyto!(candidate, current)
        _update_circle_radii!(candidate, ordered_neighbors_idx, boundary_mask, positions)
        _normalize_circle_radii!(candidate, boundary; min_radius=min_radius)

        if relaxation != 1.0
            candidate .= exp.((1 - relaxation) .* log.(max.(current, min_radius)) .+ relaxation .* log.(max.(candidate, min_radius)))
            _normalize_circle_radii!(candidate, boundary; min_radius=min_radius)
        end

        max_radius_delta = maximum(abs.(candidate .- current) ./ max.(current, min_radius))
        residual_inf_norm = _circle_packing_interior_residual(candidate, ordered_neighbors_idx, boundary_mask)
        current, candidate = candidate, current
        iter_count = iter

        if max_radius_delta <= tol && residual_inf_norm <= residual_tol
            converged = true
            break
        end
    end

    copyto!(radii, current)
    return Dict{String,Any}(
        "converged" => converged,
        "iterations" => iter_count,
        "max_radius_delta" => max_radius_delta,
        "residual_inf_norm" => residual_inf_norm,
        "boundary_mask" => boundary_mask,
        "solver" => "fixed_point",
    )
end

function _circle_packing_warm_start!(
    radii::AbstractVector{Float64},
    boundary,
    ordered_neighbors_idx,
    boundary_mask::AbstractVector{Bool};
    maxiter::Int=0,
    tol::Float64=1.0e-8,
    relaxation::Float64=0.85,
    min_radius::Float64=1.0e-8,
)
    maxiter <= 0 && return Dict{String,Any}(
        "converged" => false,
        "iterations" => 0,
        "max_radius_delta" => Inf,
        "residual_inf_norm" => Inf,
        "boundary_mask" => boundary_mask,
        "solver" => "warm_start_skipped",
    )
    return _circle_packing_fixed_point_solve!(
        radii,
        boundary,
        ordered_neighbors_idx,
        boundary_mask;
        maxiter=maxiter,
        tol=max(1.0e-4, sqrt(tol)),
        residual_tol=max(1.0e-3, tol^(1 / 3)),
        relaxation=relaxation,
        min_radius=min_radius,
    )
end

function _circle_packing_gradient_descent_solve!(
    radii::AbstractVector{Float64},
    boundary,
    ordered_neighbors_idx;
    maxiter::Int=200,
    tol::Float64=1.0e-8,
    min_radius::Float64=1.0e-8,
)
    n = length(radii)
    boundary_mask = _boundary_mask(boundary, n)
    boundary_prev_idx, boundary_next_idx = _boundary_neighbor_indices(boundary, n)
    residual = zeros(Float64, n)
    trial_residual = similar(residual)
    u = log.(max.(radii, min_radius))
    trial_u = similar(u)
    trial_radii = similar(radii)
    J, affected_rows, jacobian_slots, _, _ = _circle_packing_jacobian_workspace(ordered_neighbors_idx)
    best_radii = copy(radii)
    best_norm = Inf
    converged = false
    iter_count = 0
    max_radius_delta = Inf

    for iter in 1:maxiter
        _circle_packing_residuals!(residual, radii, ordered_neighbors_idx, boundary_mask, boundary_prev_idx, boundary_next_idx)
        residual_norm = maximum(abs.(residual))
        iter_count = iter
        if residual_norm < best_norm
            best_norm = residual_norm
            copyto!(best_radii, radii)
        end
        if !(isfinite(residual_norm))
            break
        elseif residual_norm <= tol
            converged = true
            max_radius_delta = 0.0
            break
        end

        _fill_circle_packing_jacobian!(
            J,
            radii,
            residual,
            ordered_neighbors_idx,
            boundary_mask,
            boundary_prev_idx,
            boundary_next_idx,
            affected_rows,
            jacobian_slots,
        )
        grad = J' * residual
        grad_scale = max(maximum(abs.(grad)), 1.0)

        accepted = false
        α = 0.2
        while α >= 1.0e-6
            trial_u .= u .- α .* grad ./ grad_scale
            @inbounds for i in eachindex(trial_radii)
                trial_radii[i] = max(exp(trial_u[i]), min_radius)
            end
            if _boundary_cycle_feasible(trial_radii, boundary)
                _circle_packing_residuals!(trial_residual, trial_radii, ordered_neighbors_idx, boundary_mask, boundary_prev_idx, boundary_next_idx)
                trial_norm = maximum(abs.(trial_residual))
                if isfinite(trial_norm) && trial_norm < residual_norm
                    max_radius_delta = maximum(abs.(trial_radii .- radii) ./ max.(radii, min_radius))
                    radii .= trial_radii
                    u .= trial_u
                    accepted = true
                    break
                end
            end
            α *= 0.5
        end

        accepted || break
    end

    if best_norm < maximum(abs.(residual))
        copyto!(radii, best_radii)
        _circle_packing_residuals!(residual, radii, ordered_neighbors_idx, boundary_mask, boundary_prev_idx, boundary_next_idx)
    end
    return Dict{String,Any}(
        "converged" => converged,
        "iterations" => iter_count,
        "max_radius_delta" => max_radius_delta,
        "residual_inf_norm" => maximum(abs.(residual)),
        "boundary_mask" => boundary_mask,
        "solver" => "gradient_descent",
    )
end

function _solve_circle_packing_positions!(
    positions::AbstractMatrix{Float64},
    rhs::AbstractMatrix{Float64},
    workspace,
    A::SparseMatrixCSC{Float64,Int},
    radii,
    ordered_neighbors_idx,
    boundary_mask,
    interior_vertices,
    interior_index,
    neighbor_slots,
    diag_slots,
    conductance_scratch;
    first_solve::Bool=false,
)
    isempty(interior_vertices) && return workspace
    workspace === nothing && (workspace = _circle_packing_position_cg_workspace(size(A, 1)))

    _update_conductance_matrix!(
        A,
        rhs,
        positions,
        radii,
        ordered_neighbors_idx,
        boundary_mask,
        interior_vertices,
        interior_index,
        neighbor_slots,
        diag_slots,
        conductance_scratch,
    )

    workspace.preconditioner = _refresh_ilu0_preconditioner(workspace.preconditioner, A)
    cg_maxiter = _circle_packing_cg_maxiter(size(A, 1); cap=640)
    _, x_history, _ = _circle_packing_cg_solve!(
        workspace.x,
        A,
        vec(view(rhs, :, 1)),
        workspace.x_statevars;
        preconditioner=workspace.preconditioner,
        reltol=1.0e-6,
        abstol=1.0e-10,
        maxiter=cg_maxiter,
    )
    _, y_history, _ = _circle_packing_cg_solve!(
        workspace.y,
        A,
        vec(view(rhs, :, 2)),
        workspace.y_statevars;
        preconditioner=workspace.preconditioner,
        reltol=1.0e-6,
        abstol=1.0e-10,
        maxiter=cg_maxiter,
    )

    if x_history === nothing || y_history === nothing || !all(isfinite, workspace.x) || !all(isfinite, workspace.y)
        factor = lu(A)
        workspace.x .= factor \ rhs[:, 1]
        workspace.y .= factor \ rhs[:, 2]
    end

    positions[interior_vertices, 1] .= workspace.x
    positions[interior_vertices, 2] .= workspace.y
    return workspace
end

function _inverse_stereographic_positions(pos2::AbstractMatrix{<:Real}; sphere_radius::Float64=1.0, projection_scale::Float64=1.0)
    sphere_radius > 0.0 || throw(ArgumentError("sphere_radius must be positive"))
    projection_scale > 0.0 || throw(ArgumentError("projection_scale must be positive"))
    size(pos2, 2) == 2 || throw(ArgumentError("inverse stereographic projection expects 2D positions"))

    out = Matrix{Float64}(undef, size(pos2, 1), 3)
    for i in 1:size(pos2, 1)
        x = projection_scale * float(pos2[i, 1])
        y = projection_scale * float(pos2[i, 2])
        rho2 = x * x + y * y
        denom = rho2 + 1.0
        out[i, 1] = sphere_radius * (2.0 * x / denom)
        out[i, 2] = sphere_radius * (2.0 * y / denom)
        out[i, 3] = sphere_radius * ((rho2 - 1.0) / denom)
    end
    return out
end

function _sphere_circle_geometry_from_disk(pos2::AbstractMatrix{<:Real}, radii::AbstractVector{<:Real}; sphere_radius::Float64=1.0, projection_scale::Float64=1.0)
    sphere_radius > 0.0 || throw(ArgumentError("sphere_radius must be positive"))
    projection_scale > 0.0 || throw(ArgumentError("projection_scale must be positive"))
    size(pos2, 2) == 2 || throw(ArgumentError("sphere circle lifting expects 2D positions"))
    length(radii) == size(pos2, 1) || throw(ArgumentError("circle radii must align with 2D positions"))

    centers = Matrix{Float64}(undef, size(pos2, 1), 3)
    normals = Matrix{Float64}(undef, size(pos2, 1), 3)
    plane_radii = Vector{Float64}(undef, size(pos2, 1))
    offsets = Vector{Float64}(undef, size(pos2, 1))

    for i in 1:size(pos2, 1)
        x = projection_scale * float(pos2[i, 1])
        y = projection_scale * float(pos2[i, 2])
        r = projection_scale * float(radii[i])
        c = x * x + y * y - r * r

        nx = 2.0 * x
        ny = 2.0 * y
        nz = c - 1.0
        nlen = sqrt(nx * nx + ny * ny + nz * nz)
        nlen > 0.0 || throw(ArgumentError("failed to lift a planar circle to the sphere"))

        h = clamp((1.0 + c) / nlen, -1.0, 1.0)
        nx /= nlen
        ny /= nlen
        nz /= nlen

        normals[i, 1] = nx
        normals[i, 2] = ny
        normals[i, 3] = nz
        offsets[i] = h

        centers[i, 1] = sphere_radius * h * nx
        centers[i, 2] = sphere_radius * h * ny
        centers[i, 3] = sphere_radius * h * nz
        plane_radii[i] = sphere_radius * sqrt(max(0.0, 1.0 - h * h))
    end

    return Dict{String,Any}(
        "centers" => centers,
        "normals" => normals,
        "radii" => plane_radii,
        "offsets" => offsets,
        "sphere_radius" => sphere_radius,
        "projection_scale" => projection_scale,
        "projection" => "inverse_stereographic",
    )
end

function _sphere_circle_center_positions(geometry)
    geometry === nothing && return Matrix{Float64}(undef, 0, 3)
    normals = Float64.(geometry["normals"])
    size(normals, 2) == 3 || throw(ArgumentError("sphere circle geometry normals must have shape (N, 3)"))
    sphere_radius = float(get(geometry, "sphere_radius", 1.0))
    return sphere_radius .* normals
end

function _sphere_circle_geometry_to_disk(geometry)
    geometry === nothing && return Matrix{Float64}(undef, 0, 2), Float64[]
    normals = Float64.(geometry["normals"])
    offsets = Float64.(geometry["offsets"])
    size(normals, 2) == 3 || throw(ArgumentError("sphere circle geometry normals must have shape (N, 3)"))
    size(normals, 1) == length(offsets) || throw(ArgumentError("sphere circle geometry normals and offsets must align"))

    pos2 = Matrix{Float64}(undef, size(normals, 1), 2)
    radii = Vector{Float64}(undef, size(normals, 1))

    for i in 1:size(normals, 1)
        nx = normals[i, 1]
        ny = normals[i, 2]
        nz = normals[i, 3]
        h = offsets[i]
        denom = h - nz
        abs(denom) > 1.0e-12 || throw(ArgumentError("failed to stereographically recover a planar circle from the sphere"))

        x = nx / denom
        y = ny / denom
        c = (h + nz) / denom
        r2 = x * x + y * y - c

        pos2[i, 1] = x
        pos2[i, 2] = y
        radii[i] = sqrt(max(0.0, r2))
    end

    return pos2, radii
end

function _sphere_balance_barycenter(entries)
    isempty(entries) && return zeros(Float64, 3)
    acc = zeros(Float64, 3)
    @inbounds for entry in entries
        acc[1] += float(entry.axis[1])
        acc[2] += float(entry.axis[2])
        acc[3] += float(entry.axis[3])
    end
    acc ./= length(entries)
    return acc
end

function _sphere_balance_parameter_steps(params::AbstractVector{<:Real})
    steps = zeros(Float64, length(params))
    @inbounds for i in eachindex(params)
        steps[i] = max(1.0e-6, 1.0e-4 * max(1.0, abs(float(params[i]))))
    end
    return steps
end

function _sphere_balance_geometry_from_affine_params(
    pos2::AbstractMatrix{<:Real},
    radii::AbstractVector{<:Real},
    params::AbstractVector{<:Real};
    sphere_radius::Float64=1.0,
    projection_scale::Float64=1.0,
)
    length(params) == 3 || throw(ArgumentError("sphere balance parameter vector must have length 3"))
    tx = float(params[1])
    ty = float(params[2])
    log_scale = clamp(float(params[3]), -12.0, 12.0)
    effective_projection_scale = float(projection_scale) * exp(log_scale)
    translated = _translate_disk_circle_centers(pos2, tx, ty)
    geometry = _sphere_circle_geometry_from_disk(
        translated,
        radii;
        sphere_radius=sphere_radius,
        projection_scale=effective_projection_scale,
    )
    geometry["outer_geometry"] = _outer_sphere_circle_geometry_for_balance(
        sphere_radius,
        effective_projection_scale;
        translation=(tx, ty),
    )
    geometry["projection_scale"] = effective_projection_scale
    geometry["sphere_balance_translation"] = Float64[tx, ty]
    return geometry
end

function _sphere_balance_affine_residual(
    pos2::AbstractMatrix{<:Real},
    radii::AbstractVector{<:Real},
    params::AbstractVector{<:Real};
    sphere_radius::Float64=1.0,
    projection_scale::Float64=1.0,
)
    geometry = _sphere_balance_geometry_from_affine_params(
        pos2,
        radii,
        params;
        sphere_radius=sphere_radius,
        projection_scale=projection_scale,
    )
    residual = _sphere_balance_barycenter(_sphere_balance_entries(geometry; include_outer=true))
    return residual, geometry
end

function _translate_disk_circle_centers(pos2::AbstractMatrix{<:Real}, tx::Real, ty::Real)
    out = Float64.(pos2)
    @inbounds for i in 1:size(out, 1)
        out[i, 1] -= float(tx)
        out[i, 2] -= float(ty)
    end
    return out
end

function _outer_sphere_circle_geometry_for_balance(
    sphere_radius::Real,
    projection_scale::Real;
    translation=(0.0, 0.0),
)
    tx = float(first(translation))
    ty = float(last(translation))
    return _sphere_circle_geometry_from_disk(
        reshape(Float64[-tx, -ty], 1, 2),
        Float64[1.0];
        sphere_radius=float(sphere_radius),
        projection_scale=float(projection_scale),
    )
end

function _sphere_circle_geometry_copy(geometry)
    geometry === nothing && return nothing
    out = Dict{String,Any}(
        "centers" => Float64.(geometry["centers"]),
        "normals" => Float64.(geometry["normals"]),
        "radii" => Float64.(geometry["radii"]),
        "offsets" => Float64.(geometry["offsets"]),
        "sphere_radius" => float(get(geometry, "sphere_radius", 1.0)),
        "projection_scale" => float(get(geometry, "projection_scale", 1.0)),
        "projection" => get(geometry, "projection", "inverse_stereographic"),
    )
    if haskey(geometry, "outer_geometry")
        out["outer_geometry"] = _sphere_circle_geometry_copy(geometry["outer_geometry"])
    end
    return out
end

@inline _vec3_dot(a::NTuple{3,Float64}, b::NTuple{3,Float64}) = a[1] * b[1] + a[2] * b[2] + a[3] * b[3]
@inline _vec3_dot(a::AbstractVector{<:Real}, b::AbstractVector{<:Real}) = float(a[1]) * float(b[1]) + float(a[2]) * float(b[2]) + float(a[3]) * float(b[3])

@inline function _vec3_cross(a::NTuple{3,Float64}, b::NTuple{3,Float64})
    return (
        a[2] * b[3] - a[3] * b[2],
        a[3] * b[1] - a[1] * b[3],
        a[1] * b[2] - a[2] * b[1],
    )
end

function _normalize_vec3(v; fallback=(0.0, 0.0, 1.0))
    x = float(v[1])
    y = float(v[2])
    z = float(v[3])
    nrm = sqrt(x * x + y * y + z * z)
    if nrm <= 1.0e-12
        fx = float(fallback[1])
        fy = float(fallback[2])
        fz = float(fallback[3])
        fnrm = sqrt(fx * fx + fy * fy + fz * fz)
        fnrm <= 1.0e-12 && return (0.0, 0.0, 1.0)
        return (fx / fnrm, fy / fnrm, fz / fnrm)
    end
    return (x / nrm, y / nrm, z / nrm)
end

function _rotation_matrix_axis_angle(axis, angle::Real)
    x, y, z = _normalize_vec3(axis)
    θ = float(angle)
    c = cos(θ)
    s = sin(θ)
    t = 1.0 - c
    return Float64[
        t*x*x+c t*x*y-s*z t*x*z+s*y;
        t*y*x+s*z t*y*y+c t*y*z-s*x;
        t*z*x-s*y t*z*y+s*x t*z*z+c
    ]
end

function _rotation_matrix_from_unit_vectors(from, to)
    f = _normalize_vec3(from)
    t = _normalize_vec3(to)
    c = clamp(_vec3_dot(f, t), -1.0, 1.0)
    if c >= 1.0 - 1.0e-12
        return Matrix{Float64}(I, 3, 3)
    elseif c <= -1.0 + 1.0e-12
        axis = abs(f[1]) < 0.9 ? _vec3_cross(f, (1.0, 0.0, 0.0)) : _vec3_cross(f, (0.0, 1.0, 0.0))
        return _rotation_matrix_axis_angle(axis, π)
    end

    v = _vec3_cross(f, t)
    vx, vy, vz = v
    s2 = max(vx * vx + vy * vy + vz * vz, 1.0e-18)
    k = Float64[
        0.0 -vz vy;
        vz 0.0 -vx;
        -vy vx 0.0
    ]
    return Matrix{Float64}(I, 3, 3) + k + ((1.0 - c) / s2) .* (k * k)
end

function _rotate_sphere_circle_geometry(geometry, rotation::AbstractMatrix{<:Real})
    geometry === nothing && return nothing
    size(rotation) == (3, 3) || throw(ArgumentError("sphere rotation matrix must have shape (3, 3)"))

    function rotate_rows(arr)
        out = Matrix{Float64}(undef, size(arr, 1), 3)
        for i in 1:size(arr, 1)
            x = float(arr[i, 1])
            y = float(arr[i, 2])
            z = float(arr[i, 3])
            out[i, 1] = rotation[1, 1] * x + rotation[1, 2] * y + rotation[1, 3] * z
            out[i, 2] = rotation[2, 1] * x + rotation[2, 2] * y + rotation[2, 3] * z
            out[i, 3] = rotation[3, 1] * x + rotation[3, 2] * y + rotation[3, 3] * z
        end
        return out
    end

    out = Dict{String,Any}(
        "centers" => rotate_rows(geometry["centers"]),
        "normals" => rotate_rows(geometry["normals"]),
        "radii" => Float64.(geometry["radii"]),
        "offsets" => Float64.(geometry["offsets"]),
        "sphere_radius" => float(get(geometry, "sphere_radius", 1.0)),
        "projection_scale" => float(get(geometry, "projection_scale", 1.0)),
        "projection" => get(geometry, "projection", "inverse_stereographic"),
    )
    if haskey(geometry, "outer_geometry")
        out["outer_geometry"] = _rotate_sphere_circle_geometry(geometry["outer_geometry"], rotation)
    end
    return out
end

function _sphere_balance_outer_geometry(
    geometry;
    include_outer::Bool=false,
    outer_geometry=nothing,
    outer_translation=(0.0, 0.0),
    outer_projection_scale::Real=get(geometry, "projection_scale", 1.0),
)
    include_outer || return nothing
    if outer_geometry !== nothing
        return _sphere_circle_geometry_copy(outer_geometry)
    elseif haskey(geometry, "outer_geometry")
        return _sphere_circle_geometry_copy(geometry["outer_geometry"])
    end
    return _outer_sphere_circle_geometry_for_balance(
        get(geometry, "sphere_radius", 1.0),
        outer_projection_scale;
        translation=outer_translation,
    )
end

function _sphere_balance_cap_area_fractions(
    geometry;
    include_outer::Bool=false,
    outer_geometry=nothing,
    outer_translation=(0.0, 0.0),
    outer_projection_scale::Real=get(geometry, "projection_scale", 1.0),
)
    geometry === nothing && return Float64[]
    offsets = Float64.(geometry["offsets"])
    fractions = Vector{Float64}(undef, length(offsets) + (include_outer ? 1 : 0))
    @inbounds for i in eachindex(offsets)
        fractions[i] = 0.5 * (1.0 - clamp(offsets[i], -1.0, 1.0))
    end
    if include_outer
        outer_geometry = _sphere_balance_outer_geometry(
            geometry;
            include_outer=true,
            outer_geometry=outer_geometry,
            outer_translation=outer_translation,
            outer_projection_scale=outer_projection_scale,
        )
        outer_offset = outer_geometry === nothing ? 0.0 : float(outer_geometry["offsets"][1])
        fractions[end] = 0.5 * (1.0 + clamp(outer_offset, -1.0, 1.0))
    end
    return fractions
end

function _sphere_balance_summary(
    geometry;
    include_outer::Bool=false,
    outer_geometry=nothing,
    outer_translation=(0.0, 0.0),
    outer_projection_scale::Real=get(geometry, "projection_scale", 1.0),
)
    fractions = _sphere_balance_cap_area_fractions(
        geometry;
        include_outer=include_outer,
        outer_geometry=outer_geometry,
        outer_translation=outer_translation,
        outer_projection_scale=outer_projection_scale,
    )
    isempty(fractions) && return (objective=0.0, max_fraction=0.0, min_fraction=0.0, outer_fraction=0.0)

    sorted_desc = sort(fractions; rev=true)
    top_count = min(length(sorted_desc), 6)
    top_mean = sum(@view(sorted_desc[1:top_count])) / top_count
    safe = max.(fractions, eps(Float64))
    logs = log.(safe)
    log_center = sum(logs) / length(logs)
    log_centered = logs .- log_center
    log_rms = sqrt(sum(abs2, log_centered) / length(log_centered))
    log_span = maximum(logs) - minimum(logs)
    objective = 2.4 * sorted_desc[1] + 0.9 * top_mean + 1.1 * log_span + 0.45 * log_rms
    outer_fraction = include_outer ? fractions[end] : 0.0
    return (
        objective=objective,
        max_fraction=sorted_desc[1],
        min_fraction=sorted_desc[end],
        outer_fraction=outer_fraction,
    )
end

function _sphere_balance_entries(
    geometry;
    include_outer::Bool=false,
    outer_geometry=nothing,
    outer_translation=(0.0, 0.0),
    outer_projection_scale::Real=get(geometry, "projection_scale", 1.0),
)
    geometry === nothing && return NamedTuple[]
    normals = Float64.(geometry["normals"])
    offsets = Float64.(geometry["offsets"])
    entries = NamedTuple[]
    for i in 1:size(normals, 1)
        offset = clamp(offsets[i], -1.0, 1.0)
        push!(entries, (
            fraction=0.5 * (1.0 - offset),
            axis=_normalize_vec3(view(normals, i, :)),
            is_outer=false,
            index=i,
        ))
    end
    if include_outer
        outer = _sphere_balance_outer_geometry(
            geometry;
            include_outer=true,
            outer_geometry=outer_geometry,
            outer_translation=outer_translation,
            outer_projection_scale=outer_projection_scale,
        )
        if outer !== nothing
            n = _normalize_vec3(view(outer["normals"], 1, :))
            offset = clamp(float(outer["offsets"][1]), -1.0, 1.0)
            push!(entries, (
                fraction=0.5 * (1.0 + offset),
                axis=(-n[1], -n[2], -n[3]),
                is_outer=true,
                index=1,
            ))
        end
    end
    return entries
end

function _apply_sphere_mobius_balance(
    pos2::AbstractMatrix{<:Real},
    radii::AbstractVector{<:Real};
    sphere_radius::Float64=1.0,
    projection_scale::Float64=1.0,
    max_rounds::Integer=50,
)
    params = zeros(Float64, 3)
    tol = 1.0e-10
    residual, geometry = _sphere_balance_affine_residual(
        pos2,
        radii,
        params;
        sphere_radius=sphere_radius,
        projection_scale=projection_scale,
    )
    residual_norm = norm(residual)
    initial_norm = residual_norm
    best_norm = residual_norm
    best_params = copy(params)
    best_geometry = geometry
    converged = residual_norm <= tol
    accepted_steps = 0
    iter_count = 0

    for iter in 1:Int(max_rounds)
        iter_count = iter
        residual_norm <= tol && (converged = true; break)

        steps = _sphere_balance_parameter_steps(params)
        J = Matrix{Float64}(undef, 3, 3)
        jacobian_ok = true
        for j in 1:3
            δ = steps[j]
            trial_plus = copy(params)
            trial_minus = copy(params)
            trial_plus[j] += δ
            trial_minus[j] -= δ
            try
                plus_residual, _ = _sphere_balance_affine_residual(
                    pos2,
                    radii,
                    trial_plus;
                    sphere_radius=sphere_radius,
                    projection_scale=projection_scale,
                )
                minus_residual, _ = _sphere_balance_affine_residual(
                    pos2,
                    radii,
                    trial_minus;
                    sphere_radius=sphere_radius,
                    projection_scale=projection_scale,
                )
                J[:, j] .= (plus_residual .- minus_residual) ./ (2.0 * δ)
            catch
                jacobian_ok = false
                break
            end
        end
        jacobian_ok || break

        g = transpose(J) * residual
        step = nothing
        for λ in (1.0e-10, 1.0e-8, 1.0e-6, 1.0e-4, 1.0e-2, 1.0, 100.0)
            try
                candidate_step = -((transpose(J) * J + λ * Matrix{Float64}(I, 3, 3)) \ g)
                if all(isfinite, candidate_step)
                    step = candidate_step
                    break
                end
            catch
            end
        end
        if step === nothing
            ng = norm(g)
            step = ng <= 1.0e-14 ? zeros(Float64, 3) : -(g / ng)
        end
        all(isfinite, step) || break

        accepted = false
        α = 1.0
        while α >= 1.0e-4
            trial_params = params .+ α .* step
            try
                trial_residual, trial_geometry = _sphere_balance_affine_residual(
                    pos2,
                    radii,
                    trial_params;
                    sphere_radius=sphere_radius,
                    projection_scale=projection_scale,
                )
                trial_norm = norm(trial_residual)
                if isfinite(trial_norm) && trial_norm + 1.0e-12 < residual_norm
                    params .= trial_params
                    residual .= trial_residual
                    residual_norm = trial_norm
                    geometry = trial_geometry
                    accepted = true
                    accepted_steps += 1
                    if trial_norm < best_norm
                        best_norm = trial_norm
                        best_params .= trial_params
                        best_geometry = trial_geometry
                    end
                    break
                end
            catch
            end
            α *= 0.5
        end

        accepted || break
    end

    if best_norm + 1.0e-14 < residual_norm
        params .= best_params
        geometry = best_geometry
        residual, _ = _sphere_balance_affine_residual(
            pos2,
            radii,
            params;
            sphere_radius=sphere_radius,
            projection_scale=projection_scale,
        )
        residual_norm = norm(residual)
    end

    converged = residual_norm <= tol
    effective_projection_scale = float(get(geometry, "projection_scale", projection_scale))
    translation = Float64.(get(geometry, "sphere_balance_translation", [0.0, 0.0]))

    return geometry, Dict{String,Any}(
        "sphere_balance_applied" => accepted_steps > 0,
        "sphere_balance_method" => "affine_barycenter_newton",
        "sphere_balance_iterations" => iter_count,
        "sphere_balance_refine_rounds" => accepted_steps,
        "sphere_balance_max_rounds" => Int(max_rounds),
        "sphere_balance_converged" => converged,
        "sphere_balance_barycenter_norm_before" => initial_norm,
        "sphere_balance_barycenter_norm_after" => residual_norm,
        "sphere_balance_translation" => translation,
        "sphere_balance_projection_scale_before" => float(projection_scale),
        "sphere_balance_projection_scale_after" => effective_projection_scale,
        "sphere_balance_projection_scale_factor" => effective_projection_scale / float(projection_scale),
    )
end

function _rotate_sphere_circle_geometry_z(geometry, angle::Real)
    geometry === nothing && return nothing
    θ = float(angle)
    abs(θ) <= 1.0e-12 && return geometry
    c = cos(θ)
    s = sin(θ)
    rotation = Float64[
        c -s 0.0;
        s c 0.0;
        0.0 0.0 1.0
    ]
    return _rotate_sphere_circle_geometry(geometry, rotation)
end

function _sphere_circle_normalization_angle(geometry, boundary_vertices)
    geometry === nothing && return 0.0
    boundary = Int.(collect(boundary_vertices))
    isempty(boundary) && return 0.0
    centers = _sphere_circle_center_positions(geometry)
    idx = boundary[1] + 1
    1 <= idx <= size(centers, 1) || return 0.0
    x = centers[idx, 1]
    y = centers[idx, 2]
    abs(x) + abs(y) <= 1.0e-12 && return 0.0
    return -atan(y, x)
end

function _normalize_sphere_circle_geometry(geometry, boundary_vertices)
    geometry === nothing && return nothing
    angle = _sphere_circle_normalization_angle(geometry, boundary_vertices)
    abs(angle) <= 1.0e-12 && return geometry
    return _rotate_sphere_circle_geometry_z(geometry, angle)
end

function compute_circle_packing_layout(
    num_vertices::Integer,
    edges,
    boundary_vertices;
    faces=nothing,
    triangles=nothing,
    triangle_edge_ids=nothing,
    maxiter::Int=200,
    tol::Float64=1.0e-8,
    relaxation::Float64=1.0,
    initial_radius::Float64=0.5,
    min_radius::Float64=1.0e-8,
    project_to_sphere::Bool=false,
    sphere_radius::Float64=1.0,
    sphere_projection_scale::Float64=1.0,
    return_metadata::Bool=false,
)
    n = Int(num_vertices)
    n >= 0 || throw(ArgumentError("num_vertices must be nonnegative"))
    relaxation > 0.0 || throw(ArgumentError("relaxation must be positive"))
    initial_radius > 0.0 || throw(ArgumentError("initial_radius must be positive"))
    min_radius > 0.0 || throw(ArgumentError("min_radius must be positive"))
    sphere_radius > 0.0 || throw(ArgumentError("sphere_radius must be positive"))
    sphere_projection_scale > 0.0 || throw(ArgumentError("sphere_projection_scale must be positive"))

    boundary = _circle_packing_boundary(boundary_vertices, n)
    tri, triangulation_source = _circle_packing_triangles(n; faces=faces, triangles=triangles)
    tri, ordered_neighbors, neighbor_order_source = _validate_disk_triangulation(n, tri, boundary; triangle_edge_ids=triangle_edge_ids)

    original_edges = sanitize_edge_array(edges)
    size(original_edges, 1) == 0 && throw(ArgumentError("circle packing requires a nonempty edge set"))
    minimum(original_edges) >= 0 || throw(ArgumentError("edge endpoints must be >= 0"))
    maximum(original_edges) < n || throw(ArgumentError("edge endpoints must lie in [0, num_vertices)"))

    packing_edges = first(collapse_undirected_edges(vcat(original_edges, _triangle_edge_array(tri), _boundary_cycle_edges(boundary)); drop_loops=true))

    ordered_neighbors_idx = _neighbor_index_lists(ordered_neighbors)
    radii, boundary_mask = _initialize_circle_packing_radii(
        ordered_neighbors_idx,
        boundary;
        initial_radius=initial_radius,
        min_radius=min_radius,
    )
    warmup_iters = n >= 512 ? min(maxiter ÷ 4, 64) : 0
    warmup_meta = _circle_packing_warm_start!(
        radii,
        boundary,
        ordered_neighbors_idx,
        boundary_mask;
        maxiter=warmup_iters,
        tol=tol,
        relaxation=min(relaxation, 0.85),
        min_radius=min_radius,
    )
    newton_meta = _circle_packing_newton_solve!(
        radii,
        boundary,
        ordered_neighbors_idx;
        maxiter=maxiter,
        tol=tol,
        min_radius=min_radius,
    )
    converged = Bool(get(newton_meta, "converged", false))
    iter_count = Int(get(warmup_meta, "iterations", 0)) + Int(get(newton_meta, "iterations", 0))
    max_radius_delta = Float64(get(newton_meta, "max_radius_delta", Inf))
    boundary_mask = get(newton_meta, "boundary_mask", boundary_mask)
    radius_solver = "newton_log_radii"
    residual_inf_norm = Float64(get(newton_meta, "residual_inf_norm", Inf))

    if !converged || !_boundary_cycle_feasible(radii, boundary)
        fallback_seed, boundary_mask = _initialize_circle_packing_radii(
            ordered_neighbors_idx,
            boundary;
            initial_radius=initial_radius,
            min_radius=min_radius,
        )
        gradient_radii = copy(fallback_seed)
        gradient_meta = _circle_packing_gradient_descent_solve!(
            gradient_radii,
            boundary,
            ordered_neighbors_idx;
            maxiter=maxiter,
            tol=tol,
            min_radius=min_radius,
        )
        fallback_radii = copy(fallback_seed)
        fallback_meta = _circle_packing_fixed_point_solve!(
            fallback_radii,
            boundary,
            ordered_neighbors_idx,
            boundary_mask;
            maxiter=maxiter,
            tol=tol,
            relaxation=relaxation,
            min_radius=min_radius,
        )
        gradient_norm = Float64(get(gradient_meta, "residual_inf_norm", Inf))
        fallback_norm = Float64(get(fallback_meta, "residual_inf_norm", Inf))
        if gradient_norm <= fallback_norm
            radii .= gradient_radii
            converged = Bool(get(gradient_meta, "converged", false))
            iter_count = Int(get(gradient_meta, "iterations", 0))
            max_radius_delta = Float64(get(gradient_meta, "max_radius_delta", max_radius_delta))
            boundary_mask = get(gradient_meta, "boundary_mask", boundary_mask)
            residual_inf_norm = gradient_norm
            newton_meta = merge(newton_meta, Dict("fallback_solver" => "gradient_descent", "warm_start_solver" => get(warmup_meta, "solver", "warm_start")))
            radius_solver = "gradient_fallback"
        else
            radii .= fallback_radii
            converged = Bool(get(fallback_meta, "converged", false))
            iter_count = Int(get(fallback_meta, "iterations", 0))
            max_radius_delta = Float64(get(fallback_meta, "max_radius_delta", max_radius_delta))
            boundary_mask = get(fallback_meta, "boundary_mask", boundary_mask)
            residual_inf_norm = fallback_norm
            newton_meta = merge(newton_meta, Dict("fallback_solver" => "fixed_point", "warm_start_solver" => get(warmup_meta, "solver", "warm_start")))
            radius_solver = "fixed_point_fallback"
        end
    end

    positions = zeros(Float64, n, 2)
    rhs = similar(positions)
    A, interior_vertices, interior_index, neighbor_slots, diag_slots, max_degree = _conductance_workspace(ordered_neighbors_idx, boundary_mask)
    conductance_scratch = zeros(Float64, max(max_degree, 1))
    rhs = zeros(Float64, size(A, 1), 2)
    _update_boundary_centers!(positions, radii, boundary)
    _solve_circle_packing_positions!(
        positions,
        rhs,
        nothing,
        A,
        radii,
        ordered_neighbors_idx,
        boundary_mask,
        interior_vertices,
        interior_index,
        neighbor_slots,
        diag_slots,
        conductance_scratch;
        first_solve=true,
    )

    metadata = Dict{String,Any}(
        "boundary_positions_mode" => project_to_sphere ? "sphere_inverse_stereographic" : "unit_outer_circle",
        "triangulation_source" => triangulation_source,
        "neighbor_order_source" => neighbor_order_source,
        "packing_num_edges" => size(packing_edges, 1),
        "packing_num_triangles" => size(tri, 1),
        "triangle_multiplicity_preserved" => true,
        "iterations" => iter_count,
        "converged" => converged,
        "max_radius_delta" => max_radius_delta,
        "radius_solver" => radius_solver,
        "radius_residual_inf_norm" => residual_inf_norm,
        "warm_start_iterations" => Int(get(warmup_meta, "iterations", 0)),
        "warm_start_solver" => get(warmup_meta, "solver", "warm_start"),
        "boundary_angle_residual" => abs(_boundary_angle_sum_error(1.0, boundary, radii)),
        "interior_angle_residual" => _circle_packing_interior_residual(radii, ordered_neighbors_idx, boundary_mask),
        "edge_tangency_residual" => _circle_packing_edge_residual(positions, radii, packing_edges),
        "boundary_vertices" => collect(boundary),
    )

    out_positions = positions
    out_radii = radii
    if project_to_sphere
        requested_projection_scale = sphere_projection_scale
        initial_geometry = _sphere_circle_geometry_from_disk(
            positions,
            radii;
            sphere_radius=sphere_radius,
            projection_scale=requested_projection_scale,
        )
        initial_geometry["outer_geometry"] = _outer_sphere_circle_geometry_for_balance(
            sphere_radius,
            requested_projection_scale,
        )
        initial_summary = _sphere_balance_summary(
            initial_geometry;
            include_outer=true,
            outer_geometry=initial_geometry["outer_geometry"],
            outer_projection_scale=requested_projection_scale,
        )
        sphere_geometry, balance_meta = _apply_sphere_mobius_balance(
            positions,
            radii;
            sphere_radius=sphere_radius,
            projection_scale=requested_projection_scale,
        )
        haskey(sphere_geometry, "outer_geometry") || (sphere_geometry["outer_geometry"] = _outer_sphere_circle_geometry_for_balance(sphere_radius, requested_projection_scale))
        effective_projection_scale = float(get(sphere_geometry, "projection_scale", requested_projection_scale))
        normalization_angle = _sphere_circle_normalization_angle(sphere_geometry, boundary)
        sphere_geometry = _rotate_sphere_circle_geometry_z(sphere_geometry, normalization_angle)
        outer_geometry = haskey(sphere_geometry, "outer_geometry") ? sphere_geometry["outer_geometry"] : nothing
        final_summary = _sphere_balance_summary(
            sphere_geometry;
            include_outer=true,
            outer_geometry=outer_geometry,
            outer_projection_scale=effective_projection_scale,
        )
        metadata["sphere_projection"] = "inverse_stereographic"
        metadata["packing_topology"] = "sphere"
        metadata["sphere_radius"] = sphere_radius
        metadata["sphere_projection_scale"] = sphere_projection_scale
        metadata["sphere_requested_projection_scale"] = requested_projection_scale
        metadata["sphere_effective_projection_scale"] = effective_projection_scale
        metadata["sphere_normalization_angle"] = normalization_angle
        merge!(metadata, balance_meta)
        metadata["sphere_balance_objective_before"] = initial_summary.objective
        metadata["sphere_balance_objective_after"] = final_summary.objective
        metadata["sphere_balance_max_cap_fraction_before"] = initial_summary.max_fraction
        metadata["sphere_balance_max_cap_fraction_after"] = final_summary.max_fraction
        metadata["sphere_balance_min_cap_fraction_before"] = initial_summary.min_fraction
        metadata["sphere_balance_min_cap_fraction_after"] = final_summary.min_fraction
        metadata["sphere_balance_outer_cap_fraction_before"] = initial_summary.outer_fraction
        metadata["sphere_balance_outer_cap_fraction_after"] = final_summary.outer_fraction
        metadata["sphere_position_mode"] = "circle_centers"
        metadata["sphere_outer_circle_geometry"] = outer_geometry
        metadata["sphere_circle_geometry"] = sphere_geometry
        out_positions = _sphere_circle_center_positions(sphere_geometry)
        out_radii = sphere_geometry["radii"]
    else
        metadata["packing_topology"] = "disk"
    end

    if return_metadata
        return out_positions, out_radii, metadata
    end
    return out_positions, out_radii
end

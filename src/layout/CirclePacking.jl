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
        raw[row + 1, 1] = b
        raw[row + 1, 2] = c
        raw[row + 2, 1] = c
        raw[row + 2, 2] = a
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
        total += _boundary_angle(radii[boundary[i] + 1], radii[boundary[j] + 1], outer_radius)
    end
    return total - 2π
end

function _find_outer_radius(boundary_vertices, radii::AbstractVector{<:Real}; tol::Float64=1.0e-12, maxiter::Int=96)
    boundary = Int.(collect(boundary_vertices))
    lower = 0.0
    for i in eachindex(boundary)
        j = i == length(boundary) ? 1 : i + 1
        lower = max(lower, float(radii[boundary[i] + 1]) + float(radii[boundary[j] + 1]))
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
    radii .= max.(Float64.(radii), min_radius)
    outer = _find_outer_radius(boundary_vertices, radii)
    radii ./= outer
    return radii
end

function _ordered_neighbors_from_triangles(num_vertices::Integer, triangles, boundary_vertices)
    n = Int(num_vertices)
    tri = sanitize_triangles(triangles; drop_degenerate=true, deduplicate=false)
    boundary = Int.(collect(boundary_vertices))
    boundary_mask = falses(n)
    boundary_mask[boundary .+ 1] .= true
    boundary_pos = Dict(v => i for (i, v) in enumerate(boundary))

    wedges = [Vector{NTuple{2,Int32}}() for _ in 1:n]

    for i in 1:size(tri, 1)
        a = tri[i, 1]
        b = tri[i, 2]
        c = tri[i, 3]
        push!(wedges[Int(a) + 1], (b, c))
        push!(wedges[Int(b) + 1], (c, a))
        push!(wedges[Int(c) + 1], (a, b))
    end

    orders = Vector{Vector{Int32}}(undef, n)
    for vertex in 0:(n - 1)
        local_wedges = wedges[vertex + 1]
        if isempty(local_wedges)
            orders[vertex + 1] = Int32[]
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

        if boundary_mask[vertex + 1]
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
            orders[vertex + 1] = order
        else
            isempty(odd_vertices) || throw(ArgumentError("circle packing expects an interior disk vertex to induce an Eulerian cycle, but vertex $vertex has odd endpoints $(Tuple(Int.(odd_vertices)))"))
            start = minimum(collect(keys(local_adj)))
            cycle = eulerian_trail(start)
            length(cycle) >= 2 && cycle[1] == cycle[end] || throw(ArgumentError("failed to close an interior neighbor cycle"))
            order = cycle[1:end-1]
            length(order) == length(local_wedges) || throw(ArgumentError("interior neighbor order has the wrong length"))
            orders[vertex + 1] = order
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
    boundary_mask[boundary .+ 1] .= true
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
    for vertex in 0:(n - 1)
        local_adj_v = local_adj[vertex + 1]
        if isempty(local_adj_v)
            orders[vertex + 1] = Int32[]
            continue
        end

        for values in values(local_adj_v)
            sort!(values; by=x -> (x[1], x[2]))
        end

        odd_edges = sort!(Int32[eid for (eid, d) in local_deg[vertex + 1] if isodd(d)])

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

        neighbor_of = local_neighbor[vertex + 1]
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

        if boundary_mask[vertex + 1]
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
                length(edge_path) == local_wedge_count[vertex + 1] + 1 || throw(ArgumentError("boundary neighbor order has the wrong length"))
                orders[vertex + 1] = Int32[neighbor_of[eid] for eid in edge_path]
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
                length(stitched) == length(local_deg[vertex + 1]) || throw(ArgumentError("boundary neighbor order has the wrong length"))
                orders[vertex + 1] = stitched
            end
        else
            if component_count == 1
                if isempty(odd_edges)
                    start_edge = minimum(collect(keys(local_adj_v)))
                    edge_cycle = eulerian_trail(local_adj_v, start_edge)
                    length(edge_cycle) >= 2 && edge_cycle[1] == edge_cycle[end] || throw(ArgumentError("failed to close an interior neighbor cycle"))
                    order_edges = edge_cycle[1:end-1]
                    length(order_edges) == local_wedge_count[vertex + 1] || throw(ArgumentError("interior neighbor order has the wrong length"))
                    orders[vertex + 1] = Int32[neighbor_of[eid] for eid in order_edges]
                elseif length(odd_edges) == 2 && neighbor_of[odd_edges[1]] == neighbor_of[odd_edges[2]]
                    edge_path = eulerian_trail(local_adj_v, odd_edges[1])
                    edge_path[end] == odd_edges[2] || throw(ArgumentError("failed to recover an interior repeated-neighbor cycle"))
                    length(edge_path) == local_wedge_count[vertex + 1] + 1 || throw(ArgumentError("interior neighbor order has the wrong length"))
                    orders[vertex + 1] = Int32[neighbor_of[eid] for eid in edge_path]
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
                length(stitched) == length(local_deg[vertex + 1]) || throw(ArgumentError("interior neighbor order has the wrong length"))
                orders[vertex + 1] = stitched
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
        prev_vertex = boundary[i - 1] + 1
        angle += _boundary_angle(radii[vertex], radii[prev_vertex], 1.0)
        radius_to_center = 1.0 - radii[vertex]
        positions[vertex, 1] = radius_to_center * cos(angle)
        positions[vertex, 2] = radius_to_center * sin(angle)
    end
    return positions
end

function _edge_conductance(r_center::Real, r_current::Real, r_prev::Real, r_next::Real)
    value = (
        sqrt(max(0.0, float(r_center) * float(r_prev) * float(r_current) / (float(r_center) + float(r_prev) + float(r_current)))) +
        sqrt(max(0.0, float(r_center) * float(r_next) * float(r_current) / (float(r_center) + float(r_next) + float(r_current))))
    ) / (float(r_center) + float(r_current))
    return isfinite(value) ? value : 0.0
end

function _conductance_matrix(radii::AbstractVector{<:Real}, ordered_neighbors, boundary_vertices)
    n = length(radii)
    boundary_mask = falses(n)
    boundary_mask[Int.(collect(boundary_vertices)) .+ 1] .= true

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
                radii[Int(nbd[j]) + 1],
                radii[Int(nbd[prev]) + 1],
                radii[Int(nbd[next]) + 1],
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
    boundary_mask[Int.(collect(boundary_vertices)) .+ 1] .= true

    for i in 1:n
        nbd = ordered_neighbors[i]
        isempty(nbd) && continue

        if boundary_mask[i]
            length(nbd) >= 2 || continue
            phi_total = 0.0
            radius_acc = 0.0
            for j in 1:(length(nbd) - 1)
                u = Int(nbd[j]) + 1
                v = Int(nbd[j + 1]) + 1
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

function _circle_packing_interior_residual(radii::AbstractVector{<:Real}, ordered_neighbors, boundary_vertices)
    boundary_mask = falses(length(radii))
    boundary_mask[Int.(collect(boundary_vertices)) .+ 1] .= true
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
                    1 / ((1 + radii[i] / radii[Int(nbd[prev]) + 1]) * (1 + radii[i] / radii[Int(nbd[j]) + 1])),
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
    return_metadata::Bool=false,
)
    n = Int(num_vertices)
    n >= 0 || throw(ArgumentError("num_vertices must be nonnegative"))
    relaxation > 0.0 || throw(ArgumentError("relaxation must be positive"))
    initial_radius > 0.0 || throw(ArgumentError("initial_radius must be positive"))
    min_radius > 0.0 || throw(ArgumentError("min_radius must be positive"))

    boundary = _circle_packing_boundary(boundary_vertices, n)
    tri, triangulation_source = _circle_packing_triangles(n; faces=faces, triangles=triangles)
    tri, ordered_neighbors, neighbor_order_source = _validate_disk_triangulation(n, tri, boundary; triangle_edge_ids=triangle_edge_ids)

    original_edges = sanitize_edge_array(edges)
    size(original_edges, 1) == 0 && throw(ArgumentError("circle packing requires a nonempty edge set"))
    minimum(original_edges) >= 0 || throw(ArgumentError("edge endpoints must be >= 0"))
    maximum(original_edges) < n || throw(ArgumentError("edge endpoints must lie in [0, num_vertices)"))

    packing_edges = first(collapse_undirected_edges(vcat(original_edges, _triangle_edge_array(tri), _boundary_cycle_edges(boundary)); drop_loops=true))

    radii = fill(float(initial_radius), n)
    _normalize_circle_radii!(radii, boundary; min_radius=min_radius)

    positions = zeros(Float64, n, 2)
    converged = false
    iter_count = 0
    max_radius_delta = Inf

    for iter in 1:maxiter
        _update_boundary_centers!(positions, radii, boundary)
        positions = Matrix{Float64}(_conductance_matrix(radii, ordered_neighbors, boundary) \ positions)

        candidate = copy(radii)
        _update_circle_radii!(candidate, ordered_neighbors, boundary, positions)
        _normalize_circle_radii!(candidate, boundary; min_radius=min_radius)

        if relaxation != 1.0
            candidate .= exp.((1 - relaxation) .* log.(max.(radii, min_radius)) .+ relaxation .* log.(max.(candidate, min_radius)))
            _normalize_circle_radii!(candidate, boundary; min_radius=min_radius)
        end

        max_radius_delta = maximum(abs.(candidate .- radii) ./ max.(radii, min_radius))
        radii = candidate
        iter_count = iter

        if max_radius_delta <= tol
            converged = true
            break
        end
    end

    _update_boundary_centers!(positions, radii, boundary)
    positions = Matrix{Float64}(_conductance_matrix(radii, ordered_neighbors, boundary) \ positions)

    metadata = Dict{String,Any}(
        "boundary_positions_mode" => "unit_outer_circle",
        "triangulation_source" => triangulation_source,
        "neighbor_order_source" => neighbor_order_source,
        "packing_num_edges" => size(packing_edges, 1),
        "packing_num_triangles" => size(tri, 1),
        "triangle_multiplicity_preserved" => true,
        "iterations" => iter_count,
        "converged" => converged,
        "max_radius_delta" => max_radius_delta,
        "boundary_angle_residual" => abs(_boundary_angle_sum_error(1.0, boundary, radii)),
        "interior_angle_residual" => _circle_packing_interior_residual(radii, ordered_neighbors, boundary),
        "edge_tangency_residual" => _circle_packing_edge_residual(positions, radii, packing_edges),
        "boundary_vertices" => collect(boundary),
    )

    if return_metadata
        return positions, radii, metadata
    end
    return positions, radii
end

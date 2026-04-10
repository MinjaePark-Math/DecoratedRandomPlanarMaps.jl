# Schnyder-type triangulation model.

const RIGHT_STEP = 0
const LEFT_STEP = 1
const UP_STEP = 2
const DOWN_STEP = 3

struct SchnyderMap <: AbstractRandomPlanarMap
    m::PlanarMap
    faces::Vector{Vector{Int32}}
    edge_u::Vector{Int32}
    edge_v::Vector{Int32}
    edge_color::Vector{String}
    nverts::Int
end

const SchnyderWoodMap = SchnyderMap

num_vertices(m::SchnyderMap) = m.nverts
num_faces(m::SchnyderMap) = length(m.faces)

function layout_edges(m::SchnyderMap; drop_loops::Bool=true)
    pairs = hcat(m.edge_u, m.edge_v)
    return first(collapse_undirected_edges(pairs; drop_loops=drop_loops))
end

random_dyck_path(halflength::Integer, rng::AbstractRNG) = sample_uniform_dyck_path(halflength, rng)

function random_subset_indicator(total::Integer, chosen::Integer, rng::AbstractRNG)
    inds = falses(Int(total))
    if chosen > 0
        perm = randperm(rng, Int(total))
        for idx in perm[1:Int(chosen)]
            inds[idx] = true
        end
    end
    return inds
end

function random_walk_in_quarter_plane(halflength::Integer, rng::AbstractRNG)
    h = Int(halflength)
    k = round(Int, randn(rng) * sqrt(h / 8.0) + 0.5 * h)
    k = clamp(k, 0, h)

    steps1 = random_dyck_path(k, rng)
    steps2 = random_dyck_path(h - k, rng)
    indicator = random_subset_indicator(2h, 2k, rng)

    num1 = 1
    num2 = 1
    steps = Int[]
    sizehint!(steps, 2h)

    for i in 1:(2h)
        if indicator[i]
            push!(steps, steps1[num1] == 1 ? UP_STEP : DOWN_STEP)
            num1 += 1
        else
            push!(steps, steps2[num2] == 1 ? RIGHT_STEP : LEFT_STEP)
            num2 += 1
        end
    end
    return steps
end

function increment_by_step(step::Integer, x::Integer, y::Integer)
    if step == UP_STEP
        return Int(x), Int(y + 1)
    elseif step == DOWN_STEP
        return Int(x), Int(y - 1)
    elseif step == RIGHT_STEP
        return Int(x + 1), Int(y)
    else
        return Int(x - 1), Int(y)
    end
end

mirror_step(step::Integer) = step == RIGHT_STEP ? UP_STEP : step == LEFT_STEP ? DOWN_STEP : step == UP_STEP ? RIGHT_STEP : LEFT_STEP

function reflect_walk(walk)
    x = 0
    y = 0
    out = Int[]
    for step in walk
        above_diag = y > x
        x, y = increment_by_step(step, x, y)
        if y > x || above_diag
            push!(out, mirror_step(step))
        else
            push!(out, Int(step))
        end
    end
    return out
end

function walk_to_path(walk)
    x = 0
    y = 0
    path = [[x, y]]
    for step in walk
        x, y = increment_by_step(step, x, y)
        push!(path, [x, y])
    end
    return path
end

function path_to_walk(path)
    walk = Int[]
    for i in 2:length(path)
        x0, y0 = path[i - 1]
        x1, y1 = path[i]
        if x1 == x0
            push!(walk, y1 > y0 ? UP_STEP : DOWN_STEP)
        else
            push!(walk, x1 > x0 ? RIGHT_STEP : LEFT_STEP)
        end
    end
    return walk
end

function try_mirror_path!(path, pt0::Integer, pt1::Integer)
    a = Int(pt0)
    b = Int(pt1)
    for i in (a + 1):(b - 1)
        if path[a][1] - path[i][2] + path[a][2] < 0
            return
        end
        if path[a][2] - path[i][1] + path[a][1] < 0
            return
        end
    end

    for i in (a + 1):(b - 1)
        tmp = path[a][1] - path[i][2] + path[a][2]
        path[i][2] = path[a][2] - path[i][1] + path[a][1]
        path[i][1] = tmp
    end
    return
end

function markov_update_path!(path, rng::AbstractRNG)
    length(path) < 4 && return

    large_move_frequency = 10
    small_move_size = length(path) < 40 ? fld(length(path), 2) : 20
    small_move_size < 2 && return

    if rand(rng, 0:(large_move_frequency - 1)) == 0
        pt0 = rand(rng, 1:length(path))
        pt1 = rand(rng, 1:length(path))
        if pt1 < pt0
            pt0, pt1 = pt1, pt0
        end
        pt1 <= pt0 + 1 && return
    else
        move_size = rand(rng, 2:small_move_size)
        pt0 = rand(rng, 1:(length(path) - move_size))
        pt1 = pt0 + move_size
    end

    if path[pt1][1] + path[pt1][2] == path[pt0][1] + path[pt0][2] && rand(rng, Bool)
        try_mirror_path!(path, pt0, pt1)
        return
    end

    num_diag_before = 0
    num_diag_after = 0
    for i in 1:(pt1 - pt0 - 1)
        new_x = path[pt0][1] + path[pt1][1] - path[pt1 - i][1]
        new_y = path[pt0][2] + path[pt1][2] - path[pt1 - i][2]
        if new_x < 0 || new_y < 0
            return
        end
        if path[pt0 + i][1] == path[pt0 + i][2]
            num_diag_before += 1
        end
        if new_x == new_y
            num_diag_after += 1
        end
    end

    if num_diag_after > num_diag_before
        delta = num_diag_after - num_diag_before
        delta > 29 && return
        rand(rng, 0:(2^delta - 1)) == 0 || return
    end

    end_i = fld(pt1 - pt0, 2)
    for i in 1:end_i
        for dim in 1:2
            tmp = path[pt0][dim] + path[pt1][dim] - path[pt1 - i][dim]
            path[pt1 - i][dim] = path[pt0][dim] + path[pt1][dim] - path[pt0 + i][dim]
            path[pt0 + i][dim] = tmp
        end
    end
    return
end

function random_walk_in_cone(halflength::Integer, steps::Integer, rng::AbstractRNG)
    walk = random_walk_in_quarter_plane(halflength, rng)
    path = walk_to_path(walk)
    for _ in 1:Int(steps)
        markov_update_path!(path, rng)
    end
    return reflect_walk(path_to_walk(path))
end

function add_first_schnyder_tree!(walk, m::PlanarMap)
    cur = Int(m.root)
    for step in walk
        if step == RIGHT_STEP || step == DOWN_STEP
            cur = Int(get_next(m, insert_edge!(m, cur; color="orange"), 1))
        else
            cur = Int(get_next(m, cur, 1))
        end
    end
    return
end

function add_second_schnyder_tree!(walk, m::PlanarMap)
    cur = Int(get_next(m, m.root, 3))
    num_heads = [0]
    for step in walk
        if step == RIGHT_STEP || step == UP_STEP
            push!(num_heads, 0)
        else
            num_heads[end] += 1
        end
    end

    tails = Int[]
    green_corners = Int[]
    size_walk = length(walk)

    for i in 1:(size_walk + 1)
        cur = Int(get_next(m, cur, 1))
        if i == size_walk + 1 || walk[i] == RIGHT_STEP || walk[i] == DOWN_STEP
            heads = popfirst!(num_heads)
            while heads > 0
                insert_edge!(m, pop!(tails), cur; color="navy")
                heads -= 1
            end
            if i < size_walk
                push!(green_corners, Int(get_previous(m, cur, 1)))
            end
        end

        if i < size_walk && (walk[i + 1] == UP_STEP || walk[i + 1] == LEFT_STEP)
            push!(tails, cur)
        end
    end

    return green_corners
end

function build_schnyder_edge_data(m::PlanarMap)
    pairs = Dict{Tuple{Int32,Int32},String}()
    for eid in eachindex(m.edges)
        adj = Int(m.edges[eid].adj)
        eid > adj && continue
        u = m.edges[eid].vertex
        v = m.edges[adj].vertex
        u == v && continue
        a = min(u, v)
        b = max(u, v)
        if !haskey(pairs, (a, b))
            pairs[(a, b)] = string(m.edges[eid].color === nothing ? "black" : m.edges[eid].color)
        end
    end

    keys_sorted = sort!(collect(keys(pairs)))
    edge_u = Int32[]
    edge_v = Int32[]
    edge_color = String[]
    for (u, v) in keys_sorted
        push!(edge_u, u)
        push!(edge_v, v)
        push!(edge_color, pairs[(u, v)])
    end
    return edge_u, edge_v, edge_color
end

function generate_schnyder_map(; size_faces::Integer=100, seed::Integer=1)
    size = Int(size_faces)
    size >= 4 || throw(ArgumentError("Schnyder maps require at least 4 faces"))
    iseven(size) || throw(ArgumentError("Schnyder maps currently require an even number of faces"))

    rng = MersenneTwister(Int(seed))
    n = size ÷ 2 - 1
    steps = 40 * n
    walk = random_walk_in_cone(n, steps, rng)

    m = PlanarMap()
    make_polygon!(m, 3; color="outer")
    add_first_schnyder_tree!(walk, m)
    green_corners = add_second_schnyder_tree!(walk, m)
    for corner in green_corners
        insert_edge!(m, get_next(m, corner, 1), get_previous(m, corner, 1); color="green")
    end

    nverts = assign_vertex_ids!(m)
    faces = extract_faces(m)
    edge_u, edge_v, edge_color = build_schnyder_edge_data(m)
    return SchnyderMap(m, faces, edge_u, edge_v, edge_color, nverts)
end

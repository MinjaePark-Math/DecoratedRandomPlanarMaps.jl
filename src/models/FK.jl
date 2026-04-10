# FK-decorated map model (with spanning-tree p=0 as a special case).

const GREEN = Int8(0)
const RED = Int8(1)
const BLUE = Int8(2)
const PURPLE = Int8(3)
const ORANGE = Int8(4)

const FK_COLOR_NAME = Dict(
    Int(GREEN) => "green",
    Int(RED) => "red",
    Int(BLUE) => "blue",
    Int(PURPLE) => "purple",
    Int(ORANGE) => "orange",
)

struct FKMap <: AbstractRandomPlanarMap
    word::String
    resolved_word::String
    production_steps::Vector{Int32}
    order_steps::Vector{Int32}
    production_symbols::String
    order_symbols::String
    fulfilled_by::String
    resolved_order_symbols::String
    vertex_color::Vector{Int8}
    vertex_of_prod::Vector{Int32}
    predecessor_vertex::Vector{Int32}
    opposite_vertex_at_production::Vector{Int32}
    opposite_vertex_at_order::Vector{Int32}
    production_face::Vector{Int32}
    order_face::Vector{Int32}
    quadrilateral_faces::Matrix{Int32}
    triangulation_faces::Matrix{Int32}
    triangle_green_edges::Matrix{Int32}
    triangle_diag_edge::Vector{Int32}
    edge_u::Vector{Int32}
    edge_v::Vector{Int32}
    edge_color::Vector{Int8}
    edge_active::BitVector
    diagonal_before::Vector{Int32}
    diagonal_after::Vector{Int32}
    root_edge::Int32
end

const HCMap = FKMap
const FKDecoratedMap = FKMap
const SpanningTreeMap = FKMap

num_vertices(m::FKMap) = length(m.vertex_color)
num_faces(m::FKMap) = size(m.quadrilateral_faces, 1)
num_triangles(m::FKMap) = size(m.triangulation_faces, 1)
num_edges(m::FKMap; active_only::Bool=false) = active_only ? count(m.edge_active) : length(m.edge_u)

function active_edge_pairs(
    m::FKMap;
    colors=nothing,
    collapse::Bool=false,
    drop_loops::Bool=false,
)
    mask = copy(m.edge_active)
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

collapsed_green_edges(m::FKMap; drop_loops::Bool=true) = active_edge_pairs(m; colors=(GREEN,), collapse=true, drop_loops=drop_loops)
collapsed_active_edges(m::FKMap; drop_loops::Bool=true) = active_edge_pairs(m; collapse=true, drop_loops=drop_loops)
layout_edges(m::FKMap; drop_loops::Bool=true) = collapsed_active_edges(m; drop_loops=drop_loops)

function collapsed_active_edges_by_color(m::FKMap; drop_loops::Bool=true)
    out = Dict{Int,Matrix{Int32}}()
    for color in (GREEN, RED, BLUE, PURPLE, ORANGE)
        out[Int(color)] = active_edge_pairs(m; colors=(color,), collapse=true, drop_loops=drop_loops)
    end
    return out
end

function raw_active_edges_by_color(m::FKMap; drop_loops::Bool=true)
    out = Dict{Int,Matrix{Int32}}()
    for color in (GREEN, RED, BLUE, PURPLE, ORANGE)
        out[Int(color)] = active_edge_pairs(m; colors=(color,), collapse=false, drop_loops=drop_loops)
    end
    return out
end

function variant(m::FKMap)
    active_colors = Set(Int.(unique(m.edge_color[m.edge_active])))
    if (Int(PURPLE) in active_colors) || (Int(ORANGE) in active_colors) || occursin('F', m.order_symbols)
        return "fk"
    end
    return "spanning_tree"
end

function maximal_gasket_span(m::FKMap, gasket::AbstractString="h")
    gasket_norm = lowercase(strip(string(gasket)))
    gasket_norm in ("h", "c", "h_gasket", "c_gasket") || throw(ArgumentError("invalid gasket"))
    kind = gasket_norm[1]
    start_type = kind == 'h' ? 'c' : 'h'

    spans = Tuple{Int,Int}[]
    for i in eachindex(m.production_steps)
        if m.production_symbols[i] == start_type && m.order_symbols[i] == 'F'
            push!(spans, (Int(m.order_steps[i]) - Int(m.production_steps[i]), i - 1))
        end
    end
    isempty(spans) && throw(ArgumentError("no maximal gasket could be extracted"))
    _, prod_idx = findmax(spans)
    chosen = spans[prod_idx][2]
    return Int32(chosen), m.production_steps[chosen+1], m.order_steps[chosen+1]
end

function maximal_gasket_boundary(m::FKMap, gasket::AbstractString="h")
    gasket_norm = lowercase(strip(string(gasket)))
    gasket_norm in ("h", "c", "h_gasket", "c_gasket") || throw(ArgumentError("invalid gasket"))
    kind = gasket_norm[1]
    boundary_order = kind == 'h' ? 'H' : 'C'

    prod_idx, start_step, end_step = maximal_gasket_span(m, string(kind))
    root = Int(prod_idx)

    boundary_vertices = Int32[m.opposite_vertex_at_production[root+1]]
    order = sortperm(m.order_steps)
    for j in order
        order_step = Int(m.order_steps[j])
        if order_step <= Int(start_step) || order_step > Int(end_step)
            continue
        end
        if m.order_symbols[j] != boundary_order
            continue
        end
        if Int(m.production_steps[j]) >= Int(start_step)
            continue
        end
        push!(boundary_vertices, m.vertex_of_prod[j])
    end
    push!(boundary_vertices, m.opposite_vertex_at_order[root+1])

    unique_boundary = Int32[]
    seen = Set{Int32}()
    for v in boundary_vertices
        if v in seen
            continue
        end
        push!(seen, v)
        push!(unique_boundary, v)
    end
    return unique_boundary
end

struct QuadrantBridgeSampler
    n::Int
    cdf::Vector{Float64}
end

function QuadrantBridgeSampler(n::Integer)
    n_int = Int(n)
    n_int >= 1 || throw(ArgumentError("n must be >= 1"))

    weights = Vector{BigFloat}(undef, n_int + 1)
    weights[1] = big(1.0)
    for a in 0:(n_int-1)
        num = big(n_int - a) * big(n_int - a + 1)
        den = big(a + 1) * big(a + 2)
        weights[a+2] = weights[a+1] * (num / den)
    end
    total = sum(weights)
    cdf = Vector{Float64}(undef, n_int + 1)
    acc = big(0.0)
    for i in eachindex(weights)
        acc += weights[i] / total
        cdf[i] = Float64(acc)
    end
    cdf[end] = 1.0
    return QuadrantBridgeSampler(n_int, cdf)
end

function sample_uniform_bridge_word(s::QuadrantBridgeSampler, rng::AbstractRNG)
    n = s.n
    a = searchsortedfirst(s.cdf, rand(rng)) - 1

    x_steps = sample_uniform_dyck_path(a, rng)
    y_steps = sample_uniform_dyck_path(n - a, rng)

    marks = falses(2n)
    if 2a > 0
        for idx in randperm(rng, 2n)[1:(2a)]
            marks[idx] = true
        end
    end

    out = Vector{Char}(undef, 2n)
    i = 1
    j = 1
    for t in 1:(2n)
        if marks[t]
            st = x_steps[i]
            out[t] = st == 1 ? 'h' : 'H'
            i += 1
        else
            st = y_steps[j]
            out[t] = st == 1 ? 'c' : 'C'
            j += 1
        end
    end
    return String(out)
end

function sample_uniform_excursion_word(s::QuadrantBridgeSampler, rng::AbstractRNG)
    n = s.n
    while true
        word = sample_uniform_bridge_word(s, rng)
        x = 0
        y = 0
        ok = true
        for (t, ch) in enumerate(codeunits(word))
            c = Char(ch)
            if c == 'h'
                x += 1
            elseif c == 'H'
                x -= 1
            elseif c == 'c'
                y += 1
            else
                y -= 1
            end
            if t < 2n && x == 0 && y == 0
                ok = false
                break
            end
        end
        ok && return word
    end
end

function analyze_typed_word(word::AbstractString; require_excursion::Bool=true)
    chars = collect(replace(String(word), r"\s+" => ""))
    n_prod = count(c -> c == 'h' || c == 'c', chars)

    prev_stack = fill(-1, n_prod)
    next_stack = fill(-1, n_prod)
    prev_same = fill(-1, n_prod)
    prod_type = fill('\0', n_prod)

    top_any = -1
    top_by_type = Dict('h' => -1, 'c' => -1)
    created = 0
    positions = Int[]

    for (pos1, ch) in enumerate(chars)
        pos = pos1 - 1
        if ch == 'h' || ch == 'c'
            idx = created
            created += 1
            prod_type[idx+1] = ch

            prev_stack[idx+1] = top_any
            if top_any != -1
                next_stack[top_any+1] = idx
            end
            top_any = idx

            prev_same[idx+1] = top_by_type[ch]
            top_by_type[ch] = idx

        elseif ch == 'H' || ch == 'C'
            target_type = ch == 'H' ? 'h' : 'c'
            idx = top_by_type[target_type]
            idx == -1 && throw(ArgumentError("typed order has no earlier burger of the same type"))

            if idx == top_any
                push!(positions, pos)
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
        else
            throw(ArgumentError("invalid typed symbol $(repr(ch))"))
        end

        if require_excursion && pos < length(chars) - 1 && top_any == -1
            throw(ArgumentError("typed word returns to the empty stack before the final time"))
        end
    end

    if top_any != -1 || top_by_type['h'] != -1 || top_by_type['c'] != -1
        throw(ArgumentError("typed word does not reduce to the empty word"))
    end

    return length(positions), positions
end

function decorate_typed_word_with_fresh_orders(typed_word::AbstractString, eligible_positions, p::Real, rng::AbstractRNG)
    0.0 <= p <= 1.0 || throw(ArgumentError("p must lie in [0, 1]"))
    chars = collect(String(typed_word))
    if p == 0.0
        return String(chars)
    elseif p == 1.0
        for pos in eligible_positions
            chars[pos+1] = 'F'
        end
        return String(chars)
    end

    q = (2.0 * float(p)) / (1.0 + float(p))
    for pos in eligible_positions
        if rand(rng) < q
            chars[pos+1] = 'F'
        end
    end
    return String(chars)
end

function sample_exact_all_fresh(n::Integer, rng::AbstractRNG)
    n_int = Int(n)
    n_int >= 1 || throw(ArgumentError("n must be >= 1"))
    path = sample_uniform_primitive_dyck_path(n_int, rng)
    out = Vector{Char}(undef, 2n_int)
    for i in eachindex(path)
        if path[i] == 1
            out[i] = rand(rng) < 0.5 ? 'h' : 'c'
        else
            out[i] = 'F'
        end
    end
    return String(out)
end

function _weighted_choice_log(log_weights, rng::AbstractRNG)
    m = maximum(log_weights)
    weights = exp.(log_weights .- m)
    total = sum(weights)
    u = rand(rng) * total
    acc = 0.0
    for (i, w) in enumerate(weights)
        acc += w
        if u <= acc
            return i
        end
    end
    return length(weights)
end


function _normalize_fk_sampling_method(method)
    method_sym = method isa Symbol ? method : Symbol(replace(lowercase(strip(string(method))), '-' => '_'))
    method_sym in (:auto, :exact_rejection, :approx) || throw(
        ArgumentError("unknown FK sampling method $(repr(method)); expected one of :auto, :exact_rejection, :approx")
    )
    return method_sym
end

_smoothstep01(t::Real) = t <= 0 ? 0.0 : t >= 1 ? 1.0 : float(t)^2 * (3.0 - 2.0 * float(t))


function recommended_fk_approx_sampler_params(n::Integer, p::Real)
    n_int = Int(n)
    n_int >= 1 || throw(ArgumentError("n must be >= 1"))
    0.0 <= p <= 1.0 || throw(ArgumentError("p must lie in [0, 1]"))

    base_pool = clamp(round(Int, sqrt(float(n_int))), 64, 256)
    tilt = _smoothstep01(float(p) / 0.5)
    pool = clamp(round(Int, base_pool * (1.0 + tilt)), 64, 512)
    mh = max(16, 4 * pool)
    return pool, mh
end

function _resolve_fk_approx_sampler_params(n::Integer, p::Real, pool_size, mh_steps)
    rec_pool, _ = recommended_fk_approx_sampler_params(n, p)

    pool = pool_size === nothing ? rec_pool : Int(pool_size)
    pool >= 1 || throw(ArgumentError("pool_size must be >= 1 when provided"))

    mh = mh_steps === nothing ? max(16, 4 * pool) : Int(mh_steps)
    mh >= 0 || throw(ArgumentError("mh_steps must be >= 0 when provided"))

    return pool, mh
end

function estimate_fk_exact_rejection_acceptance(n::Integer, p::Real; seed=nothing, samples::Integer=64)
    n_int = Int(n)
    n_int >= 1 || throw(ArgumentError("n must be >= 1"))
    0.0 <= p <= 1.0 || throw(ArgumentError("p must lie in [0, 1]"))
    samples_int = max(1, Int(samples))
    rng = seed === nothing ? Random.default_rng() : MersenneTwister(Int(seed))

    if p == 0.0 || p == 1.0
        return 1.0
    end

    base = QuadrantBridgeSampler(n_int)
    log_r = log1p(float(p)) - log1p(-float(p))
    acc = 0.0
    for _ in 1:samples_int
        word = sample_uniform_excursion_word(base, rng)
        count, _ = analyze_typed_word(word; require_excursion=true)
        acc += exp((count - n_int) * log_r)
    end
    return acc / samples_int
end

function _sample_fk_word_exact_rejection(n::Integer, p::Real, rng::AbstractRNG; max_tries::Integer=1_000_000)
    n_int = Int(n)
    n_int >= 1 || throw(ArgumentError("n must be >= 1"))
    0.0 <= p <= 1.0 || throw(ArgumentError("p must lie in [0, 1]"))
    tries_int = max(1, Int(max_tries))

    if p == 1.0
        return sample_exact_all_fresh(n_int, rng)
    end

    base = QuadrantBridgeSampler(n_int)
    if p == 0.0
        return sample_uniform_excursion_word(base, rng)
    end

    log_r = log1p(float(p)) - log1p(-float(p))
    for _ in 1:tries_int
        typed_word = sample_uniform_excursion_word(base, rng)
        count, positions = analyze_typed_word(typed_word; require_excursion=true)
        if count == n_int || log(rand(rng)) <= (count - n_int) * log_r
            return decorate_typed_word_with_fresh_orders(typed_word, positions, p, rng)
        end
    end

    throw(ErrorException(
        "exact FK rejection sampler exceeded max_tries=$(tries_int). " *
        "Try a smaller size / smaller p, increase max_exact_tries, or use sampling_method=:approx."
    ))
end

function sample_fk_word_exact_rejection(n::Integer, p::Real; seed=nothing, max_tries::Integer=1_000_000)
    rng = seed === nothing ? Random.default_rng() : MersenneTwister(Int(seed))
    return _sample_fk_word_exact_rejection(n, p, rng; max_tries=max_tries)
end

function _sample_fk_word_approx(n::Integer, p::Real, rng::AbstractRNG; pool_size=nothing, mh_steps=nothing)
    n_int = Int(n)
    n_int >= 1 || throw(ArgumentError("n must be >= 1"))
    0.0 <= p <= 1.0 || throw(ArgumentError("p must lie in [0, 1]"))

    if p == 1.0
        return sample_exact_all_fresh(n_int, rng)
    end

    base = QuadrantBridgeSampler(n_int)
    if p == 0.0
        return sample_uniform_excursion_word(base, rng)
    end

    log_r = log1p(float(p)) - log1p(-float(p))
    pool_size_int, mh_steps_int = _resolve_fk_approx_sampler_params(n_int, p, pool_size, mh_steps)
    pool_words = String[]
    pool_counts = Int[]
    pool_positions = Vector{Vector{Int}}()

    for _ in 1:pool_size_int
        word = sample_uniform_excursion_word(base, rng)
        count, positions = analyze_typed_word(word; require_excursion=true)
        push!(pool_words, word)
        push!(pool_counts, count)
        push!(pool_positions, positions)
    end

    if pool_size_int == 1
        word = pool_words[1]
        top_count = pool_counts[1]
        eligible_positions = pool_positions[1]
    else
        idx = _weighted_choice_log([count * log_r for count in pool_counts], rng)
        word = pool_words[idx]
        top_count = pool_counts[idx]
        eligible_positions = pool_positions[idx]
    end

    for _ in 1:max(0, mh_steps_int)
        cand = sample_uniform_excursion_word(base, rng)
        cand_count, cand_positions = analyze_typed_word(cand; require_excursion=true)
        delta = cand_count - top_count
        if delta >= 0 || rand(rng) < exp(delta * log_r)
            word = cand
            top_count = cand_count
            eligible_positions = cand_positions
        end
    end

    return decorate_typed_word_with_fresh_orders(word, eligible_positions, p, rng)
end

function sample_fk_word(
    n::Integer,
    p::Real;
    seed=nothing,
    method=:auto,
    pool_size=nothing,
    mh_steps=nothing,
    exact_pilot_samples::Integer=64,
    exact_min_acceptance::Real=1.0e-4,
    max_exact_tries::Integer=1_000_000,
)
    n_int = Int(n)
    n_int >= 1 || throw(ArgumentError("n must be >= 1"))
    0.0 <= p <= 1.0 || throw(ArgumentError("p must lie in [0, 1]"))
    rng = seed === nothing ? Random.default_rng() : MersenneTwister(Int(seed))

    method_sym = _normalize_fk_sampling_method(method)
    if method_sym == :auto && 0.0 < p < 1.0
        est = estimate_fk_exact_rejection_acceptance(n_int, p; seed=seed, samples=exact_pilot_samples)
        method_sym = est >= float(exact_min_acceptance) ? :exact_rejection : :approx
    end

    if method_sym == :exact_rejection
        if method == :auto || lowercase(strip(string(method))) == "auto"
            try
                return _sample_fk_word_exact_rejection(n_int, p, rng; max_tries=max_exact_tries)
            catch err
                @warn string("exact FK rejection sampler failed in auto mode; falling back to approximate sampler: ", err)
                return _sample_fk_word_approx(n_int, p, rng; pool_size=pool_size, mh_steps=mh_steps)
            end
        end
        return _sample_fk_word_exact_rejection(n_int, p, rng; max_tries=max_exact_tries)
    end

    return _sample_fk_word_approx(n_int, p, rng; pool_size=pool_size, mh_steps=mh_steps)
end

const sample_hc_word = sample_fk_word
const sample_hc_word_exact_rejection = sample_fk_word_exact_rejection

function _resolve_and_pair_arrays(word::AbstractString; require_excursion::Bool=true)
    chars = collect(replace(String(word), r"\s+" => ""))
    all(c -> c in ['h', 'c', 'H', 'C', 'F'], chars) || throw(ArgumentError("word must use only h, c, H, C, F"))

    m = length(chars)
    n_prod = count(c -> c == 'h' || c == 'c', chars)
    m == 2n_prod || throw(ArgumentError("a valid finite map word must have as many productions as orders"))

    prod_steps = Vector{Int32}(undef, n_prod)
    order_steps = Vector{Int32}(undef, n_prod)
    prod_symbols = fill('\0', n_prod)
    order_symbols = fill('\0', n_prod)
    fulfilled_by = fill('\0', n_prod)
    resolved_order_symbols = fill('\0', n_prod)

    prev_stack = fill(-1, n_prod)
    next_stack = fill(-1, n_prod)
    prev_same = fill(-1, n_prod)
    top_any = -1
    top_by_type = Dict('h' => -1, 'c' => -1)
    created = 0
    resolved_chars = Char[]

    for (pos1, ch) in enumerate(chars)
        pos = pos1 - 1
        if ch == 'h' || ch == 'c'
            idx = created
            created += 1
            prod_steps[idx+1] = Int32(pos)
            prod_symbols[idx+1] = ch

            prev_stack[idx+1] = top_any
            if top_any != -1
                next_stack[top_any+1] = idx
            end
            top_any = idx

            prev_same[idx+1] = top_by_type[ch]
            top_by_type[ch] = idx
            push!(resolved_chars, ch)

        elseif ch == 'H' || ch == 'C'
            target_type = ch == 'H' ? 'h' : 'c'
            idx = top_by_type[target_type]
            idx == -1 && throw(ArgumentError("order $(repr(ch)) at step $pos cannot be fulfilled"))

            order_steps[idx+1] = Int32(pos)
            order_symbols[idx+1] = ch
            fulfilled_by[idx+1] = target_type
            resolved_order_symbols[idx+1] = ch

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
            push!(resolved_chars, ch)

        else
            idx = top_any
            idx == -1 && throw(ArgumentError("fresh order at step $pos cannot be fulfilled"))

            target_type = prod_symbols[idx+1]
            resolved_ch = target_type == 'h' ? 'H' : 'C'

            order_steps[idx+1] = Int32(pos)
            order_symbols[idx+1] = 'F'
            fulfilled_by[idx+1] = target_type
            resolved_order_symbols[idx+1] = resolved_ch

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
            push!(resolved_chars, resolved_ch)
        end

        if require_excursion && pos < m - 1 && top_any == -1
            throw(ArgumentError("inventory becomes empty before the final step; this word is not an excursion word"))
        end
    end

    top_any == -1 || throw(ArgumentError("the reduced word is not empty at the end"))

    return (
        String(chars),
        String(resolved_chars),
        prod_steps,
        order_steps,
        String(prod_symbols),
        String(order_symbols),
        String(fulfilled_by),
        String(resolved_order_symbols),
    )
end

function build_fk_map_from_word(word::AbstractString; require_excursion::Bool=true)
    word_s, resolved, production_steps, order_steps, production_symbols, order_symbols,
    fulfilled_by, resolved_order_symbols = _resolve_and_pair_arrays(word; require_excursion=require_excursion)

    m = length(word_s)
    n_prod = length(production_steps)

    tri_faces = fill(Int32(-1), m, 3)
    tri_green_edges = fill(Int32(-1), m, 2)
    tri_diag_edge = fill(Int32(-1), m)

    edge_u = Int32[]
    edge_v = Int32[]
    edge_color = Int8[]
    edge_active = BitVector()

    function new_edge(u::Integer, v::Integer, color::Integer, active::Bool=true)
        push!(edge_u, Int32(u))
        push!(edge_v, Int32(v))
        push!(edge_color, Int8(color))
        push!(edge_active, active)
        return Int32(length(edge_u) - 1)
    end

    root_edge = new_edge(0, 1, Int(GREEN), true)

    red_vertices = Int32[0]
    blue_vertices = Int32[1]
    red_edges = Int32[]
    blue_edges = Int32[]
    red_prod_indices = Int32[]
    blue_prod_indices = Int32[]
    special_edge = root_edge

    vertex_color = Int8[0, 1]
    next_vid = 2

    vertex_of_prod = fill(Int32(-1), n_prod)
    predecessor_vertex = fill(Int32(-1), n_prod)
    opposite_vertex_at_production = fill(Int32(-1), n_prod)
    opposite_vertex_at_order = fill(Int32(-1), n_prod)
    production_face = fill(Int32(-1), n_prod)
    order_face = fill(Int32(-1), n_prod)
    diagonal_before = fill(Int32(-1), n_prod)
    diagonal_after = fill(Int32(-1), n_prod)

    old_special_at_production = fill(Int32(-1), n_prod)
    new_special_at_production = fill(Int32(-1), n_prod)
    old_special_at_order = fill(Int32(-1), n_prod)
    new_special_at_order = fill(Int32(-1), n_prod)

    next_prod_idx = 0

    resolved_chars = collect(resolved)
    for (step1, ch) in enumerate(resolved_chars)
        step = step1 - 1
        if ch == 'h'
            prod_idx = next_prod_idx
            next_prod_idx += 1
            pred = red_vertices[end]
            opp = blue_vertices[end]
            v = next_vid
            next_vid += 1
            push!(vertex_color, Int8(0))

            diag = new_edge(pred, v, Int(RED), true)
            new_special = new_edge(v, opp, Int(GREEN), true)

            tri_faces[step1, :] = Int32[pred, opp, v]
            tri_green_edges[step1, :] = Int32[special_edge, new_special]
            tri_diag_edge[step1] = diag

            push!(red_vertices, Int32(v))
            push!(red_edges, diag)
            push!(red_prod_indices, Int32(prod_idx))

            vertex_of_prod[prod_idx+1] = Int32(v)
            predecessor_vertex[prod_idx+1] = Int32(pred)
            opposite_vertex_at_production[prod_idx+1] = Int32(opp)
            production_face[prod_idx+1] = Int32(step)
            old_special_at_production[prod_idx+1] = special_edge
            new_special_at_production[prod_idx+1] = new_special
            diagonal_before[prod_idx+1] = diag
            special_edge = new_special

        elseif ch == 'c'
            prod_idx = next_prod_idx
            next_prod_idx += 1
            pred = blue_vertices[end]
            opp = red_vertices[end]
            v = next_vid
            next_vid += 1
            push!(vertex_color, Int8(1))

            diag = new_edge(pred, v, Int(BLUE), true)
            new_special = new_edge(opp, v, Int(GREEN), true)

            tri_faces[step1, :] = Int32[opp, pred, v]
            tri_green_edges[step1, :] = Int32[special_edge, new_special]
            tri_diag_edge[step1] = diag

            push!(blue_vertices, Int32(v))
            push!(blue_edges, diag)
            push!(blue_prod_indices, Int32(prod_idx))

            vertex_of_prod[prod_idx+1] = Int32(v)
            predecessor_vertex[prod_idx+1] = Int32(pred)
            opposite_vertex_at_production[prod_idx+1] = Int32(opp)
            production_face[prod_idx+1] = Int32(step)
            old_special_at_production[prod_idx+1] = special_edge
            new_special_at_production[prod_idx+1] = new_special
            diagonal_before[prod_idx+1] = diag
            special_edge = new_special

        elseif ch == 'H'
            prod_idx = Int(pop!(red_prod_indices))
            _v = pop!(red_vertices)
            pred = red_vertices[end]
            opp = blue_vertices[end]
            _diag = pop!(red_edges)

            new_special = (length(red_vertices) == 1 && length(blue_vertices) == 1) ? root_edge : new_edge(opp, pred, Int(GREEN), true)

            tri_faces[step1, :] = Int32[opp, pred, _v]
            tri_green_edges[step1, :] = Int32[special_edge, new_special]
            tri_diag_edge[step1] = _diag

            opposite_vertex_at_order[prod_idx+1] = Int32(opp)
            order_face[prod_idx+1] = Int32(step)
            old_special_at_order[prod_idx+1] = special_edge
            new_special_at_order[prod_idx+1] = new_special
            special_edge = new_special

        elseif ch == 'C'
            prod_idx = Int(pop!(blue_prod_indices))
            _v = pop!(blue_vertices)
            pred = blue_vertices[end]
            opp = red_vertices[end]
            _diag = pop!(blue_edges)

            new_special = (length(red_vertices) == 1 && length(blue_vertices) == 1) ? root_edge : new_edge(opp, pred, Int(GREEN), true)

            tri_faces[step1, :] = Int32[opp, pred, _v]
            tri_green_edges[step1, :] = Int32[special_edge, new_special]
            tri_diag_edge[step1] = _diag

            opposite_vertex_at_order[prod_idx+1] = Int32(opp)
            order_face[prod_idx+1] = Int32(step)
            old_special_at_order[prod_idx+1] = special_edge
            new_special_at_order[prod_idx+1] = new_special
            special_edge = new_special
        end
    end

    quad_faces = Matrix{Int32}(undef, n_prod, 4)
    for i in 1:n_prod
        quad_faces[i, :] = Int32[
            opposite_vertex_at_production[i],
            vertex_of_prod[i],
            opposite_vertex_at_order[i],
            predecessor_vertex[i],
        ]

        if order_symbols[i] == 'F'
            new_color = fulfilled_by[i] == 'h' ? PURPLE : ORANGE
            new_diag = new_edge(opposite_vertex_at_production[i], opposite_vertex_at_order[i], Int(new_color), true)
            diagonal_after[i] = new_diag
            edge_active[Int(diagonal_before[i])+1] = false

            tri_faces[Int(production_face[i])+1, :] = Int32[
                vertex_of_prod[i],
                opposite_vertex_at_production[i],
                opposite_vertex_at_order[i],
            ]
            tri_faces[Int(order_face[i])+1, :] = Int32[
                predecessor_vertex[i],
                opposite_vertex_at_production[i],
                opposite_vertex_at_order[i],
            ]
            tri_green_edges[Int(production_face[i])+1, :] = Int32[new_special_at_production[i], old_special_at_order[i]]
            tri_green_edges[Int(order_face[i])+1, :] = Int32[old_special_at_production[i], new_special_at_order[i]]
            tri_diag_edge[Int(production_face[i])+1] = new_diag
            tri_diag_edge[Int(order_face[i])+1] = new_diag
        else
            diagonal_after[i] = diagonal_before[i]
        end
    end

    return FKMap(
        word_s,
        resolved,
        production_steps,
        order_steps,
        production_symbols,
        order_symbols,
        fulfilled_by,
        resolved_order_symbols,
        vertex_color,
        vertex_of_prod,
        predecessor_vertex,
        opposite_vertex_at_production,
        opposite_vertex_at_order,
        production_face,
        order_face,
        quad_faces,
        tri_faces,
        tri_green_edges,
        tri_diag_edge,
        edge_u,
        edge_v,
        edge_color,
        edge_active,
        diagonal_before,
        diagonal_after,
        root_edge,
    )
end

const build_hc_map_from_word = build_fk_map_from_word

function _resolve_p_q(; p=nothing, q=nothing)
    if q !== nothing
        qf = float(q)
        if isinf(qf)
            pval = 1.0
        else
            qf >= 0.0 || throw(ArgumentError("q must be nonnegative"))
            s = sqrt(qf)
            pval = s / (2.0 + s)
        end
        if p !== nothing
            abs(float(p) - pval) <= 1e-12 || throw(ArgumentError("p and q specify inconsistent FK parameters"))
        end
        return pval
    end
    return p === nothing ? 0.25 : float(p)
end

function build_fk_map(
    ;
    faces::Integer,
    p=nothing,
    q=nothing,
    seed::Integer=1,
    sampling_method=:auto,
    pool_size=nothing,
    mh_steps=nothing,
    exact_pilot_samples::Integer=64,
    exact_min_acceptance::Real=1.0e-4,
    max_exact_tries::Integer=1_000_000,
)
    pval = _resolve_p_q(; p=p, q=q)
    word = sample_fk_word(
        Int(faces),
        pval;
        seed=seed,
        method=sampling_method,
        pool_size=pool_size,
        mh_steps=mh_steps,
        exact_pilot_samples=exact_pilot_samples,
        exact_min_acceptance=exact_min_acceptance,
        max_exact_tries=max_exact_tries,
    )
    return build_fk_map_from_word(word; require_excursion=true)
end

const build_hc_map = build_fk_map

build_spanning_tree_map(; faces::Integer, seed::Integer=1) = build_fk_map(; faces=faces, p=0.0, seed=seed)

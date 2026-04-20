# Mated-CRT maps sampled from correlated Brownian excursions / bridges.

struct MatedCRTMap <: AbstractRandomPlanarMap
    topology::Symbol
    nverts::Int
    gamma::Float64
    gamma_prime::Float64
    kappa::Float64
    kappa_prime::Float64
    brownian_correlation::Float64
    endpoint::NTuple{2,Float64}
    refinement::Int
    sampler::String
    time_grid::Vector{Float64}
    L::Vector{Float64}
    R::Vector{Float64}
    faces::Vector{Vector{Int32}}
    face_edge_ids::Vector{Vector{Int32}}
    boundary_vertices::Vector{Int32}
    outer_face_index::Int32
    edge_u::Vector{Int32}
    edge_v::Vector{Int32}
    edge_kind::Vector{Symbol}
    sphere_layout_map::Union{Nothing,FKMap}
end

const MatedCRTPlanarMap = MatedCRTMap

num_vertices(m::MatedCRTMap) = m.nverts
num_faces(m::MatedCRTMap) = length(m.faces)

function layout_edges(m::MatedCRTMap; drop_loops::Bool=true)
    pairs = hcat(m.edge_u, m.edge_v)
    return first(collapse_undirected_edges(pairs; drop_loops=drop_loops))
end

function _mated_crt_topology(x)
    raw = lowercase(strip(string(x)))
    raw = replace(raw, '-' => '_')
    if raw in ("disk", "disc")
        return :disk
    elseif raw in ("sphere", "spherical")
        return :sphere
    end
    throw(ArgumentError("mated_crt topology must be `disk` or `sphere`"))
end

function _mated_crt_gamma_from_correlation(correlation::Real)
    rho = float(correlation)
    -1.0 <= rho <= 1.0 || throw(ArgumentError("correlation must lie in [-1, 1]"))
    return sqrt((4.0 / π) * acos(-rho))
end

function _resolve_mated_crt_parameters(;
    gamma=nothing,
    gamma_prime=nothing,
    kappa=nothing,
    kappa_prime=nothing,
    correlation=nothing,
)
    gamma_candidates = Float64[]

    if gamma !== nothing
        γ = float(gamma)
        0.0 < γ < 2.0 || throw(ArgumentError("gamma must lie in (0, 2)"))
        push!(gamma_candidates, γ)
    end

    if gamma_prime !== nothing
        γp = float(gamma_prime)
        γp > 2.0 || throw(ArgumentError("gamma_prime must be > 2"))
        push!(gamma_candidates, 4.0 / γp)
    end

    if kappa !== nothing
        κ = float(kappa)
        0.0 < κ < 4.0 || throw(ArgumentError("kappa must lie in (0, 4); use kappa_prime for the space-filling parameter"))
        push!(gamma_candidates, sqrt(κ))
    end

    if kappa_prime !== nothing
        κp = float(kappa_prime)
        κp > 4.0 || throw(ArgumentError("kappa_prime must be > 4"))
        push!(gamma_candidates, 4.0 / sqrt(κp))
    end

    if correlation !== nothing
        push!(gamma_candidates, _mated_crt_gamma_from_correlation(correlation))
    end

    γ = isempty(gamma_candidates) ? sqrt(2.0) : gamma_candidates[1]
    for candidate in gamma_candidates
        abs(candidate - γ) <= 1.0e-8 || throw(ArgumentError("mated_crt parameters specify inconsistent values for gamma / kappa / kappa_prime / gamma_prime / correlation"))
    end

    0.0 < γ < 2.0 || throw(ArgumentError("resolved gamma must lie in (0, 2)"))
    κ = γ^2
    κp = 16.0 / κ
    γp = 4.0 / γ
    rho = -cos(π * κ / 4.0)

    return Dict{String,Float64}(
        "gamma" => γ,
        "gamma_prime" => γp,
        "kappa" => κ,
        "kappa_prime" => κp,
        "correlation" => rho,
    )
end

function _mated_crt_initial_bridge(num_steps::Integer, endpoint::NTuple{2,Float64})
    m = Int(num_steps)
    times = collect(range(0.0, 1.0; length=m + 1))
    path = zeros(Float64, m + 1, 2)

    amp1 = max(0.45, 0.6 + 0.35 * endpoint[1])
    amp2 = max(0.45, 0.65 + 0.35 * endpoint[2])

    for i in eachindex(times)
        t = times[i]
        bump = sin(π * t)
        path[i, 1] = endpoint[1] * t + amp1 * bump
        path[i, 2] = endpoint[2] * t + amp2 * bump
    end

    path[1, :] .= 0.0
    path[end, 1] = endpoint[1]
    path[end, 2] = endpoint[2]
    return path
end

function _sample_truncated_bivariate_normal(
    mean::AbstractVector{<:Real},
    chol::AbstractMatrix{<:Real},
    rng::AbstractRNG;
    lower::Real=1.0e-9,
    max_tries::Integer=4_096,
)
    μ1 = float(mean[1])
    μ2 = float(mean[2])
    l11 = float(chol[1, 1])
    l21 = float(chol[2, 1])
    l22 = float(chol[2, 2])

    for _ in 1:Int(max_tries)
        z1 = randn(rng)
        z2 = randn(rng)
        x1 = μ1 + l11 * z1
        x2 = μ2 + l21 * z1 + l22 * z2
        if x1 > lower && x2 > lower
            return x1, x2
        end
    end

    # A conservative fallback keeps the Gibbs chain moving if rejection is unlucky.
    return max(μ1, lower), max(μ2, lower)
end

function _sample_positive_bridge_path(
    endpoint::NTuple{2,Float64},
    correlation::Real,
    num_steps::Integer,
    rng::AbstractRNG;
    burnin_sweeps::Integer=64,
    gibbs_sweeps::Integer=128,
    lower::Real=1.0e-9,
)
    m = Int(num_steps)
    m >= 2 || throw(ArgumentError("num_steps must be >= 2"))
    rho = float(correlation)
    abs(rho) < 1.0 || throw(ArgumentError("correlation must lie in (-1, 1)"))

    path = _mated_crt_initial_bridge(m, endpoint)
    cov = Symmetric((0.5 / m) * Float64[1.0 rho; rho 1.0])
    chol = cholesky(cov).L
    total_sweeps = max(0, Int(burnin_sweeps)) + max(0, Int(gibbs_sweeps))

    for sweep in 1:total_sweeps
        if isodd(sweep)
            for idx in 2:m
                mean = 0.5 .* (path[idx - 1, :] .+ path[idx + 1, :])
                x1, x2 = _sample_truncated_bivariate_normal(mean, chol, rng; lower=lower)
                path[idx, 1] = x1
                path[idx, 2] = x2
            end
        else
            for idx in m:-1:2
                mean = 0.5 .* (path[idx - 1, :] .+ path[idx + 1, :])
                x1, x2 = _sample_truncated_bivariate_normal(mean, chol, rng; lower=lower)
                path[idx, 1] = x1
                path[idx, 2] = x2
            end
        end
    end

    path[1, :] .= 0.0
    path[end, 1] = endpoint[1]
    path[end, 2] = endpoint[2]
    return path
end

function _mated_crt_sphere_sampler(x)
    raw = lowercase(strip(string(x)))
    raw = replace(raw, '-' => '_')
    if raw in ("exact", "cone_walk_exact", "exact_cone_walk")
        return :exact
    elseif raw in ("approx", "approximate", "hybrid", "approx_hybrid", "hybrid_cone_walk")
        return :approx
    end
    throw(ArgumentError("mated_crt sphere sampler must be `exact` or `approx`"))
end

function _mated_crt_sample_weighted_index(weights::AbstractVector{<:Real}, rng::AbstractRNG)
    total = 0.0
    for weight in weights
        total += max(float(weight), 0.0)
    end
    total > 0.0 || throw(ArgumentError("weighted sampling requires a positive total weight"))

    threshold = rand(rng) * total
    partial = 0.0
    for (idx, weight) in enumerate(weights)
        partial += max(float(weight), 0.0)
        if threshold <= partial
            return idx
        end
    end
    return length(weights)
end

function _mated_crt_logsumexp(values::AbstractVector{<:Real})
    isempty(values) && return -Inf
    m = maximum(float.(values))
    isfinite(m) || return m
    total = 0.0
    for value in values
        total += exp(float(value) - m)
    end
    return m + log(total)
end

function _mated_crt_lazy_correlated_steps(correlation::Real)
    rho = float(correlation)
    abs(rho) < 1.0 || throw(ArgumentError("exact cone-walk sphere sampler requires correlation in (-1, 1)"))

    q = abs(rho) / (2.0 - abs(rho))
    cardinal_weight = (1.0 - q) / 4.0
    diagonal_weight = q / 2.0

    steps = NamedTuple{(:dx, :dy, :prod, :cost, :weight, :logw),Tuple{Int,Int,Int,Int,Float64,Float64}}[]

    if cardinal_weight > 0.0
        for (dx, dy) in ((1, 0), (-1, 0), (0, 1), (0, -1))
            prod = max(dx, 0) + max(dy, 0)
            push!(steps, (
                dx=dx,
                dy=dy,
                prod=prod,
                cost=abs(dx) + abs(dy),
                weight=cardinal_weight,
                logw=log(cardinal_weight),
            ))
        end
    end

    if diagonal_weight > 0.0
        diagonals = rho > 0.0 ? ((1, 1), (-1, -1)) : ((1, -1), (-1, 1))
        for (dx, dy) in diagonals
            prod = max(dx, 0) + max(dy, 0)
            push!(steps, (
                dx=dx,
                dy=dy,
                prod=prod,
                cost=abs(dx) + abs(dy),
                weight=diagonal_weight,
                logw=log(diagonal_weight),
            ))
        end
    end

    isempty(steps) && throw(ArgumentError("exact cone-walk sphere sampler needs at least one admissible step"))
    return steps
end

function _mated_crt_exact_cone_logweight!(
    cache::Dict{NTuple{3,Int},Float64},
    steps,
    remaining_productions::Integer,
    x::Integer,
    y::Integer,
)
    p = Int(remaining_productions)
    xx = Int(x)
    yy = Int(y)
    xx >= 0 || return -Inf
    yy >= 0 || return -Inf
    p >= 0 || return -Inf

    key = (p, xx, yy)
    if haskey(cache, key)
        return cache[key]
    end

    if p == 0 && xx == 0 && yy == 0
        cache[key] = 0.0
        return 0.0
    end

    log_terms = Float64[]
    for step in steps
        p_next = p - step.prod
        p_next < 0 && continue

        x_next = xx + step.dx
        y_next = yy + step.dy
        if x_next < 0 || y_next < 0
            continue
        end
        if x_next == 0 && y_next == 0 && p_next != 0
            continue
        end

        child = _mated_crt_exact_cone_logweight!(cache, steps, p_next, x_next, y_next)
        isfinite(child) || continue
        push!(log_terms, step.logw + child)
    end

    value = isempty(log_terms) ? -Inf : _mated_crt_logsumexp(log_terms)
    cache[key] = value
    return value
end

function _mated_crt_sample_exact_cone_step!(
    cache::Dict{NTuple{3,Int},Float64},
    steps,
    remaining_productions::Integer,
    x::Integer,
    y::Integer,
    rng::AbstractRNG,
)
    p = Int(remaining_productions)
    xx = Int(x)
    yy = Int(y)

    step_ids = Int[]
    log_weights = Float64[]
    next_states = NTuple{3,Int}[]

    for (idx, step) in enumerate(steps)
        p_next = p - step.prod
        p_next < 0 && continue

        x_next = xx + step.dx
        y_next = yy + step.dy
        if x_next < 0 || y_next < 0
            continue
        end
        if x_next == 0 && y_next == 0 && p_next != 0
            continue
        end

        child = _mated_crt_exact_cone_logweight!(cache, steps, p_next, x_next, y_next)
        isfinite(child) || continue

        push!(step_ids, idx)
        push!(log_weights, step.logw + child)
        push!(next_states, (p_next, x_next, y_next))
    end

    isempty(step_ids) && throw(ArgumentError("failed to continue the exact cone-walk sphere excursion"))

    max_log = maximum(log_weights)
    weights = Float64[exp(logw - max_log) for logw in log_weights]
    chosen = _mated_crt_sample_weighted_index(weights, rng)
    return steps[step_ids[chosen]], next_states[chosen]
end

function _mated_crt_append_exact_suffix!(
    word::Vector{UInt8},
    cache::Dict{NTuple{3,Int},Float64},
    steps,
    remaining_productions::Integer,
    x::Integer,
    y::Integer,
    rng::AbstractRNG,
)
    p = Int(remaining_productions)
    xx = Int(x)
    yy = Int(y)
    while !(p == 0 && xx == 0 && yy == 0)
        step, (p_next, x_next, y_next) = _mated_crt_sample_exact_cone_step!(cache, steps, p, xx, yy, rng)
        _mated_crt_append_exact_step_block!(word, xx, yy, step.dx, step.dy, rng)
        p = p_next
        xx = x_next
        yy = y_next
    end
    return word
end

function _mated_crt_append_exact_step_block!(
    word::Vector{UInt8},
    x::Integer,
    y::Integer,
    dx::Integer,
    dy::Integer,
    rng::AbstractRNG,
)
    xx = Int(x)
    yy = Int(y)

    if dx == 1 && dy == 0
        push!(word, UInt8('h'))
        return
    elseif dx == -1 && dy == 0
        push!(word, UInt8('H'))
        return
    elseif dx == 0 && dy == 1
        push!(word, UInt8('c'))
        return
    elseif dx == 0 && dy == -1
        push!(word, UInt8('C'))
        return
    elseif dx == 1 && dy == 1
        if rand(rng, Bool)
            push!(word, UInt8('h'))
            push!(word, UInt8('c'))
        else
            push!(word, UInt8('c'))
            push!(word, UInt8('h'))
        end
        return
    elseif dx == -1 && dy == -1
        if rand(rng, Bool)
            push!(word, UInt8('H'))
            push!(word, UInt8('C'))
        else
            push!(word, UInt8('C'))
            push!(word, UInt8('H'))
        end
        return
    elseif dx == 1 && dy == -1
        choose_reverse = !(xx == 0 && yy == 1) && rand(rng, Bool)
        if choose_reverse
            push!(word, UInt8('C'))
            push!(word, UInt8('h'))
        else
            push!(word, UInt8('h'))
            push!(word, UInt8('C'))
        end
        return
    elseif dx == -1 && dy == 1
        choose_reverse = !(xx == 1 && yy == 0) && rand(rng, Bool)
        if choose_reverse
            push!(word, UInt8('H'))
            push!(word, UInt8('c'))
        else
            push!(word, UInt8('c'))
            push!(word, UInt8('H'))
        end
        return
    end

    throw(ArgumentError("unsupported exact cone-walk step ($dx, $dy)"))
end

function _mated_crt_exact_word_coordinates(word::AbstractString)
    chars = codeunits(word)
    L = Vector{Float64}(undef, length(chars) + 1)
    R = Vector{Float64}(undef, length(chars) + 1)
    x = 0
    y = 0
    L[1] = 0.0
    R[1] = 0.0

    for (idx, raw) in enumerate(chars)
        ch = Char(raw)
        if ch == 'h'
            x += 1
        elseif ch == 'H'
            x -= 1
        elseif ch == 'c'
            y += 1
        elseif ch == 'C'
            y -= 1
        else
            throw(ArgumentError("invalid exact cone-walk letter $(repr(ch))"))
        end
        x >= 0 || throw(ArgumentError("exact cone-walk word left the quadrant in the first coordinate"))
        y >= 0 || throw(ArgumentError("exact cone-walk word left the quadrant in the second coordinate"))
        L[idx + 1] = float(x)
        R[idx + 1] = float(y)
    end

    x == 0 || throw(ArgumentError("exact cone-walk word does not return to zero in the first coordinate"))
    y == 0 || throw(ArgumentError("exact cone-walk word does not return to zero in the second coordinate"))
    return L, R
end

function _mated_crt_sample_exact_sphere_word(
    target_productions::Integer,
    correlation::Real,
    rng::AbstractRNG,
)
    target = Int(target_productions)
    target >= 1 || throw(ArgumentError("exact cone-walk sphere sampler needs at least one production"))

    rho = float(correlation)
    if abs(rho) <= 1.0e-12
        return sample_uniform_excursion_word(QuadrantBridgeSampler(target), rng)
    end

    steps = _mated_crt_lazy_correlated_steps(rho)
    cache = Dict{NTuple{3,Int},Float64}()
    root = _mated_crt_exact_cone_logweight!(cache, steps, target, 0, 0)
    isfinite(root) || throw(ArgumentError("failed to sample an exact correlated cone excursion"))

    word = UInt8[]
    sizehint!(word, 2 * target)
    _mated_crt_append_exact_suffix!(word, cache, steps, target, 0, 0, rng)
    length(word) == 2 * target || throw(ArgumentError("exact cone-walk sphere word has the wrong length"))
    return String(Char.(word))
end

_mated_crt_remaining_letters(p::Integer, x::Integer, y::Integer) = 2 * Int(p) + Int(x) + Int(y)

function _mated_crt_cone_model(steps; harmonic_shift::NTuple{2,Float64}=(1.0, 1.0))
    cov = zeros(Float64, 2, 2)
    mean_cost = 0.0
    total_weight = 0.0
    for step in steps
        total_weight += step.weight
        mean_cost += step.weight * step.cost
        cov[1, 1] += step.weight * step.dx * step.dx
        cov[1, 2] += step.weight * step.dx * step.dy
        cov[2, 1] += step.weight * step.dy * step.dx
        cov[2, 2] += step.weight * step.dy * step.dy
    end
    total_weight > 0.0 || throw(ArgumentError("approximate cone model requires positive step weight"))
    mean_cost > 0.0 || throw(ArgumentError("approximate cone model requires positive mean step cost"))
    cov ./= total_weight
    mean_cost /= total_weight

    eig = eigen(Symmetric(cov))
    minimum(eig.values) > 0.0 || throw(ArgumentError("approximate cone model requires positive-definite covariance"))
    inv_sqrt_cov = eig.vectors * Diagonal(1.0 ./ sqrt.(eig.values)) * eig.vectors'
    inv_cov = eig.vectors * Diagonal(1.0 ./ eig.values) * eig.vectors'

    ray_x = inv_sqrt_cov * [1.0, 0.0]
    ray_y = inv_sqrt_cov * [0.0, 1.0]
    cross_xy = ray_x[1] * ray_y[2] - ray_x[2] * ray_y[1]
    dot_xy = dot(ray_x, ray_y)
    α = atan(abs(cross_xy), dot_xy)
    0.0 < α < π || throw(ArgumentError("failed to construct a valid cone angle for the approximate sphere sampler"))
    ν = π / α

    return (
        cov=cov,
        inv_cov=inv_cov,
        inv_sqrt_cov=inv_sqrt_cov,
        ray_x=ray_x,
        α=α,
        ν=ν,
        mean_cost=mean_cost,
        harmonic_shift=harmonic_shift,
    )
end

function _mated_crt_approx_logharmonic(model, x::Integer, y::Integer)
    shift = model.harmonic_shift
    shifted = model.inv_sqrt_cov * [Int(x) + shift[1], Int(y) + shift[2]]
    r = hypot(shifted[1], shifted[2])
    r > 0.0 || return -Inf

    ray = model.ray_x
    cross_ru = ray[1] * shifted[2] - ray[2] * shifted[1]
    dot_ru = dot(ray, shifted)
    θ = atan(max(cross_ru, 0.0), dot_ru)
    θ = clamp(θ, 1.0e-12, model.α - 1.0e-12)
    sin_term = sin(model.ν * θ)
    sin_term > 0.0 || return -Inf
    return model.ν * log(r) + log(sin_term)
end

function _mated_crt_approx_cone_logscore(
    exact_cache::Dict{NTuple{3,Int},Float64},
    steps,
    model,
    remaining_productions::Integer,
    x::Integer,
    y::Integer,
    exact_tail_cutoff::Integer,
)
    p = Int(remaining_productions)
    xx = Int(x)
    yy = Int(y)
    if p < 0 || xx < 0 || yy < 0
        return -Inf
    elseif p == 0 && xx == 0 && yy == 0
        return 0.0
    elseif xx == 0 && yy == 0
        return -Inf
    end

    if p <= Int(exact_tail_cutoff)
        return _mated_crt_exact_cone_logweight!(exact_cache, steps, p, xx, yy)
    end

    remaining_letters = _mated_crt_remaining_letters(p, xx, yy)
    remaining_letters > 0 || return -Inf
    effective_steps = remaining_letters / model.mean_cost
    effective_steps > 0.0 || return -Inf

    logharm = _mated_crt_approx_logharmonic(model, xx, yy)
    isfinite(logharm) || return -Inf

    z = [float(xx), float(yy)]
    quad = dot(z, model.inv_cov * z)
    return logharm - (model.ν + 1.0) * log(effective_steps) - quad / (2.0 * effective_steps)
end

function _mated_crt_sample_approx_cone_step!(
    exact_cache::Dict{NTuple{3,Int},Float64},
    steps,
    model,
    remaining_productions::Integer,
    x::Integer,
    y::Integer,
    exact_tail_cutoff::Integer,
    rng::AbstractRNG,
)
    p = Int(remaining_productions)
    xx = Int(x)
    yy = Int(y)

    step_ids = Int[]
    log_weights = Float64[]
    next_states = NTuple{3,Int}[]

    for (idx, step) in enumerate(steps)
        p_next = p - step.prod
        p_next < 0 && continue

        x_next = xx + step.dx
        y_next = yy + step.dy
        if x_next < 0 || y_next < 0
            continue
        end
        if x_next == 0 && y_next == 0 && p_next != 0
            continue
        end

        child = _mated_crt_approx_cone_logscore(
            exact_cache,
            steps,
            model,
            p_next,
            x_next,
            y_next,
            exact_tail_cutoff,
        )
        isfinite(child) || continue

        push!(step_ids, idx)
        push!(log_weights, step.logw + child)
        push!(next_states, (p_next, x_next, y_next))
    end

    isempty(step_ids) && return nothing

    max_log = maximum(log_weights)
    weights = Float64[exp(logw - max_log) for logw in log_weights]
    chosen = _mated_crt_sample_weighted_index(weights, rng)
    return steps[step_ids[chosen]], next_states[chosen]
end

function _mated_crt_sample_approx_sphere_word(
    target_productions::Integer,
    correlation::Real,
    rng::AbstractRNG;
    exact_tail_cutoff::Integer=64,
)
    target = Int(target_productions)
    target >= 1 || throw(ArgumentError("approximate sphere sampler needs at least one production"))
    tail_cutoff = max(0, Int(exact_tail_cutoff))

    rho = float(correlation)
    if abs(rho) <= 1.0e-12
        return sample_uniform_excursion_word(QuadrantBridgeSampler(target), rng)
    end

    steps = _mated_crt_lazy_correlated_steps(rho)
    model = _mated_crt_cone_model(steps)
    exact_cache = Dict{NTuple{3,Int},Float64}()

    word = UInt8[]
    sizehint!(word, 2 * target)
    p = target
    x = 0
    y = 0

    while !(p == 0 && x == 0 && y == 0)
        if p <= tail_cutoff
            _mated_crt_append_exact_suffix!(word, exact_cache, steps, p, x, y, rng)
            p = 0
            x = 0
            y = 0
            break
        end

        sampled = _mated_crt_sample_approx_cone_step!(
            exact_cache,
            steps,
            model,
            p,
            x,
            y,
            tail_cutoff,
            rng,
        )

        if sampled === nothing
            _mated_crt_append_exact_suffix!(word, exact_cache, steps, p, x, y, rng)
            p = 0
            x = 0
            y = 0
            break
        end

        step, (p_next, x_next, y_next) = sampled
        _mated_crt_append_exact_step_block!(word, x, y, step.dx, step.dy, rng)
        p = p_next
        x = x_next
        y = y_next
    end

    length(word) == 2 * target || throw(ArgumentError("approximate sphere word has the wrong length"))
    return String(Char.(word))
end

function _mated_crt_edge_kind_from_fk_color(color::Integer)
    if Int(color) == Int(RED)
        return :upper
    elseif Int(color) == Int(BLUE)
        return :lower
    elseif Int(color) == Int(GREEN)
        return :spine
    end
    return :generic
end

function _mated_crt_triangle_edge_ids_from_fk(map_data::FKMap)
    tri_faces = Int32.(map_data.triangulation_faces)
    tri_green = Int32.(map_data.triangle_green_edges)
    tri_diag = Int32.(map_data.triangle_diag_edge)
    nrows = size(tri_faces, 1)

    size(tri_green, 1) == nrows || throw(ArgumentError("triangle_green_edges must align with triangulation_faces"))
    length(tri_diag) == nrows || throw(ArgumentError("triangle_diag_edge must align with triangulation_faces"))

    face_edge_ids = Vector{Vector{Int32}}(undef, nrows)
    for i in 1:nrows
        a = tri_faces[i, 1]
        b = tri_faces[i, 2]
        c = tri_faces[i, 3]
        side_pairs = (
            (min(a, b), max(a, b)),
            (min(b, c), max(b, c)),
            (min(c, a), max(c, a)),
        )
        triangle_edges = Int32[tri_green[i, 1], tri_green[i, 2], tri_diag[i]]
        assigned = falses(3)
        side_ids = fill(Int32(-1), 3)

        for edge_id in triangle_edges
            edge_id >= 0 || throw(ArgumentError("triangle edge ids must be nonnegative"))
            pair = (
                min(map_data.edge_u[Int(edge_id) + 1], map_data.edge_v[Int(edge_id) + 1]),
                max(map_data.edge_u[Int(edge_id) + 1], map_data.edge_v[Int(edge_id) + 1]),
            )

            matched = 0
            for j in 1:3
                if !assigned[j] && side_pairs[j] == pair
                    matched = j
                    break
                end
            end
            matched != 0 || throw(ArgumentError("failed to match a spanning-tree triangle edge to a triangle side"))
            side_ids[matched] = edge_id
            assigned[matched] = true
        end

        all(assigned) || throw(ArgumentError("incomplete spanning-tree triangle edge-id assignment"))
        face_edge_ids[i] = side_ids
    end

    return face_edge_ids
end

function _build_mated_crt_sphere_map_from_word(
    params::Dict{String,Float64},
    word::AbstractString;
    sampler::AbstractString="exact_correlated_cone_walk",
)
    st_map = build_fk_map_from_word(word; require_excursion=true)
    faces = [Int32.(st_map.triangulation_faces[i, :]) for i in 1:size(st_map.triangulation_faces, 1)]
    face_edge_ids = _mated_crt_triangle_edge_ids_from_fk(st_map)
    edge_kind = [_mated_crt_edge_kind_from_fk_color(st_map.edge_color[i]) for i in eachindex(st_map.edge_color)]
    L, R = _mated_crt_exact_word_coordinates(word)
    time_grid = collect(range(0.0, 1.0; length=length(L)))

    return MatedCRTMap(
        :sphere,
        num_vertices(st_map),
        params["gamma"],
        params["gamma_prime"],
        params["kappa"],
        params["kappa_prime"],
        params["correlation"],
        (0.0, 0.0),
        1,
        String(sampler),
        time_grid,
        L,
        R,
        faces,
        face_edge_ids,
        Int32[],
        Int32(-1),
        copy(st_map.edge_u),
        copy(st_map.edge_v),
        edge_kind,
        st_map,
    )
end

function _build_mated_crt_exact_sphere_map(
    params::Dict{String,Float64},
    target_vertices::Integer,
    rng::AbstractRNG,
)
    n = Int(target_vertices)
    n >= 3 || throw(ArgumentError("mated_crt sphere topology requires at least 3 vertices"))
    target_productions = n - 2

    sampler = if abs(params["correlation"]) <= 1.0e-12
        "exact_spanning_tree_word"
    else
        "exact_lazy_correlated_cone_walk"
    end
    word = _mated_crt_sample_exact_sphere_word(target_productions, params["correlation"], rng)
    return _build_mated_crt_sphere_map_from_word(
        params,
        word;
        sampler=sampler * "+spanning_tree_word",
    )
end

function _build_mated_crt_approx_sphere_map(
    params::Dict{String,Float64},
    target_vertices::Integer,
    rng::AbstractRNG;
    exact_tail_cutoff::Integer=64,
)
    n = Int(target_vertices)
    n >= 3 || throw(ArgumentError("mated_crt sphere topology requires at least 3 vertices"))
    target_productions = n - 2
    tail_cutoff = max(0, Int(exact_tail_cutoff))

    sampler = if abs(params["correlation"]) <= 1.0e-12
        "approx_hybrid_zero_corr_exact_tail$(tail_cutoff)"
    else
        "approx_hybrid_cone_walk_tail$(tail_cutoff)"
    end
    word = _mated_crt_sample_approx_sphere_word(
        target_productions,
        params["correlation"],
        rng;
        exact_tail_cutoff=tail_cutoff,
    )
    return _build_mated_crt_sphere_map_from_word(
        params,
        word;
        sampler=sampler * "+spanning_tree_word",
    )
end

function _mated_crt_face_cycles(
    num_vertices::Integer,
    edge_u::Vector{Int32},
    edge_v::Vector{Int32},
    edge_kind::Vector{Symbol},
    boundary_vertices::Vector{Int32}=Int32[],
)
    n = Int(num_vertices)
    medges = length(edge_u)
    length(edge_v) == medges || throw(ArgumentError("edge_v must align with edge_u"))
    length(edge_kind) == medges || throw(ArgumentError("edge_kind must align with edge_u"))

    half_tail = Vector{Int32}(undef, 2 * medges)
    half_head = Vector{Int32}(undef, 2 * medges)
    half_edge = Vector{Int32}(undef, 2 * medges)
    incident = [Int[] for _ in 1:n]

    for i in 1:medges
        u = Int(edge_u[i])
        v = Int(edge_v[i])
        (0 <= u < n && 0 <= v < n) || throw(ArgumentError("edge endpoints must lie in [0, n)"))
        h_uv = 2 * i - 1
        h_vu = 2 * i
        half_tail[h_uv] = Int32(u)
        half_head[h_uv] = Int32(v)
        half_tail[h_vu] = Int32(v)
        half_head[h_vu] = Int32(u)
        half_edge[h_uv] = Int32(i - 1)
        half_edge[h_vu] = Int32(i - 1)
        push!(incident[u + 1], h_uv)
        push!(incident[v + 1], h_vu)
    end

    boundary_next = Dict{Int32,Int32}()
    boundary_prev = Dict{Int32,Int32}()
    boundary_edge_ids = Dict{Tuple{Int32,Int32},Int32}()

    if !isempty(boundary_vertices)
        edge_lookup = Dict{Tuple{Int32,Int32},Vector{Int32}}()
        for i in eachindex(edge_u)
            pair = (
                min(edge_u[i], edge_v[i]),
                max(edge_u[i], edge_v[i]),
            )
            push!(get!(edge_lookup, pair, Int32[]), Int32(i - 1))
        end

        m = length(boundary_vertices)
        for i in 1:m
            a = boundary_vertices[i]
            b = boundary_vertices[i == m ? 1 : i + 1]
            boundary_next[a] = b
            boundary_prev[b] = a

            pair = (min(a, b), max(a, b))
            candidates = get(edge_lookup, pair, Int32[])
            !isempty(candidates) || throw(ArgumentError("failed to identify a mated_crt boundary edge between vertices $a and $b"))

            chosen = Int32(-1)
            for eid in candidates
                if edge_kind[Int(eid) + 1] == :lower
                    chosen = eid
                    break
                end
            end
            if chosen < 0
                for eid in candidates
                    if edge_kind[Int(eid) + 1] == :spine
                        chosen = eid
                        break
                    end
                end
            end
            if chosen < 0
                chosen = candidates[1]
            end

            boundary_edge_ids[(a, b)] = chosen
            boundary_edge_ids[(b, a)] = chosen
        end
    end

    function halfedge_order_key(h::Int)
        tail = half_tail[h]
        head = half_head[h]
        edge_id = half_edge[h]
        if !isempty(boundary_edge_ids) && get(boundary_edge_ids, (tail, head), Int32(-1)) == edge_id
            if get(boundary_next, tail, Int32(-1)) == head
                return (0, 0)
            elseif get(boundary_prev, tail, Int32(-1)) == head
                return (7, 0)
            end
        end

        edge_idx = fld(h + 1, 2)
        kind = edge_kind[edge_idx]
        u = Int(tail)
        v = Int(head)
        span = abs(v - u)

        if kind == :spine && v > u
            return (1, 0)
        elseif kind == :upper && v > u
            return (2, span)
        elseif kind == :upper && v < u
            return (3, -span)
        elseif kind == :spine && v < u
            return (4, 0)
        elseif kind == :lower && v < u
            return (5, span)
        elseif kind == :lower && v > u
            return (6, -span)
        end
        return (8, span)
    end

    for v in eachindex(incident)
        sort!(incident[v]; by=halfedge_order_key)
    end

    position_in_rotation = zeros(Int, 2 * medges)
    for v in eachindex(incident)
        for (idx, h) in enumerate(incident[v])
            position_in_rotation[h] = idx
        end
    end

    twin(h::Int) = isodd(h) ? h + 1 : h - 1
    next_half = zeros(Int, 2 * medges)
    for h in 1:(2 * medges)
        head = Int(half_head[h]) + 1
        rot = incident[head]
        isempty(rot) && throw(ArgumentError("rotation system may not have isolated vertices"))
        t = twin(h)
        idx = position_in_rotation[t]
        prev_idx = idx == 1 ? length(rot) : idx - 1
        next_half[h] = rot[prev_idx]
    end

    seen = falses(2 * medges)
    faces = Vector{Vector{Int32}}()
    face_edge_ids = Vector{Vector{Int32}}()

    for h0 in 1:(2 * medges)
        seen[h0] && continue
        cycle_vertices = Int32[]
        cycle_edges = Int32[]
        h = h0
        while !seen[h]
            seen[h] = true
            push!(cycle_vertices, half_tail[h])
            push!(cycle_edges, half_edge[h])
            h = next_half[h]
        end
        isempty(cycle_vertices) && continue
        push!(faces, cycle_vertices)
        push!(face_edge_ids, cycle_edges)
    end

    return faces, face_edge_ids
end

function _mated_crt_split_pinched_face_walk(face::Vector{Int32}, edge_ids::Vector{Int32})
    length(face) == length(edge_ids) ||
        throw(ArgumentError("face edge ids must align with the face walk"))

    pending = Tuple{Vector{Int32},Vector{Int32}}[(copy(face), copy(edge_ids))]
    simple_faces = Vector{Vector{Int32}}()
    simple_face_edge_ids = Vector{Vector{Int32}}()

    while !isempty(pending)
        current_face, current_ids = pop!(pending)
        positions = Dict{Int32,Int}()
        split_face = false

        for j in eachindex(current_face)
            v = current_face[j]
            if haskey(positions, v)
                i = positions[v]
                if j - i >= 2
                    first_piece = current_face[i:(j - 1)]
                    first_ids = current_ids[i:(j - 1)]
                    second_piece = vcat(current_face[1:(i - 1)], current_face[j:end])
                    second_ids = vcat(current_ids[1:(i - 1)], current_ids[j:end])

                    length(first_piece) >= 2 && push!(pending, (first_piece, first_ids))
                    length(second_piece) >= 2 && push!(pending, (second_piece, second_ids))
                    split_face = true
                    break
                end
            else
                positions[v] = j
            end
        end

        split_face && continue
        push!(simple_faces, current_face)
        push!(simple_face_edge_ids, current_ids)
    end

    return simple_faces, simple_face_edge_ids
end

function _collapse_mated_crt_parallel_edges(
    edge_u::Vector{Int32},
    edge_v::Vector{Int32},
    edge_kind::Vector{Symbol},
    faces::Vector{Vector{Int32}},
    face_edge_ids::Vector{Vector{Int32}},
)
    length(edge_v) == length(edge_u) || throw(ArgumentError("edge_v must align with edge_u"))
    length(edge_kind) == length(edge_u) || throw(ArgumentError("edge_kind must align with edge_u"))
    length(face_edge_ids) == length(faces) || throw(ArgumentError("face_edge_ids must align with faces"))

    pair_to_old_ids = Dict{Tuple{Int32,Int32},Vector{Int32}}()
    for i in eachindex(edge_u)
        pair = (
            min(edge_u[i], edge_v[i]),
            max(edge_u[i], edge_v[i]),
        )
        push!(get!(pair_to_old_ids, pair, Int32[]), Int32(i - 1))
    end

    sorted_pairs = sort!(collect(keys(pair_to_old_ids)); by=pair -> minimum(pair_to_old_ids[pair]))
    old_to_new = Dict{Int32,Int32}()
    new_edge_u = Int32[]
    new_edge_v = Int32[]
    new_edge_kind = Symbol[]

    for (new_idx, pair) in enumerate(sorted_pairs)
        push!(new_edge_u, pair[1])
        push!(new_edge_v, pair[2])

        kinds = unique(Symbol[edge_kind[Int(eid) + 1] for eid in pair_to_old_ids[pair]])
        chosen_kind = if length(kinds) == 1
            kinds[1]
        elseif :upper in kinds && :lower in kinds
            :merged
        elseif :spine in kinds
            :spine
        elseif :upper in kinds
            :upper
        elseif :lower in kinds
            :lower
        else
            kinds[1]
        end
        push!(new_edge_kind, chosen_kind)

        for old_id in pair_to_old_ids[pair]
            old_to_new[old_id] = Int32(new_idx - 1)
        end
    end

    new_faces = Vector{Vector{Int32}}()
    new_face_edge_ids = Vector{Vector{Int32}}()
    for (face, ids) in zip(faces, face_edge_ids)
        length(face) == length(ids) || throw(ArgumentError("face edge ids must align with the face vertices"))
        length(face) == 2 && continue
        push!(new_faces, Int32.(face))
        push!(new_face_edge_ids, Int32[old_to_new[eid] for eid in ids])
    end

    return new_edge_u, new_edge_v, new_edge_kind, new_faces, new_face_edge_ids
end

function _mated_crt_outermost_intervals(
    topology::Symbol,
    raw_num_vertices::Integer,
    edge_u::Vector{Int32},
    edge_v::Vector{Int32},
    edge_kind::Vector{Symbol},
)
    n = Int(raw_num_vertices)
    intervals = Tuple{Int32,Int32}[]
    for i in eachindex(edge_u)
        edge_kind[i] == :merged || continue
        a = min(edge_u[i], edge_v[i])
        b = max(edge_u[i], edge_v[i])
        Int(b - a) > 1 || continue
        topology == :sphere && a == 0 && b == n - 1 && continue
        push!(intervals, (a, b))
    end
    intervals = unique(sort(intervals))

    outer = Tuple{Int32,Int32}[]
    for iv in intervals
        contained = false
        for jv in intervals
            iv == jv && continue
            if jv[1] <= iv[1] && iv[2] <= jv[2] && (jv[1] < iv[1] || iv[2] < jv[2])
                contained = true
                break
            end
        end
        contained || push!(outer, iv)
    end
    return unique(sort(outer))
end

function _mated_crt_edge_kind_from_candidates(candidates::Vector{Symbol})
    kinds = Set(candidates)
    if :spine in kinds
        return :spine
    elseif :upper in kinds && :lower in kinds
        return :merged
    elseif :merged in kinds
        return :merged
    elseif :upper in kinds
        return :upper
    elseif :lower in kinds
        return :lower
    end
    return isempty(candidates) ? :generic : candidates[1]
end

_mated_crt_triangle_key(a::Int32, b::Int32, c::Int32) =
    Tuple(sort(Int32[a, b, c]))

function _mated_crt_triangle_pairs(a::Int32, b::Int32, c::Int32)
    return (
        (min(a, b), max(a, b)),
        (min(b, c), max(b, c)),
        (min(c, a), max(c, a)),
    )
end

function _mated_crt_collapsed_face_components(
    discarded_face_vertices::Vector{Vector{Int32}},
    discarded_removed_vertices::Vector{Vector{Int32}},
)
    length(discarded_removed_vertices) == length(discarded_face_vertices) ||
        throw(ArgumentError("discarded mated_crt face metadata must align"))
    isempty(discarded_face_vertices) && return Vector{Vector{Int32}}()

    removed_to_faces = Dict{Int32,Vector{Int}}()
    for (idx, removed_vertices) in enumerate(discarded_removed_vertices)
        for v in removed_vertices
            push!(get!(removed_to_faces, v, Int[]), idx)
        end
    end

    visited = falses(length(discarded_face_vertices))
    components = Vector{Vector{Int32}}()
    for start in eachindex(discarded_face_vertices)
        visited[start] && continue
        queue = Int[start]
        visited[start] = true
        boundary_vertices = Set(discarded_face_vertices[start])
        while !isempty(queue)
            idx = pop!(queue)
            union!(boundary_vertices, discarded_face_vertices[idx])
            for v in discarded_removed_vertices[idx]
                for nbr in get(removed_to_faces, v, Int[])
                    visited[nbr] && continue
                    visited[nbr] = true
                    push!(queue, nbr)
                end
            end
        end
        push!(components, sort!(collect(boundary_vertices)))
    end
    return components
end

function _mated_crt_candidate_missing_triangles(
    topology::Symbol,
    collapsed_components::Vector{Vector{Int32}},
    available_pairs::Set{Tuple{Int32,Int32}},
    triangle_keys::Set{NTuple{3,Int32}},
    boundary_vertices::Vector{Int32},
    edge_slack::Dict{Tuple{Int32,Int32},Int},
)
    isempty(collapsed_components) && return NTuple{3,Int32}[]

    outer_key = if topology == :disk && length(boundary_vertices) == 3
        _mated_crt_triangle_key(boundary_vertices[1], boundary_vertices[2], boundary_vertices[3])
    else
        nothing
    end

    candidates = NTuple{3,Int32}[]
    for component in collapsed_components
        length(component) >= 3 || continue
        for ia in 1:(length(component) - 2)
            for ib in (ia + 1):(length(component) - 1)
                for ic in (ib + 1):length(component)
                    key = (component[ia], component[ib], component[ic])
                    key in triangle_keys && continue
                    outer_key !== nothing && key == outer_key && continue

                    pairs = _mated_crt_triangle_pairs(key...)
                    all(pair -> pair in available_pairs, pairs) || continue
                    all(get(edge_slack, pair, 0) > 0 for pair in pairs) || continue
                    push!(candidates, key)
                end
            end
        end
    end
    return unique(sort(candidates))
end

function _solve_mated_crt_missing_triangles(
    edge_slack::Dict{Tuple{Int32,Int32},Int},
    candidates::Vector{NTuple{3,Int32}},
    missing_faces::Integer,
)
    missing = Int(missing_faces)
    missing >= 0 || throw(ArgumentError("missing_faces must be nonnegative"))
    missing == 0 && return NTuple{3,Int32}[]

    working_slack = copy(edge_slack)
    chosen = Int[]

    candidate_pairs = [_mated_crt_triangle_pairs(tri...) for tri in candidates]

    function recurse(start_idx::Int, remaining::Int)
        remaining == 0 && return true
        length(candidates) - start_idx + 1 >= remaining || return false

        for idx in start_idx:length(candidates)
            all(get(working_slack, pair, 0) > 0 for pair in candidate_pairs[idx]) || continue
            valid = true
            for pair in candidate_pairs[idx]
                working_slack[pair] = get(working_slack, pair, 0) - 1
                if working_slack[pair] < 0
                    valid = false
                end
            end
            if valid
                push!(chosen, idx)
                recurse(idx + 1, remaining - 1) && return true
                pop!(chosen)
            end
            for pair in candidate_pairs[idx]
                working_slack[pair] = get(working_slack, pair, 0) + 1
            end
        end
        return false
    end

    recurse(1, missing) || return nothing
    return [candidates[idx] for idx in chosen]
end

function _recover_mated_crt_missing_triangles!(
    topology::Symbol,
    internal_faces::Vector{Vector{Int32}},
    triangle_keys::Set{NTuple{3,Int32}},
    available_pairs::Set{Tuple{Int32,Int32}},
    boundary_vertices::Vector{Int32},
    collapsed_components::Vector{Vector{Int32}},
    expected_internal_faces::Integer,
)
    target_faces = Int(expected_internal_faces)
    target_faces >= length(internal_faces) ||
        throw(ArgumentError("collapsed mated_crt quotient produced too many internal faces"))

    missing_faces = target_faces - length(internal_faces)
    missing_faces == 0 && return

    boundary_pairs = Set{Tuple{Int32,Int32}}()
    if topology == :disk
        for i in eachindex(boundary_vertices)
            j = i == length(boundary_vertices) ? 1 : i + 1
            push!(boundary_pairs, (min(boundary_vertices[i], boundary_vertices[j]), max(boundary_vertices[i], boundary_vertices[j])))
        end
    end

    edge_counts = Dict{Tuple{Int32,Int32},Int}()
    for face in internal_faces
        a, b, c = face
        for pair in _mated_crt_triangle_pairs(a, b, c)
            edge_counts[pair] = get(edge_counts, pair, 0) + 1
        end
    end

    edge_slack = Dict{Tuple{Int32,Int32},Int}()
    for pair in available_pairs
        target = topology == :disk && pair in boundary_pairs ? 1 : 2
        slack = target - get(edge_counts, pair, 0)
        slack < 0 && throw(ArgumentError("collapsed mated_crt quotient overuses an edge"))
        slack > 0 && (edge_slack[pair] = slack)
    end

    candidates = _mated_crt_candidate_missing_triangles(
        topology,
        collapsed_components,
        available_pairs,
        triangle_keys,
        boundary_vertices,
        edge_slack,
    )
    solution = _solve_mated_crt_missing_triangles(edge_slack, candidates, missing_faces)
    solution === nothing &&
        throw(ArgumentError("failed to recover the missing collapsed mated_crt quotient triangles"))

    for tri in solution
        push!(internal_faces, Int32[tri[1], tri[2], tri[3]])
        push!(triangle_keys, tri)
    end
end

function _validate_mated_crt_internal_edge_counts(
    topology::Symbol,
    internal_faces::Vector{Vector{Int32}},
    boundary_vertices::Vector{Int32},
)
    edge_counts = Dict{Tuple{Int32,Int32},Int}()
    for face in internal_faces
        a, b, c = face
        for pair in _mated_crt_triangle_pairs(a, b, c)
            edge_counts[pair] = get(edge_counts, pair, 0) + 1
        end
    end

    boundary_pairs = Set{Tuple{Int32,Int32}}()
    if topology == :disk
        for i in eachindex(boundary_vertices)
            j = i == length(boundary_vertices) ? 1 : i + 1
            pair = (min(boundary_vertices[i], boundary_vertices[j]), max(boundary_vertices[i], boundary_vertices[j]))
            push!(boundary_pairs, pair)
            edge_counts[pair] = get(edge_counts, pair, 0) + 1
        end
    end

    for (pair, count) in edge_counts
        expected = topology == :disk && pair in boundary_pairs ? 2 : 2
        count == expected ||
            throw(ArgumentError("collapsed mated_crt quotient edge incidences do not form a triangulation"))
    end
end

function _mated_crt_collapse_intervals(
    topology::Symbol,
    raw_num_vertices::Integer,
    boundary_vertices::Vector{Int32},
    outer_face_index::Int32,
    faces::Vector{Vector{Int32}},
    face_edge_ids::Vector{Vector{Int32}},
    edge_u::Vector{Int32},
    edge_v::Vector{Int32},
    edge_kind::Vector{Symbol},
)
    n_raw = Int(raw_num_vertices)
    intervals = _mated_crt_outermost_intervals(topology, n_raw, edge_u, edge_v, edge_kind)
    isempty(intervals) && return (
        nverts=n_raw,
        boundary_vertices=copy(boundary_vertices),
        outer_face_index=outer_face_index,
        faces=copy(faces),
        face_edge_ids=copy([Int32.(ids) for ids in face_edge_ids]),
        edge_u=copy(edge_u),
        edge_v=copy(edge_v),
        edge_kind=copy(edge_kind),
    )

    removed = falses(n_raw)
    for (a, b) in intervals
        for v in (Int(a) + 1):(Int(b) - 1)
            removed[v + 1] = true
        end
    end

    keep_vertices = Int32[v for v in 0:(n_raw - 1) if !removed[v + 1]]
    old_to_new = Dict{Int,Int}(Int(v) => i - 1 for (i, v) in enumerate(keep_vertices))

    pair_to_rank = Dict{Tuple{Int32,Int32},Int}()
    pair_to_kinds = Dict{Tuple{Int32,Int32},Vector{Symbol}}()
    for i in eachindex(edge_u)
        if removed[Int(edge_u[i]) + 1] || removed[Int(edge_v[i]) + 1]
            continue
        end
        pair = (min(edge_u[i], edge_v[i]), max(edge_u[i], edge_v[i]))
        if !haskey(pair_to_rank, pair)
            pair_to_rank[pair] = i - 1
            pair_to_kinds[pair] = Symbol[]
        else
            pair_to_rank[pair] = min(pair_to_rank[pair], i - 1)
        end
        push!(pair_to_kinds[pair], edge_kind[i])
    end

    boundary_new = Int32[]
    if topology == :disk
        for v in boundary_vertices
            haskey(old_to_new, Int(v)) || continue
            mapped = Int32(old_to_new[Int(v)])
            isempty(boundary_new) || boundary_new[end] != mapped || continue
            push!(boundary_new, mapped)
        end
        length(boundary_new) >= 2 && boundary_new[1] == boundary_new[end] && pop!(boundary_new)
    end

    mapped_pair_to_rank = Dict{Tuple{Int32,Int32},Int}()
    mapped_pair_to_kinds = Dict{Tuple{Int32,Int32},Vector{Symbol}}()
    for (pair, rank) in pair_to_rank
        mapped_pair = (
            Int32(old_to_new[Int(pair[1])]),
            Int32(old_to_new[Int(pair[2])]),
        )
        mapped_pair = (min(mapped_pair[1], mapped_pair[2]), max(mapped_pair[1], mapped_pair[2]))
        mapped_pair_to_rank[mapped_pair] = min(get(mapped_pair_to_rank, mapped_pair, typemax(Int)), rank)
        append!(get!(mapped_pair_to_kinds, mapped_pair, Symbol[]), get(pair_to_kinds, pair, Symbol[]))
    end

    internal_faces = Vector{Vector{Int32}}()
    seen_face_keys = Set{NTuple{3,Int32}}()
    discarded_face_vertices = Vector{Vector{Int32}}()
    discarded_removed_vertices = Vector{Vector{Int32}}()
    outer_old = Int(outer_face_index)
    for i in eachindex(faces)
        topology == :disk && outer_old >= 0 && i == outer_old + 1 && continue
        mapped = Int32[]
        for v in faces[i]
            haskey(old_to_new, Int(v)) || continue
            push!(mapped, Int32(old_to_new[Int(v)]))
        end
        mapped_unique = unique(mapped)
        if length(mapped_unique) != 3
            push!(discarded_face_vertices, sort(mapped_unique))
            push!(discarded_removed_vertices, unique(Int32[v for v in faces[i] if !haskey(old_to_new, Int(v))]))
            continue
        end
        key = Tuple(Int32.(sort(mapped)))
        key in seen_face_keys && continue
        push!(seen_face_keys, key)
        push!(internal_faces, mapped)
    end

    expected_internal_faces = if topology == :disk
        2 * length(keep_vertices) - 2 - length(boundary_new)
    else
        2 * length(keep_vertices) - 4
    end

    available_pairs = Set(keys(mapped_pair_to_rank))
    collapsed_components = _mated_crt_collapsed_face_components(
        discarded_face_vertices,
        discarded_removed_vertices,
    )
    _recover_mated_crt_missing_triangles!(
        topology,
        internal_faces,
        seen_face_keys,
        available_pairs,
        boundary_new,
        collapsed_components,
        expected_internal_faces,
    )
    _validate_mated_crt_internal_edge_counts(topology, internal_faces, boundary_new)

    face_edge_pairs = Vector{Vector{Tuple{Int32,Int32}}}()
    pair_set = Set{Tuple{Int32,Int32}}()
    for face in internal_faces
        a, b, c = face
        pairs = [
            (min(a, b), max(a, b)),
            (min(b, c), max(b, c)),
            (min(c, a), max(c, a)),
        ]
        push!(face_edge_pairs, pairs)
        union!(pair_set, pairs)
    end

    outer_pairs = Tuple{Int32,Int32}[]
    if topology == :disk
        for i in eachindex(boundary_new)
            j = i == length(boundary_new) ? 1 : i + 1
            pair = (min(boundary_new[i], boundary_new[j]), max(boundary_new[i], boundary_new[j]))
            push!(outer_pairs, pair)
            push!(pair_set, pair)
        end
    end

    ordered_pairs = sort!(collect(pair_set); by=pair -> get(mapped_pair_to_rank, pair, typemax(Int)))
    pair_to_new_id = Dict{Tuple{Int32,Int32},Int32}()
    new_edge_u = Int32[]
    new_edge_v = Int32[]
    new_edge_kind = Symbol[]
    for (i, pair) in enumerate(ordered_pairs)
        pair_to_new_id[pair] = Int32(i - 1)
        push!(new_edge_u, pair[1])
        push!(new_edge_v, pair[2])
        push!(new_edge_kind, _mated_crt_edge_kind_from_candidates(get(mapped_pair_to_kinds, pair, Symbol[])))
    end

    final_faces = Vector{Vector{Int32}}()
    final_face_edge_ids = Vector{Vector{Int32}}()
    for (face, pairs) in zip(internal_faces, face_edge_pairs)
        push!(final_faces, face)
        push!(final_face_edge_ids, Int32[pair_to_new_id[p] for p in pairs])
    end

    new_outer_face_index = Int32(-1)
    if topology == :disk
        push!(final_faces, copy(boundary_new))
        push!(final_face_edge_ids, Int32[pair_to_new_id[p] for p in outer_pairs])
        new_outer_face_index = Int32(length(final_faces) - 1)
    end

    n_new = length(keep_vertices)
    if topology == :disk
        b = length(boundary_new)
        expected_faces = 2 * n_new - 2 - b
        expected_edges = 3 * n_new - 3 - b
        length(boundary_new) >= 3 || throw(ArgumentError("collapsed mated_crt disk boundary must have length at least 3"))
        length(final_faces) == expected_faces + 1 || throw(ArgumentError("collapsed mated_crt disk faces do not form a triangulation"))
        length(new_edge_u) == expected_edges || throw(ArgumentError("collapsed mated_crt disk edges do not form a triangulation"))
    else
        expected_faces = 2 * n_new - 4
        expected_edges = 3 * n_new - 6
        length(final_faces) == expected_faces || throw(ArgumentError("collapsed mated_crt sphere faces do not form a triangulation"))
        length(new_edge_u) == expected_edges || throw(ArgumentError("collapsed mated_crt sphere edges do not form a triangulation"))
    end

    return (
        nverts=n_new,
        boundary_vertices=boundary_new,
        outer_face_index=new_outer_face_index,
        faces=final_faces,
        face_edge_ids=final_face_edge_ids,
        edge_u=new_edge_u,
        edge_v=new_edge_v,
        edge_kind=new_edge_kind,
    )
end

function _same_cyclic_cycle(a, b)
    length(a) == length(b) || return false
    isempty(a) && return true
    aa = Int32.(collect(a))
    bb = Int32.(collect(b))
    n = length(aa)

    function _match(candidate)
        doubled = vcat(candidate, candidate)
        for shift in 1:n
            ok = true
            for i in 1:n
                if doubled[shift + i - 1] != bb[i]
                    ok = false
                    break
                end
            end
            ok && return true
        end
        return false
    end

    _match(aa) && return true
    return _match(reverse(aa))
end

function _matrix_from_face_list(faces::Vector{Vector{Int32}})
    isempty(faces) && return Matrix{Int32}(undef, 0, 3)
    width = maximum(length(face) for face in faces)
    out = fill(Int32(-1), length(faces), width)
    for i in eachindex(faces)
        for j in eachindex(faces[i])
            out[i, j] = faces[i][j]
        end
    end
    return out
end

function _mated_crt_boundary_vertices(L::Vector{Float64}, refinement::Integer, nverts::Integer)
    r = Int(refinement)
    n = Int(nverts)
    suffix = zeros(Float64, n)
    current = L[n * r + 1]
    suffix[end] = current
    for p in (n - 1):-1:1
        current = min(current, minimum(@view L[p * r + 1:(p + 1) * r + 1]))
        suffix[p] = current
    end

    boundary = Int32[]
    for p in 1:n
        cell_min = minimum(@view L[(p - 1) * r + 1:p * r + 1])
        if cell_min <= suffix[p] + 1.0e-10
            push!(boundary, Int32(p - 1))
        end
    end
    return boundary
end

function _build_mated_crt_edges(
    L::Vector{Float64},
    R::Vector{Float64},
    nverts::Integer,
    refinement::Integer,
)
    n = Int(nverts)
    r = Int(refinement)

    cell_min_L = [minimum(@view L[(p - 1) * r + 1:p * r + 1]) for p in 1:n]
    cell_min_R = [minimum(@view R[(p - 1) * r + 1:p * r + 1]) for p in 1:n]

    edge_u = Int32[]
    edge_v = Int32[]
    edge_kind = Symbol[]

    for p in 1:(n - 1)
        push!(edge_u, Int32(p - 1))
        push!(edge_v, Int32(p))
        push!(edge_kind, :spine)
    end

    for p in 1:(n - 1)
        mid_min_L = L[p * r + 1]
        mid_min_R = R[p * r + 1]
        for q in (p + 1):n
            if q > p + 1
                mid_min_L = min(mid_min_L, cell_min_L[q - 1])
                mid_min_R = min(mid_min_R, cell_min_R[q - 1])
            end

            if q > p + 1
                if max(cell_min_L[p], cell_min_L[q]) <= mid_min_L + 1.0e-10
                    push!(edge_u, Int32(p - 1))
                    push!(edge_v, Int32(q - 1))
                    push!(edge_kind, :lower)
                end
                if max(cell_min_R[p], cell_min_R[q]) <= mid_min_R + 1.0e-10
                    push!(edge_u, Int32(p - 1))
                    push!(edge_v, Int32(q - 1))
                    push!(edge_kind, :upper)
                end
            end
        end
    end

    return edge_u, edge_v, edge_kind
end

function _build_mated_crt_map_from_paths(
    topology::Symbol,
    params::Dict{String,Float64},
    L::Vector{Float64},
    R::Vector{Float64},
    refinement::Integer;
    sampler::AbstractString="discrete_positive_bridge_gibbs",
)
    length(L) == length(R) || throw(ArgumentError("L and R must have the same length"))
    r = Int(refinement)
    r >= 1 || throw(ArgumentError("refinement must be >= 1"))
    m = length(L) - 1
    m >= r || throw(ArgumentError("path must contain at least one vertex interval"))
    n = fld(m, r)
    n >= 3 || throw(ArgumentError("mated_crt requires at least 3 vertices"))
    n * r == m || throw(ArgumentError("path length must be divisible by refinement"))

    boundary_vertices = topology == :disk ? _mated_crt_boundary_vertices(L, r, n) : Int32[]
    edge_u, edge_v, edge_kind = _build_mated_crt_edges(L, R, n, r)
    raw_faces, raw_face_edge_ids = _mated_crt_face_cycles(n, edge_u, edge_v, edge_kind, boundary_vertices)

    faces = Vector{Vector{Int32}}()
    face_edge_ids = Vector{Vector{Int32}}()
    for (face, ids) in zip(raw_faces, raw_face_edge_ids)
        split_faces, split_ids = _mated_crt_split_pinched_face_walk(face, ids)
        append!(faces, split_faces)
        append!(face_edge_ids, split_ids)
    end

    edge_u, edge_v, edge_kind, faces, face_edge_ids =
        _collapse_mated_crt_parallel_edges(edge_u, edge_v, edge_kind, faces, face_edge_ids)

    outer_face_index = Int32(-1)

    if topology == :disk
        for i in eachindex(faces)
            if _same_cyclic_cycle(faces[i], boundary_vertices)
                outer_face_index = Int32(i - 1)
                break
            end
        end
        outer_face_index >= 0 || throw(ArgumentError("failed to identify the mated_crt outer face"))

        outer_boundary = faces[Int(outer_face_index) + 1]
        Set(outer_boundary) == Set(boundary_vertices) ||
            throw(ArgumentError("paper boundary definition does not match the recovered outer face"))
    end

    for i in eachindex(faces)
        topology == :disk && i == Int(outer_face_index) + 1 && continue
        length(faces[i]) == 3 || throw(ArgumentError("mated_crt internal faces must be triangular after collapsing digons"))
    end

    quotient = _mated_crt_collapse_intervals(
        topology,
        n,
        boundary_vertices,
        outer_face_index,
        faces,
        face_edge_ids,
        edge_u,
        edge_v,
        edge_kind,
    )

    time_grid = collect(range(0.0, 1.0; length=length(L)))
    endpoint = (L[end], R[end])

    return MatedCRTMap(
        topology,
        quotient.nverts,
        params["gamma"],
        params["gamma_prime"],
        params["kappa"],
        params["kappa_prime"],
        params["correlation"],
        endpoint,
        r,
        String(sampler) * "+interval_quotient",
        time_grid,
        copy(L),
        copy(R),
        quotient.faces,
        quotient.face_edge_ids,
        quotient.boundary_vertices,
        quotient.outer_face_index,
        quotient.edge_u,
        quotient.edge_v,
        quotient.edge_kind,
        nothing,
    )
end

function build_mated_crt_map(;
    vertices::Integer,
    topology::Union{Symbol,AbstractString}="disk",
    gamma=nothing,
    gamma_prime=nothing,
    kappa=nothing,
    kappa_prime=nothing,
    correlation=nothing,
    seed::Integer=1,
    refinement::Integer=4,
    burnin_sweeps::Integer=64,
    gibbs_sweeps::Integer=128,
    max_tries::Integer=12,
    sphere_sampler::Union{Symbol,AbstractString}=:exact,
    sphere_exact_tail_cutoff::Integer=64,
)
    n = Int(vertices)
    n >= 3 || throw(ArgumentError("vertices must be >= 3"))
    topo = _mated_crt_topology(topology)
    r = Int(refinement)
    r >= 1 || throw(ArgumentError("refinement must be >= 1"))
    attempts = max(1, Int(max_tries))
    params = _resolve_mated_crt_parameters(
        ;
        gamma=gamma,
        gamma_prime=gamma_prime,
        kappa=kappa,
        kappa_prime=kappa_prime,
        correlation=correlation,
    )

    if topo == :sphere
        rng = MersenneTwister(Int(seed))
        sphere_mode = _mated_crt_sphere_sampler(sphere_sampler)
        if sphere_mode == :exact
            return _build_mated_crt_exact_sphere_map(params, n, rng)
        end
        return _build_mated_crt_approx_sphere_map(
            params,
            n,
            rng;
            exact_tail_cutoff=sphere_exact_tail_cutoff,
        )
    end

    endpoint = (1.0, 0.0)
    total_steps = n * r

    last_error = nothing
    for attempt in 1:attempts
        rng = MersenneTwister(Int(seed) + attempt - 1)
        path = _sample_positive_bridge_path(
            endpoint,
            params["correlation"],
            total_steps,
            rng;
            burnin_sweeps=burnin_sweeps,
            gibbs_sweeps=gibbs_sweeps,
        )
        try
            return _build_mated_crt_map_from_paths(
                topo,
                params,
                vec(path[:, 1]),
                vec(path[:, 2]),
                r;
                sampler="discrete_positive_bridge_gibbs",
            )
        catch err
            last_error = err
        end
    end

    if last_error === nothing
        throw(ErrorException("failed to construct a valid mated_crt map"))
    end
    throw(last_error)
end

generate_mated_crt_map(; kwargs...) = build_mated_crt_map(; kwargs...)

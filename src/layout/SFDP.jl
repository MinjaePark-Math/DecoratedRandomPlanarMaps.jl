# Graphviz SFDP layout wrapper with a built-in force-layout fallback.

function _parse_dot_attrs(attr_text::AbstractString)
    attrs = Dict{String,String}()
    for m in eachmatch(r"(\w+)\s*=\s*(\".*?\"|[^,\]]+)"s, attr_text)
        key = m.captures[1]
        value = strip(m.captures[2])
        if startswith(value, "\"") && endswith(value, "\"") && ncodeunits(value) >= 2
            value = value[2:end-1]
        end
        attrs[key] = value
    end
    return attrs
end

function _parse_pos(raw_pos::AbstractString)
    parts = [rstrip(strip(part), '!') for part in split(raw_pos, ",") if !isempty(strip(part))]
    vals = Float64[parse(Float64, p) for p in parts]
    if length(vals) >= 3
        return vals[1:3]
    elseif length(vals) == 2
        return [vals[1], vals[2], 0.0]
    elseif length(vals) == 1
        return [vals[1], 0.0, 0.0]
    else
        return [0.0, 0.0, 0.0]
    end
end

function _normalize_pointcloud!(pos::Matrix{Float64}, scale::Float64)
    if size(pos, 1) == 0
        return pos
    end
    pos .-= mean(pos; dims=1)
    radius = maximum(sqrt.(sum(pos .^ 2; dims=2)))
    if radius > 0
        pos .*= scale / radius
    end
    return pos
end

function _string_key_dict(options)
    out = Dict{String,Any}()
    for (k, v) in pairs(options)
        out[string(k)] = v
    end
    return out
end

function _compute_force_layout_3d(
    num_vertices::Int,
    edges;
    scale::Float64=1.0,
    seed::Int=7,
    iterations::Int=200,
    K=nothing,
    repulsiveforce::Real=1.0,
)
    n = Int(num_vertices)
    pos = zeros(Float64, n, 3)
    n == 0 && return pos
    n == 1 && return pos

    rng = MersenneTwister(seed)
    pos .= randn(rng, n, 3)
    _normalize_pointcloud!(pos, 1.0)

    edge_array = sanitize_edge_array(edges)
    k = K === nothing ? 1.0 / cbrt(max(n, 1)) : float(K)
    k > 0 || throw(ArgumentError("K must be positive"))
    repulse = float(repulsiveforce)
    repulse >= 0 || throw(ArgumentError("repulsiveforce must be nonnegative"))
    iters = max(1, Int(iterations))

    for it in 1:iters
        disp = zeros(Float64, n, 3)

        for i in 1:(n - 1)
            for j in (i + 1):n
                δ = pos[i, :] .- pos[j, :]
                dist = max(norm(δ), 1.0e-6)
                force = repulse * k^2 / dist
                vecf = (δ ./ dist) .* force
                disp[i, :] .+= vecf
                disp[j, :] .-= vecf
            end
        end

        for r in 1:size(edge_array, 1)
            u = Int(edge_array[r, 1]) + 1
            v = Int(edge_array[r, 2]) + 1
            u == v && continue
            δ = pos[u, :] .- pos[v, :]
            dist = max(norm(δ), 1.0e-6)
            force = dist^2 / k
            vecf = (δ ./ dist) .* force
            disp[u, :] .-= vecf
            disp[v, :] .+= vecf
        end

        t = 0.15 * (1.0 - (it - 1) / iters)
        for i in 1:n
            d = norm(disp[i, :])
            if d > 0
                pos[i, :] .+= disp[i, :] ./ d .* min(d, t)
            end
        end
        pos .*= 0.98
    end

    return _normalize_pointcloud!(pos, scale)
end

function compute_sfdp_layout(
    num_vertices::Integer,
    edges;
    scale::Real=1.0,
    seed::Integer=7,
    K=nothing,
    repulsiveforce=nothing,
    iterations::Union{Nothing,Integer}=nothing,
    overlap=nothing,
    options=Dict{String,Any}(),
    timeout::Union{Nothing,Real}=120.0,
)
    n = Int(num_vertices)
    n >= 0 || throw(ArgumentError("num_vertices must be nonnegative"))
    float(scale) > 0 || throw(ArgumentError("scale must be positive"))

    edge_array = sanitize_edge_array(edges)
    option_dict = _string_key_dict(options)
    K_val = K === nothing ? get(option_dict, "K", nothing) : K
    repulsive_val = repulsiveforce === nothing ? get(option_dict, "repulsiveforce", nothing) : repulsiveforce
    iterations_val = iterations === nothing ? Int(get(option_dict, "iterations", 200)) : Int(iterations)
    overlap_val = overlap === nothing ? get(option_dict, "overlap", "prism") : overlap

    exe = Sys.which("sfdp")
    if exe !== nothing
        graph_options = Dict(
            "K" => "0.3",
            "repulsiveforce" => "1.0",
            "overlap" => string(overlap_val),
            "dim" => "3",
            "dimen" => "3",
        )
        for (k, v) in pairs(option_dict)
            k == "iterations" && continue
            graph_options[string(k)] = string(v)
        end
        K_val !== nothing && (graph_options["K"] = string(K_val))
        repulsive_val !== nothing && (graph_options["repulsiveforce"] = string(repulsive_val))
        overlap_val !== nothing && (graph_options["overlap"] = string(overlap_val))
        attr_str = join(["$k=\"$v\"" for (k, v) in sort!(collect(graph_options); by=first)], ", ")

        dot = IOBuffer()
        println(dot, "graph G { graph [$attr_str, start=\"$(Int(seed))\"]; node [shape=point, width=0.001];")
        for i in 0:(n - 1)
            println(dot, "  n$i;")
        end
        for r in 1:size(edge_array, 1)
            println(dot, "  n$(Int(edge_array[r,1])) -- n$(Int(edge_array[r,2]));")
        end
        println(dot, "}")
        dot_text = String(take!(dot))

        tmp = tempname() * ".dot"
        write(tmp, dot_text)
        try
            output = read(`$exe -Tdot $tmp`, String)
            pos = zeros(Float64, n, 3)
            rx = r"^\s*(n\d+)\s+\[(.*?)\];\s*$"ms
            for m in eachmatch(rx, output)
                index = parse(Int, m.captures[1][2:end]) + 1
                attrs = _parse_dot_attrs(m.captures[2])
                if haskey(attrs, "pos")
                    pos[index, :] .= _parse_pos(attrs["pos"])
                end
            end
            rm(tmp; force=true)
            return _normalize_pointcloud!(pos, float(scale))
        catch
            @warn "Graphviz sfdp failed; falling back to the built-in O(n^2) force layout. Installing/updating Graphviz from graphviz.org is recommended for faster 3D layouts."
            isfile(tmp) && rm(tmp; force=true)
        end
    else
        @warn "Graphviz sfdp was not found on PATH; falling back to the built-in O(n^2) force layout. Installing Graphviz from graphviz.org is recommended for much faster 3D layouts."
    end

    return _compute_force_layout_3d(
        n,
        edge_array;
        scale=float(scale),
        seed=Int(seed),
        iterations=iterations_val,
        K=K_val === nothing ? nothing : float(K_val),
        repulsiveforce=repulsive_val === nothing ? 1.0 : float(repulsive_val),
    )
end

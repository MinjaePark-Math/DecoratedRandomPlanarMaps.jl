# Lightweight SVG preview export for 2D layouts.

function _svg_escape(s::AbstractString)
    return replace(replace(replace(String(s), "&" => "&amp;"), "<" => "&lt;"), ">" => "&gt;")
end

function _finite_row_mask(pos::AbstractMatrix)
    mask = BitVector(undef, size(pos, 1))
    @inbounds for i in axes(pos, 1)
        ok = true
        for j in axes(pos, 2)
            if !isfinite(pos[i, j])
                ok = false
                break
            end
        end
        mask[i] = ok
    end
    return mask
end

function _svg_transform_spec(pos2; width::Int=1200, height::Int=1000, margin::Float64=30.0)
    mask = _finite_row_mask(pos2)
    if !any(mask)
        return 1.0, width / 2, height / 2
    end
    xs = pos2[mask, 1]
    ys = pos2[mask, 2]
    minx, maxx = extrema(xs)
    miny, maxy = extrema(ys)
    spanx = max(maxx - minx, 1e-6)
    spany = max(maxy - miny, 1e-6)
    scale = min((width - 2margin) / spanx, (height - 2margin) / spany)
    tx = (width - scale * (minx + maxx)) / 2
    ty = (height + scale * (miny + maxy)) / 2
    return scale, tx, ty
end

function _svg_apply_transform(pos2, scale::Float64, tx::Float64, ty::Float64)
    pts = Matrix{Float64}(undef, size(pos2, 1), 2)
    pts[:, 1] .= scale .* pos2[:, 1] .+ tx
    pts[:, 2] .= -scale .* pos2[:, 2] .+ ty
    return pts
end

function _svg_transform_points(pos2; width::Int=1200, height::Int=1000, margin::Float64=30.0)
    scale, tx, ty = _svg_transform_spec(pos2; width=width, height=height, margin=margin)
    return _svg_apply_transform(pos2, scale, tx, ty)
end

function _filter_triangles_by_mask(triangles::AbstractMatrix{<:Integer}, validmask::AbstractVector{Bool})
    tri = sanitize_triangles(triangles; drop_degenerate=true, deduplicate=true)
    size(tri, 1) == 0 && return tri
    keep = BitVector(undef, size(tri, 1))
    @inbounds for i in 1:size(tri, 1)
        a = Int(tri[i, 1]) + 1
        b = Int(tri[i, 2]) + 1
        c = Int(tri[i, 3]) + 1
        keep[i] = (1 <= a <= length(validmask)) && (1 <= b <= length(validmask)) && (1 <= c <= length(validmask)) &&
                  validmask[a] && validmask[b] && validmask[c]
    end
    return tri[keep, :]
end

function _filter_edge_array_by_mask(edge_array, validmask::AbstractVector{Bool})
    arr = sanitize_edge_array(edge_array)
    size(arr, 1) == 0 && return arr
    keep = BitVector(undef, size(arr, 1))
    @inbounds for i in 1:size(arr, 1)
        u = Int(arr[i, 1]) + 1
        v = Int(arr[i, 2]) + 1
        keep[i] = (1 <= u <= length(validmask)) && (1 <= v <= length(validmask)) && validmask[u] && validmask[v]
    end
    return arr[keep, :]
end

function _filter_exploration_pairs(exploration::AbstractMatrix)
    arr = Float64.(exploration)
    size(arr, 1) == 0 && return zeros(Float64, 0, size(arr, 2))
    d = size(arr, 2)
    rows = Float64[]
    pair_count = fld(size(arr, 1), 2)
    @inbounds for k in 1:pair_count
        i = 2k - 1
        j = 2k
        ok = true
        for col in 1:d
            if !isfinite(arr[i, col]) || !isfinite(arr[j, col])
                ok = false
                break
            end
        end
        ok || continue
        append!(rows, view(arr, i, :))
        append!(rows, view(arr, j, :))
    end
    isempty(rows) && return zeros(Float64, 0, d)
    return Matrix{Float64}(permutedims(reshape(rows, d, :)))
end

function export_svg_preview(
    map_data,
    pos,
    output_path::AbstractString;
    edge_groups=nothing,
    faces=nothing,
    triangles=nothing,
    show_faces=nothing,
    show_vertices::Bool=false,
    title=nothing,
    metadata=nothing,
)
    path = abspath(String(output_path))
    mkpath(dirname(path))

    pos_arr = Float64.(pos)
    ndims(pos_arr) == 2 || throw(ArgumentError("positions must be a matrix"))
    pos2 = size(pos_arr, 2) == 2 ? pos_arr : pos_arr[:, 1:2]
    validmask = _finite_row_mask(pos2)
    safe_pos2 = copy(pos2)
    for i in axes(safe_pos2, 1)
        validmask[i] && continue
        safe_pos2[i, :] .= 0
    end

    scale, tx, ty = _svg_transform_spec(pos2)
    screen = _svg_apply_transform(safe_pos2, scale, tx, ty)

    tri = surface_triangles(map_data; faces=faces, triangles=triangles)
    tri = _filter_triangles_by_mask(tri, validmask)
    if show_faces === nothing
        show_faces = size(tri, 1) <= 150_000
    end

    groups = grouped_edges(map_data; edge_groups=edge_groups, drop_loops=true)
    exploration = exploration_segment_points(map_data, pos2; edge_groups=edge_groups, faces=faces, triangles=tri, metadata=metadata)
    exploration = _filter_exploration_pairs(exploration)

    open(path, "w") do io
        println(io, """<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 1200 1000\" width=\"1200\" height=\"1000\">""")
        println(io, """<rect x=\"0\" y=\"0\" width=\"1200\" height=\"1000\" fill=\"#ffffff\"/>""")
        if title !== nothing
            println(io, """<text x=\"20\" y=\"30\" font-family=\"system-ui, sans-serif\" font-size=\"20\" fill=\"#111827\">$(_svg_escape(String(title)))</text>""")
        end

        if show_faces && size(tri, 1) > 0
            println(io, """<g fill=\"#fff6bf\" fill-opacity=\"0.24\" stroke=\"none\">""")
            for i in 1:size(tri, 1)
                a, b, c = Int.(tri[i, :]) .+ 1
                println(io, """<polygon points=\"$(screen[a,1]),$(screen[a,2]) $(screen[b,1]),$(screen[b,2]) $(screen[c,1]),$(screen[c,2])\"/>""")
            end
            println(io, "</g>")
        end

        for (name, edge_array) in sort!(collect(groups); by=first)
            filtered = _filter_edge_array_by_mask(edge_array, validmask)
            size(filtered, 1) == 0 && continue
            color = get(EDGE_COLOR_HINTS, name, "#2c3e50")
            width = size(filtered, 1) > 250_000 ? 0.25 : 0.7
            if name in ("red", "blue", "purple", "orange", "green", "navy")
                width *= 1.8
            elseif name == "generic"
                width *= 1.2
            end
            println(io, """<g stroke=\"$color\" stroke-width=\"$width\" stroke-opacity=\"0.9\" fill=\"none\">""")
            for i in 1:size(filtered, 1)
                u = Int(filtered[i, 1]) + 1
                v = Int(filtered[i, 2]) + 1
                println(io, """<line x1=\"$(screen[u,1])\" y1=\"$(screen[u,2])\" x2=\"$(screen[v,1])\" y2=\"$(screen[v,2])\"/>""")
            end
            println(io, "</g>")
        end

        if size(exploration, 1) > 0
            println(io, """<g stroke=\"#27ae60\" stroke-width=\"0.9\" stroke-opacity=\"0.8\" fill=\"none\">""")
            ex2 = _svg_apply_transform(Float64.(exploration[:, 1:2]), scale, tx, ty)
            for i in 1:2:size(ex2, 1)
                println(io, """<line x1=\"$(ex2[i,1])\" y1=\"$(ex2[i,2])\" x2=\"$(ex2[i + 1,1])\" y2=\"$(ex2[i + 1,2])\"/>""")
            end
            println(io, "</g>")
        end

        if show_vertices
            println(io, """<g fill=\"#111827\" fill-opacity=\"0.65\" stroke=\"none\">""")
            radius = size(pos2, 1) > 50_000 ? 1.0 : 2.0
            for i in 1:size(screen, 1)
                validmask[i] || continue
                println(io, """<circle cx=\"$(screen[i,1])\" cy=\"$(screen[i,2])\" r=\"$radius\"/>""")
            end
            println(io, "</g>")
        end

        println(io, "</svg>")
    end

    return path
end

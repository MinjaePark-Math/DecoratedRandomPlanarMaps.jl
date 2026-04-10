module DecoratedRandomPlanarMapsMakieExt

using DecoratedRandomPlanarMaps
using GLMakie
using GeometryBasics

const DRPM = DecoratedRandomPlanarMaps

function _finite_vertex_mask(pos::AbstractMatrix)
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

function _zero_invalid_rows(pos::AbstractMatrix, validmask::AbstractVector{Bool})
    out = Float32.(pos)
    @inbounds for i in axes(out, 1)
        validmask[i] && continue
        out[i, :] .= 0
    end
    return out
end

function _as_points_2d(pos::AbstractMatrix)
    return [Point2f(Float32(pos[i, 1]), Float32(pos[i, 2])) for i in axes(pos, 1)]
end

function _as_points_3d(pos::AbstractMatrix)
    pos3 = DRPM.pad_positions_3d(pos)
    return [Point3f(pos3[i, 1], pos3[i, 2], pos3[i, 3]) for i in axes(pos3, 1)]
end

function _triangle_faces(triangles)
    tri = DRPM.sanitize_triangles(triangles; drop_degenerate=true, deduplicate=true)
    out = TriangleFace{Int}[]
    size(tri, 1) == 0 && return out
    sizehint!(out, size(tri, 1))
    for i in 1:size(tri, 1)
        push!(out, TriangleFace{Int}(Int(tri[i, 1]) + 1, Int(tri[i, 2]) + 1, Int(tri[i, 3]) + 1))
    end
    return out
end

function _filter_triangles_for_valid_rows(triangles, validmask::AbstractVector{Bool})
    tri = DRPM.sanitize_triangles(triangles; drop_degenerate=true, deduplicate=true)
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

function _line_points(points, edges, validmask::AbstractVector{Bool})
    arr = DRPM.sanitize_edge_array(edges)
    if isempty(points)
        return eltype(points)[]
    end
    segs = typeof(points[1])[]
    size(arr, 1) == 0 && return segs
    sizehint!(segs, 2 * size(arr, 1))
    @inbounds for i in 1:size(arr, 1)
        u = Int(arr[i, 1]) + 1
        v = Int(arr[i, 2]) + 1
        if u == v || !(1 <= u <= length(validmask)) || !(1 <= v <= length(validmask))
            continue
        end
        (validmask[u] && validmask[v]) || continue
        push!(segs, points[u], points[v])
    end
    return segs
end

function _filter_exploration_pairs(exploration::AbstractMatrix)
    arr = Float32.(exploration)
    size(arr, 1) == 0 && return zeros(Float32, 0, size(arr, 2))
    size(arr, 2) in (2, 3) || throw(ArgumentError("exploration positions must have shape (N, 2) or (N, 3)"))

    rows = Float32[]
    d = size(arr, 2)
    sizehint!(rows, length(arr))

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

    isempty(rows) && return zeros(Float32, 0, d)
    return Matrix{Float32}(permutedims(reshape(rows, d, :)))
end

function _exploration_points(exploration::AbstractMatrix, dim::Int)
    return dim == 3 ? _as_points_3d(exploration) : _as_points_2d(exploration)
end

function _family(map_data, groups; metadata=nothing)
    if map_data isa DRPM.UniformMap
        return :uniform
    elseif map_data isa DRPM.SchnyderMap
        return :schnyder
    elseif DRPM.should_include_exploration(map_data; edge_groups=groups, metadata=metadata)
        return :fk
    else
        return :generic
    end
end

function _geometry(map_data, pos; edge_groups=nothing, faces=nothing, triangles=nothing, metadata=nothing, dim::Int)
    pos_arr = Float32.(pos)
    validmask = _finite_vertex_mask(pos_arr)
    safe_pos = _zero_invalid_rows(pos_arr, validmask)

    tri = DRPM.surface_triangles(map_data; faces=faces, triangles=triangles, drop_degenerate=true, deduplicate=true)
    tri = _filter_triangles_for_valid_rows(tri, validmask)

    groups = DRPM.grouped_edges(map_data; edge_groups=edge_groups, drop_loops=true)
    exploration = DRPM.exploration_segment_points(map_data, pos_arr; edge_groups=edge_groups, faces=faces, triangles=tri, metadata=metadata)
    exploration = _filter_exploration_pairs(exploration)

    points = dim == 3 ? _as_points_3d(safe_pos) : _as_points_2d(safe_pos)
    expl_points = size(exploration, 1) == 0 ? (dim == 3 ? Point3f[] : Point2f[]) : _exploration_points(exploration, dim)
    return points, tri, groups, expl_points, validmask
end

function _set_visibility!(items, is_visible::Bool)
    for item in items
        if hasproperty(item, :blockscene)
            item.blockscene.visible[] = is_visible
        elseif hasproperty(item, :visible)
            item.visible[] = is_visible
        end
    end
    return nothing
end

function _default_toggle_on(label::AbstractString, family::Symbol, dim::Int)
    label == "exploration" && return false
    family == :schnyder && dim == 2 && label == "generic" && return false
    return true
end

function _control_panel!(fig, entries::Vector{Pair{String,String}}; title::AbstractString="Controls", default_alpha::Float64=0.55, family::Symbol=:generic, dim::Int=2)
    ui_width = 260
    panel = fig[1, 2] = GridLayout(tellheight=false, valign=:top, width=ui_width)
    row = 1

    panel_items = Any[]

    title_label = Label(panel[row, 1:2], title; fontsize=18, halign=:left, font=:bold)
    push!(panel_items, title_label)
    row += 1

    mesh_toggle = Toggle(panel[row, 1], active=true)
    mesh_label = Label(panel[row, 2], "Show faces"; halign=:left)
    append!(panel_items, (mesh_toggle, mesh_label))
    row += 1

    group_toggles = Dict{String,Toggle}()
    for (key, label) in entries
        t = Toggle(panel[row, 1], active=_default_toggle_on(label, family, dim))
        lbl = Label(panel[row, 2], label; halign=:left)
        group_toggles[key] = t
        append!(panel_items, (t, lbl))
        row += 1
    end

    alpha_label = Label(panel[row, 1:2], "Mesh opacity"; halign=:left)
    push!(panel_items, alpha_label)
    row += 1

    alpha = Slider(panel[row, 1:2], range=0.0:0.01:1.0, startvalue=default_alpha)
    push!(panel_items, alpha)
    row += 1

    reset = Button(panel[row, 1:2], label="Reset view", halign=:center)
    push!(panel_items, reset)

    hide_button = Button(fig[1, 1], label="Hide panel", halign=:right, valign=:top, tellwidth=false, tellheight=false)

    panel_hidden = Observable(false)
    on(hide_button.clicks) do _
        panel_hidden[] = !panel_hidden[]
    end
    on(panel_hidden) do hidden
        _set_visibility!(panel_items, !hidden)
        colsize!(fig.layout, 2, hidden ? 0 : ui_width)
        hide_button.label[] = hidden ? "Show panel" : "Hide panel"
    end

    return mesh_toggle, group_toggles, alpha, reset, hide_button
end

function _mesh_color(alpha_obs)
    return lift(alpha_obs) do a
        RGBAf(0.92f0, 0.89f0, 0.82f0, Float32(a))
    end
end

function _line_style(family::Symbol, name::AbstractString, dim::Int)
    color = get(DRPM.EDGE_COLOR_HINTS, name, "#2c3e50")
    linewidth = dim == 2 ? 1.6 : 2.0
    if family == :uniform
        if name == "generic"
            color = "#7f8c8d"
            linewidth = dim == 2 ? 1.4 : 1.6
        end
    elseif family == :schnyder
        if name == "generic"
            color = "#b0b7bf"
            linewidth = dim == 2 ? 0.8 : 1.2
        elseif name == "green"
            color = "#1f7a3a"
            linewidth = dim == 2 ? 2.6 : 2.4
        elseif name in ("orange", "navy")
            linewidth = dim == 2 ? 2.6 : 2.4
        end
    elseif family == :fk
        if name == "generic"
            color = "#b0b7bf"
            linewidth = dim == 2 ? 1.2 : 1.4
        end
    end
    return (color=color, linewidth=linewidth)
end

function _exploration_style(dim::Int)
    return (color=get(DRPM.EDGE_COLOR_HINTS, "exploration", "#27ae60"), linewidth=(dim == 2 ? 1.2 : 2.4))
end


function _resolved_title(title::AbstractString, map_data, metadata, dim::Int)
    if strip(String(title)) in ("2D map", "3D map")
        return DRPM.default_render_title(map_data; metadata=metadata, dimension=dim, engine=(dim == 2 ? "Tutte" : "SFDP"))
    end
    params = DRPM.model_parameter_string(map_data; metadata=metadata)
    isempty(params) && return String(title)
    occursin(params, String(title)) && return String(title)
    return "$(title) · $(params)"
end

function DecoratedRandomPlanarMaps.render_makie_3d(map_data, pos; edge_groups=nothing, faces=nothing, triangles=nothing, metadata=nothing, title::AbstractString="3D map")
    pos3 = DRPM.pad_positions_3d(pos)
    points, tri, groups, expl_points, validmask = _geometry(map_data, pos3; edge_groups=edge_groups, faces=faces, triangles=triangles, metadata=metadata, dim=3)
    family = _family(map_data, groups; metadata=metadata)

    full_title = _resolved_title(title, map_data, metadata, 3)

    fig = Figure(size=(1280, 840), backgroundcolor=:white)
    ax = LScene(fig[1, 1], show_axis=false)

    group_keys = sort!(collect(keys(groups)))
    entries = Pair{String,String}[]
    for key in group_keys
        push!(entries, key => key)
    end
    if DRPM.should_include_exploration(map_data; edge_groups=groups, metadata=metadata)
        push!(entries, "__exploration__" => "exploration")
    end

    mesh_toggle, group_toggles, alpha, reset, _ = _control_panel!(fig, entries; title=full_title, default_alpha=0.55, family=family, dim=3)

    if !isempty(tri)
        mesh_obj = GeometryBasics.Mesh(points, _triangle_faces(tri))
        mesh!(ax, mesh_obj; color=_mesh_color(alpha.value), visible=mesh_toggle.active, shading=true, backlight=1.0f0, transparency=true)
    else
        mesh_toggle.active[] = false
    end

    for name in group_keys
        segs = _line_points(points, groups[name], validmask)
        isempty(segs) && continue
        style = _line_style(family, name, 3)
        linesegments!(ax, segs; color=style.color, linewidth=style.linewidth, visible=group_toggles[name].active)
    end

    if haskey(group_toggles, "__exploration__") && !isempty(expl_points)
        estyle = _exploration_style(3)
        linesegments!(ax, expl_points; color=estyle.color, linewidth=estyle.linewidth, visible=group_toggles["__exploration__"].active)
        group_toggles["__exploration__"].active[] = false
    elseif haskey(group_toggles, "__exploration__")
        group_toggles["__exploration__"].active[] = false
    end

    on(reset.clicks) do _
        center!(ax.scene)
    end

    display(fig)
    return fig
end

function DecoratedRandomPlanarMaps.render_makie_2d(map_data, pos; edge_groups=nothing, faces=nothing, triangles=nothing, metadata=nothing, title::AbstractString="2D map")
    size(pos, 2) == 2 || throw(ArgumentError("render_makie_2d expects a position matrix with two columns"))
    points, tri, groups, expl_points, validmask = _geometry(map_data, pos; edge_groups=edge_groups, faces=faces, triangles=triangles, metadata=metadata, dim=2)
    family = _family(map_data, groups; metadata=metadata)

    full_title = _resolved_title(title, map_data, metadata, 2)

    fig = Figure(size=(1280, 840), backgroundcolor=:white)
    ax = Axis(fig[1, 1], aspect=DataAspect(), title=full_title)
    hidedecorations!(ax)
    hidespines!(ax)

    group_keys = sort!(collect(keys(groups)))
    entries = Pair{String,String}[]
    for key in group_keys
        push!(entries, key => key)
    end
    if DRPM.should_include_exploration(map_data; edge_groups=groups, metadata=metadata)
        push!(entries, "__exploration__" => "exploration")
    end

    mesh_toggle, group_toggles, alpha, reset, _ = _control_panel!(fig, entries; title=full_title, default_alpha=0.25, family=family, dim=2)

    if !isempty(tri)
        mesh_obj = GeometryBasics.Mesh(points, _triangle_faces(tri))
        mesh!(ax, mesh_obj; color=_mesh_color(alpha.value), visible=mesh_toggle.active, shading=NoShading, transparency=true)
    else
        mesh_toggle.active[] = false
    end

    for name in group_keys
        segs = _line_points(points, groups[name], validmask)
        isempty(segs) && continue
        style = _line_style(family, name, 2)
        line_plot = linesegments!(ax, segs; color=style.color, linewidth=style.linewidth, visible=group_toggles[name].active)
        translate!(line_plot, 0, 0, 1)
        if family == :schnyder && name == "generic"
            group_toggles[name].active[] = false
        end
    end

    if haskey(group_toggles, "__exploration__") && !isempty(expl_points)
        estyle = _exploration_style(2)
        expl_plot = linesegments!(ax, expl_points; color=estyle.color, linewidth=estyle.linewidth, visible=group_toggles["__exploration__"].active)
        translate!(expl_plot, 0, 0, 2)
        group_toggles["__exploration__"].active[] = false
    elseif haskey(group_toggles, "__exploration__")
        group_toggles["__exploration__"].active[] = false
    end

    on(reset.clicks) do _
        autolimits!(ax)
    end

    display(fig)
    return fig
end

end

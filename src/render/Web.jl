# Three.js-style binary export helpers.

function _write_binary(path::AbstractString, array::AbstractVector)
    open(path, "w") do io
        write(io, array)
    end
end

# Julia matrices are column-major, but WebGL/Three.js vertex/index buffers are
# expected in interleaved row-major order.
function _write_binary(path::AbstractString, array::AbstractMatrix)
    open(path, "w") do io
        write(io, vec(copy(permutedims(array))))
    end
end

_write_vertex_positions(path::AbstractString, pos::AbstractMatrix{<:AbstractFloat}) =
    _write_binary(path, pos)
_write_line_positions(path::AbstractString, pos::AbstractMatrix{<:AbstractFloat}) =
    _write_binary(path, pos)
_write_edge_indices(path::AbstractString, edges::AbstractMatrix{<:Integer}) =
    _write_binary(path, UInt32.(edges))
_write_face_indices(path::AbstractString, faces::AbstractMatrix{<:Integer}) =
    _write_binary(path, UInt32.(faces))

function _hex_to_rgb01(hex::AbstractString)
    r, g, b = _hex_to_rgb(hex)
    return (Float32(r / 255), Float32(g / 255), Float32(b / 255))
end

function _web_visual_hex(name::AbstractString; background::AbstractString="#F6F8FB")
    style = _edge_group_style(name)
    if _is_meandric_fade_group(name)
        return _hex_blend_over(style.color, background, style.opacity)
    end
    return style.color
end

function _is_meandric_fade_group(name)
    return startswith(_normalize_group_string(name), MEANDRIC_LOOP_FADE_PREFIX)
end

function _append_color_rows!(dest::Vector{Float32}, color::NTuple{3,Float32}, count::Integer)
    n = max(0, Int(count))
    sizehint!(dest, length(dest) + 3n)
    for _ in 1:n
        push!(dest, color[1], color[2], color[3])
    end
    return dest
end

function _merge_meandric_fade_groups(groups)
    merged_edges = Matrix{Int32}(undef, 0, 2)
    merged_colors = Matrix{Float32}(undef, 0, 3)

    fade_entries = Tuple{String,Matrix{Int32}}[]
    for (name, edge_array) in groups
        _is_meandric_fade_group(name) || continue
        arr = sanitize_edge_array(edge_array)
        size(arr, 1) == 0 && continue
        push!(fade_entries, (String(name), Int32.(arr)))
    end

    isempty(fade_entries) && return merged_edges, merged_colors

    total_edges = sum(size(arr, 1) for (_, arr) in fade_entries)
    merged_edges = Matrix{Int32}(undef, total_edges, 2)
    color_rows = Float32[]
    sizehint!(color_rows, 3 * total_edges)

    row = 1
    for (name, arr) in sort!(fade_entries; by=first)
        count = size(arr, 1)
        merged_edges[row:(row + count - 1), :] .= arr
        _append_color_rows!(color_rows, _hex_to_rgb01(_web_visual_hex(name)), count)
        row += count
    end

    merged_colors = reshape(color_rows, 3, :)'
    return merged_edges, merged_colors
end

function _slugify(name::AbstractString)
    slug = replace(lowercase(strip(string(name))), r"[^0-9a-zA-Z_\-]+" => "_")
    slug = strip(slug, '_')
    return isempty(slug) ? "group" : slug
end

function _web_model_name(map_data, metadata)
    return model_display_name(map_data; metadata=metadata)
end

function _web_parameter_dict(map_data, metadata)
    return Dict(k => v for (k, v) in model_parameter_pairs(map_data; metadata=metadata))
end

function _display_edge_group_name(name::AbstractString)
    s = string(name)
    m = match(r"^zz_loop_(\d+)$", s)
    m !== nothing && return "loop " * string(parse(Int, m.captures[1]))
    m = match(r"^zy_loop_fade_(\d+)$", s)
    m !== nothing && return "loop fade " * string(parse(Int, m.captures[1]))
    return s
end

function _edge_group_legend_name(name::AbstractString)
    s = string(name)
    startswith(s, "zy_loop_fade_") && return "loop fade"
    return _display_edge_group_name(s)
end

function _edge_group_default_visible(name::AbstractString)
    name in ("exploration", "tree_upper", "tree_lower", "glue_upper", "glue_lower") && return false
    return true
end

function _merge_metadata_preserving_generated!(meta_dict::Dict{String,Any}, metadata)
    metadata === nothing && return meta_dict
    for (k, v) in pairs(metadata)
        key = string(k)
        if key == "edge_groups"
            continue
        elseif key == "sphere_circle_geometry"
            continue
        elseif key == "exploration_green_edge_pairs"
            continue
        elseif key == "files"
            if v !== nothing
                files_dict = meta_dict["files"]
                for (fk, fv) in pairs(v)
                    files_dict[string(fk)] = fv
                end
            end
        else
            meta_dict[key] = v
        end
    end
    return meta_dict
end

function export_web_binaries(
    map_data,
    pos,
    output_dir::AbstractString;
    edge_groups=nothing,
    faces=nothing,
    triangles=nothing,
    circle_radii=nothing,
    sphere_circle_geometry=nothing,
    metadata=nothing,
    write_index_html::Bool=true,
)
    out_path = abspath(string(output_dir))
    mkpath(out_path)

    pos_arr = Float32.(pos)
    ndims(pos_arr) == 2 || throw(ArgumentError("positions must be a matrix"))
    input_dim = size(pos_arr, 2)
    input_dim in (2, 3) || throw(ArgumentError("positions must have shape (N, 2) or (N, 3)"))
    _write_vertex_positions(joinpath(out_path, "data_verts.bin"), pos_arr)

    radii_file = nothing
    if circle_radii !== nothing
        radii = Float32.(collect(circle_radii))
        length(radii) == size(pos_arr, 1) || throw(ArgumentError("circle_radii must have length equal to the number of vertices"))
        radii_file = "data_radii.bin"
        _write_binary(joinpath(out_path, radii_file), radii)
    end

    tri_array = surface_triangles(map_data; faces=faces, triangles=triangles, drop_degenerate=true, deduplicate=true)
    face_file = nothing
    if size(tri_array, 1) > 0
        face_file = "data_faces.bin"
        _write_face_indices(joinpath(out_path, face_file), tri_array)
    end

    groups = grouped_edges(map_data; edge_groups=edge_groups, drop_loops=true)
    edge_meta = Vector{Any}()
    total_edges = 0

    fade_edges, fade_colors = _merge_meandric_fade_groups(groups)

    for (name, edge_array) in sort!(collect(groups); by=first)
        _is_meandric_fade_group(name) && continue
        size(edge_array, 1) == 0 && continue
        display_name = _display_edge_group_name(name)
        slug = _slugify(display_name)
        filename = "data_edges_$(slug).bin"
        _write_edge_indices(joinpath(out_path, filename), edge_array)
        push!(edge_meta, Dict(
            "name" => display_name,
            "legend_name" => _edge_group_legend_name(String(name)),
            "file" => filename,
            "color" => edge_group_color(String(name)),
            "opacity" => edge_group_opacity(String(name)),
            "linewidth" => 1.15 * edge_group_width_scale(String(name)),
            "wide" => edge_group_force_wide(String(name)),
            "num_edges" => size(edge_array, 1),
            "kind" => "index_pairs",
            "default_visible" => _edge_group_default_visible(display_name),
        ))
        total_edges += size(edge_array, 1)
    end

    if size(fade_edges, 1) > 0
        fade_name = "loop fade"
        fade_slug = _slugify(fade_name)
        fade_file = "data_edges_$(fade_slug).bin"
        fade_color_file = "data_edge_colors_$(fade_slug).bin"
        _write_edge_indices(joinpath(out_path, fade_file), fade_edges)
        _write_binary(joinpath(out_path, fade_color_file), fade_colors)
        push!(edge_meta, Dict(
            "name" => fade_name,
            "legend_name" => fade_name,
            "file" => fade_file,
            "color_file" => fade_color_file,
            "color" => _web_visual_hex(_meandric_fade_group_name(length(MEANDRIC_HIGHLIGHT_GROUPS) + 1)),
            "opacity" => 1.0,
            "linewidth" => 1.15 * edge_group_width_scale(_meandric_fade_group_name(length(MEANDRIC_HIGHLIGHT_GROUPS) + 1)),
            "wide" => false,
            "num_edges" => size(fade_edges, 1),
            "kind" => "colored_index_pairs",
            "default_visible" => true,
        ))
        total_edges += size(fade_edges, 1)
    end

    include_exploration = should_include_exploration(map_data; edge_groups=edge_groups, metadata=metadata)
    exploration_file = nothing
    if include_exploration
        exploration = exploration_segment_points(
            map_data,
            pos_arr;
            edge_groups=edge_groups,
            faces=faces,
            triangles=tri_array,
            metadata=metadata,
        )
        exploration_file = "data_exploration.bin"
        _write_line_positions(joinpath(out_path, exploration_file), exploration)
        push!(edge_meta, Dict(
            "name" => "exploration",
            "file" => exploration_file,
            "color" => edge_group_color("exploration"),
            "opacity" => edge_group_opacity("exploration"),
            "linewidth" => 1.15 * edge_group_width_scale("exploration"),
            "wide" => edge_group_force_wide("exploration"),
            "num_edges" => Int(size(exploration, 1) ÷ 2),
            "kind" => "positions",
            "default_visible" => false,
        ))
    end

    sphere_circle_file = nothing
    sphere_circle_caps_file = nothing
    if sphere_circle_geometry !== nothing
        circle_segments = sphere_circle_segment_points(sphere_circle_geometry)
        if size(circle_segments, 1) > 0
            sphere_circle_file = "data_sphere_circles.bin"
            _write_line_positions(joinpath(out_path, sphere_circle_file), circle_segments)
            push!(edge_meta, Dict(
                "name" => "circles",
                "file" => sphere_circle_file,
                "color" => edge_group_color("circles"),
                "opacity" => edge_group_opacity("circles"),
                "linewidth" => 1.15 * edge_group_width_scale("circles"),
                "wide" => edge_group_force_wide("circles"),
                "num_edges" => Int(size(circle_segments, 1) ÷ 2),
                "kind" => "positions",
                "default_visible" => true,
            ))
        end

        centers_f = Float64.(sphere_circle_geometry["centers"])
        normals_f = Float64.(sphere_circle_geometry["normals"])
        radii_f = Float64.(collect(sphere_circle_geometry["radii"]))
        offsets_f = Float64.(collect(sphere_circle_geometry["offsets"]))

        caps = hcat(
            Float32.(centers_f),
            Float32.(normals_f),
            Float32.(reshape(radii_f, :, 1)),
            Float32.(reshape(offsets_f, :, 1)),
        )
        if size(caps, 1) > 0
            sphere_circle_caps_file = "data_sphere_circle_caps.bin"
            _write_binary(joinpath(out_path, sphere_circle_caps_file), caps)
        end
    end

    files_dict = Dict{String,Any}(
        "vertices" => "data_verts.bin",
        "faces" => face_file,
        "exploration" => exploration_file,
        "circle_radii" => radii_file,
        "sphere_circles" => sphere_circle_file,
        "sphere_circle_caps" => sphere_circle_caps_file,
    )

    param_pairs = _web_parameter_dict(map_data, metadata)
    param_summary = model_parameter_string(map_data; metadata=metadata)

    meta_dict = Dict{String,Any}(
        "model" => _web_model_name(map_data, metadata),
        "parameters" => param_pairs,
        "parameters_summary" => param_summary,
        "dimension" => input_dim,
        "num_vertices" => size(pos_arr, 1),
        "num_edges" => total_edges,
        "edge_groups" => edge_meta,
        "files" => files_dict,
    )
    _merge_metadata_preserving_generated!(meta_dict, metadata)
    if sphere_circle_geometry !== nothing
        meta_dict["sphere_radius"] = sphere_circle_geometry["sphere_radius"]
        meta_dict["sphere_effective_projection_scale"] = get(sphere_circle_geometry, "projection_scale", get(meta_dict, "sphere_effective_projection_scale", get(meta_dict, "sphere_projection_scale", 1.0)))
        if !haskey(meta_dict, "sphere_projection_scale")
            meta_dict["sphere_projection_scale"] = meta_dict["sphere_effective_projection_scale"]
        end
    end

    open(joinpath(out_path, "web_meta.json"), "w") do io
        print(io, JSON3.write(meta_dict; allow_inf=true))
    end

    if write_index_html
        src = joinpath(@__DIR__, "..", "..", "assets", "index.html")
        cp(src, joinpath(out_path, "index.html"); force=true)
    end

    return out_path
end

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

_display_edge_group_name(name::AbstractString) = string(name)

function _edge_group_default_visible(name::AbstractString)
    return name == "exploration" ? false : true
end

function _merge_metadata_preserving_generated!(meta_dict::Dict{String,Any}, metadata)
    metadata === nothing && return meta_dict
    for (k, v) in pairs(metadata)
        key = string(k)
        if key == "edge_groups"
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

    tri_array = surface_triangles(map_data; faces=faces, triangles=triangles, drop_degenerate=true, deduplicate=true)
    face_file = nothing
    if size(tri_array, 1) > 0
        face_file = "data_faces.bin"
        _write_face_indices(joinpath(out_path, face_file), tri_array)
    end

    groups = grouped_edges(map_data; edge_groups=edge_groups, drop_loops=true)
    edge_meta = Vector{Any}()
    total_edges = 0
    for (name, edge_array) in sort!(collect(groups); by=first)
        size(edge_array, 1) == 0 && continue
        display_name = _display_edge_group_name(name)
        slug = _slugify(display_name)
        filename = "data_edges_$(slug).bin"
        _write_edge_indices(joinpath(out_path, filename), edge_array)
        push!(edge_meta, Dict(
            "name" => display_name,
            "file" => filename,
            "color" => get(EDGE_COLOR_HINTS, display_name, "#2c3e50"),
            "num_edges" => size(edge_array, 1),
            "kind" => "index_pairs",
            "default_visible" => _edge_group_default_visible(display_name),
        ))
        total_edges += size(edge_array, 1)
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
            "color" => get(EDGE_COLOR_HINTS, "exploration", "#27ae60"),
            "num_edges" => Int(size(exploration, 1) ÷ 2),
            "kind" => "positions",
            "default_visible" => false,
        ))
    end

    files_dict = Dict{String,Any}(
        "vertices" => "data_verts.bin",
        "faces" => face_file,
        "exploration" => exploration_file,
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

    open(joinpath(out_path, "web_meta.json"), "w") do io
        print(io, JSON3.write(meta_dict; allow_inf=true))
    end

    if write_index_html
        src = joinpath(@__DIR__, "..", "..", "assets", "index.html")
        cp(src, joinpath(out_path, "index.html"); force=true)
    end

    return out_path
end

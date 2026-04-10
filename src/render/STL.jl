# Binary STL export.

@inline function _triangle_area_cross(a, b, c)
    ab = b .- a
    ac = c .- a
    return cross(Float32.(ab), Float32.(ac))
end

function export_stl_binary(map_data, pos, output_path::AbstractString; faces=nothing, triangles=nothing, area_eps::Float32=1.0f-8)
    path = abspath(String(output_path))
    mkpath(dirname(path))

    tri = surface_triangles(map_data; faces=faces, triangles=triangles, drop_degenerate=true, deduplicate=true)
    pos3 = Float32.(pad_positions_3d(pos))

    kept_pts = Vector{NTuple{3,NTuple{3,Float32}}}()
    kept_normals = Vector{NTuple{3,Float32}}()
    sizehint!(kept_pts, size(tri, 1))
    sizehint!(kept_normals, size(tri, 1))

    for i in 1:size(tri, 1)
        ia = Int(tri[i, 1]) + 1
        ib = Int(tri[i, 2]) + 1
        ic = Int(tri[i, 3]) + 1
        a = pos3[ia, :]
        b = pos3[ib, :]
        c = pos3[ic, :]
        all(isfinite, a) && all(isfinite, b) && all(isfinite, c) || continue
        n = _triangle_area_cross(a, b, c)
        len = sqrt(sum(n .^ 2))
        len > area_eps || continue
        push!(kept_normals, (n[1] / len, n[2] / len, n[3] / len))
        push!(kept_pts, ((a[1], a[2], a[3]), (b[1], b[2], b[3]), (c[1], c[2], c[3])))
    end

    open(path, "w") do io
        write(io, zeros(UInt8, 80))
        write(io, UInt32(length(kept_pts)))
        for i in eachindex(kept_pts)
            n = kept_normals[i]
            tri_pts = kept_pts[i]
            write(io, Float32[n[1], n[2], n[3]])
            write(io, Float32[tri_pts[1][1], tri_pts[1][2], tri_pts[1][3]])
            write(io, Float32[tri_pts[2][1], tri_pts[2][2], tri_pts[2][3]])
            write(io, Float32[tri_pts[3][1], tri_pts[3][2], tri_pts[3][3]])
            write(io, UInt16(0))
        end
    end

    return path
end

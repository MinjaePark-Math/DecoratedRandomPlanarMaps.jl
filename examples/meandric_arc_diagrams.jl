using DecoratedRandomPlanarMaps

function adjacency_pairs(adj)
    pairs = Tuple{Int,Int}[]
    unmatched = Int[]
    for i in eachindex(adj)
        j = Int(adj[i])
        if j > i - 1
            push!(pairs, (i - 1, j))
        elseif j < 0
            push!(unmatched, i - 1)
        end
    end
    return pairs, unmatched
end

function arc_path(x1, x2, y0; upper::Bool, curvature::Float64=0.45)
    span = abs(x2 - x1)
    height = max(18.0, curvature * span)
    ctrl_y = upper ? (y0 - height) : (y0 + height)
    mid_x = 0.5 * (x1 + x2)
    return "M $(round(x1; digits=2)) $(round(y0; digits=2)) Q $(round(mid_x; digits=2)) $(round(ctrl_y; digits=2)) $(round(x2; digits=2)) $(round(y0; digits=2))"
end

function ray_path(x, y0, length)
    return "M $(round(x; digits=2)) $(round(y0; digits=2)) L $(round(x; digits=2)) $(round(y0 + length; digits=2))"
end

function write_arc_diagram_svg(path::AbstractString, upper_adj, lower_adj;
    title::AbstractString,
    lower_unmatched_label::AbstractString="boundary",
)
    upper_pairs, _ = adjacency_pairs(upper_adj)
    lower_pairs, lower_unmatched = adjacency_pairs(lower_adj)
    npoints = length(upper_adj)
    npoints == length(lower_adj) || throw(ArgumentError("upper and lower adjacencies must have the same length"))

    spacing = 24.0
    margin_x = 36.0
    margin_y = 36.0
    ray_length = 42.0
    label_gap = 14.0
    function xcoord(v0)
        return margin_x + spacing * v0
    end

    upper_span = isempty(upper_pairs) ? 0 : maximum(j - i for (i, j) in upper_pairs)
    lower_span = isempty(lower_pairs) ? 0 : maximum(j - i for (i, j) in lower_pairs)
    upper_height = max(24.0, 0.45 * spacing * upper_span)
    lower_height = max(24.0, 0.45 * spacing * lower_span)

    baseline_y = margin_y + upper_height + 16.0
    width = 2margin_x + spacing * (npoints - 1)
    height = baseline_y + lower_height + (isempty(lower_unmatched) ? 0.0 : ray_length) + margin_y + 18.0

    io = IOBuffer()
    println(io, """<svg xmlns="http://www.w3.org/2000/svg" width="$(round(Int, ceil(width)))" height="$(round(Int, ceil(height)))" viewBox="0 0 $(round(Int, ceil(width))) $(round(Int, ceil(height)))">""")
    println(io, """  <rect x="0" y="0" width="100%" height="100%" fill="white"/>""")
    println(io, """  <text x="$(margin_x)" y="22" font-family="monospace" font-size="16" fill="#111111">$(title)</text>""")
    println(io, """  <text x="$(margin_x)" y="$(height - 10)" font-family="monospace" font-size="11" fill="#666666">upper = blue, lower = red, unmatched lower openings = $(lower_unmatched_label)</text>""")

    for i in 0:(npoints - 2)
        println(
            io,
            """  <line x1="$(round(xcoord(i); digits=2))" y1="$(round(baseline_y; digits=2))" x2="$(round(xcoord(i + 1); digits=2))" y2="$(round(baseline_y; digits=2))" stroke="#b0b7bf" stroke-width="1.3"/>""",
        )
    end

    for (u, v) in upper_pairs
        println(
            io,
            """  <path d="$(arc_path(xcoord(u), xcoord(v), baseline_y; upper=true))" fill="none" stroke="#2980b9" stroke-width="2.0" stroke-linecap="round"/>""",
        )
    end

    for (u, v) in lower_pairs
        println(
            io,
            """  <path d="$(arc_path(xcoord(u), xcoord(v), baseline_y; upper=false))" fill="none" stroke="#c0392b" stroke-width="2.0" stroke-linecap="round"/>""",
        )
    end

    for v in lower_unmatched
        println(
            io,
            """  <path d="$(ray_path(xcoord(v), baseline_y, ray_length))" fill="none" stroke="#c0392b" stroke-width="2.0" stroke-linecap="round"/>""",
        )
    end

    for i in 0:(npoints - 1)
        x = xcoord(i)
        println(
            io,
            """  <text x="$(round(x; digits=2))" y="$(round(baseline_y + label_gap; digits=2))" text-anchor="middle" font-family="monospace" font-size="9" fill="#444444">$(i + 1)</text>""",
        )
    end

    println(io, "</svg>")

    mkpath(dirname(path))
    write(path, String(take!(io)))
    return abspath(path)
end

order = 24
boundary_half_length = 6
seed = 17
outroot = joinpath(@__DIR__, "out", "meandric_arc_diagrams")
mkpath(outroot)

uniform_system = build_uniform_meandric_system(; order=order, seed=seed)
half_plane_system = build_half_plane_meandric_system(; order=order, boundary_half_length=boundary_half_length, seed=seed)

uniform_path = write_arc_diagram_svg(
    joinpath(outroot, "uniform_meandric_arc_diagram.svg"),
    uniform_system.upper_adj,
    uniform_system.lower_adj;
    title="uniform_meandric · order=$(order) · seed=$(seed)",
)

half_plane_path = write_arc_diagram_svg(
    joinpath(outroot, "half_plane_meandric_arc_diagram.svg"),
    half_plane_system.upper_adj,
    half_plane_system.lower_adj;
    title="half_plane_meandric · order=$(order) · boundary_half_length=$(boundary_half_length) · seed=$(seed)",
    lower_unmatched_label="downward boundary rays",
)

println("wrote ", uniform_path)
println("wrote ", half_plane_path)

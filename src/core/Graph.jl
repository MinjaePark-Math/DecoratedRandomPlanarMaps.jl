# Sparse graph utilities shared by layout and rendering code.

struct WeightedGraph
    num_vertices::Int
    edges::Matrix{Int32}            # 0-based vertex ids
    weights::Vector{Float64}
    adjacency::SparseMatrixCSC{Float64,Int}
    laplacian::SparseMatrixCSC{Float64,Int}
    degrees::Vector{Float64}
end

function _empty_edges()
    return Matrix{Int32}(undef, 0, 2)
end

function sanitize_edge_array(edges)::Matrix{Int32}
    if edges isa AbstractMatrix
        arr = Int64.(edges)
        if isempty(arr)
            return _empty_edges()
        end
        size(arr, 2) == 2 || throw(ArgumentError("edge array must have shape (E, 2)"))
        return Int32.(arr)
    elseif edges isa AbstractVector
        isempty(edges) && return _empty_edges()
        arr = Matrix{Int32}(undef, length(edges), 2)
        for (i, e) in enumerate(edges)
            arr[i, 1] = Int32(Int(e[1]))
            arr[i, 2] = Int32(Int(e[2]))
        end
        return arr
    else
        throw(ArgumentError("unsupported edge container"))
    end
end

function collapse_undirected_edges(
    edges;
    weights=nothing,
    drop_loops::Bool=true,
)
    edge_array = sanitize_edge_array(edges)
    if size(edge_array, 1) == 0
        return _empty_edges(), Float64[]
    end

    if weights === nothing
        weight_array = ones(Float64, size(edge_array, 1))
    else
        weight_array = Float64.(collect(weights))
        length(weight_array) == size(edge_array, 1) || throw(ArgumentError("weights must match edges"))
    end

    acc = Dict{Tuple{Int32,Int32},Float64}()
    for i in axes(edge_array, 1)
        u = edge_array[i, 1]
        v = edge_array[i, 2]
        if drop_loops && u == v
            continue
        end
        a = min(u, v)
        b = max(u, v)
        acc[(a, b)] = get(acc, (a, b), 0.0) + weight_array[i]
    end

    if isempty(acc)
        return _empty_edges(), Float64[]
    end

    keys_sorted = sort!(collect(keys(acc)))
    out_edges = Matrix{Int32}(undef, length(keys_sorted), 2)
    out_weights = Vector{Float64}(undef, length(keys_sorted))
    for (i, (u, v)) in enumerate(keys_sorted)
        out_edges[i, 1] = u
        out_edges[i, 2] = v
        out_weights[i] = acc[(u, v)]
    end
    return out_edges, out_weights
end

function build_weighted_graph(
    num_vertices::Integer,
    edges;
    weights=nothing,
    drop_loops::Bool=true,
)::WeightedGraph
    n = Int(num_vertices)
    n >= 0 || throw(ArgumentError("num_vertices must be nonnegative"))

    collapsed_edges, collapsed_weights = collapse_undirected_edges(edges; weights=weights, drop_loops=drop_loops)
    if size(collapsed_edges, 1) == 0
        adjacency = spzeros(Float64, n, n)
        laplacian = spzeros(Float64, n, n)
        degrees = zeros(Float64, n)
        return WeightedGraph(n, collapsed_edges, collapsed_weights, adjacency, laplacian, degrees)
    end

    minimum(collapsed_edges) >= 0 || throw(ArgumentError("edge endpoints must be >= 0"))
    maximum(collapsed_edges) < n || throw(ArgumentError("edge endpoints must lie in [0, num_vertices)"))
    all(collapsed_weights .>= 0.0) || throw(ArgumentError("weights must be nonnegative"))

    u = Int.(collapsed_edges[:, 1]) .+ 1
    v = Int.(collapsed_edges[:, 2]) .+ 1
    w = collapsed_weights

    rows = vcat(u, v)
    cols = vcat(v, u)
    data = vcat(w, w)
    adjacency = sparse(rows, cols, data, n, n)
    degrees = vec(sum(adjacency, dims=2))
    laplacian = spdiagm(0 => degrees) - adjacency
    return WeightedGraph(n, collapsed_edges, collapsed_weights, adjacency, laplacian, degrees)
end

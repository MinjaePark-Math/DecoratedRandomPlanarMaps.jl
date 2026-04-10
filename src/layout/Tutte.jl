# Tutte-style 2D embeddings based on sparse Laplacian solves.

mutable struct DirichletSolveInfo
    solver::String
    fixed_vertices::Int
    free_vertices::Int
    iterations::Union{Nothing,Vector{Int}}
    source_vertex::Union{Nothing,Int}
end

DirichletSolveInfo(; solver::String, fixed_vertices::Int, free_vertices::Int, iterations=nothing, source_vertex=nothing) =
    DirichletSolveInfo(solver, fixed_vertices, free_vertices, iterations, source_vertex)

function _as_2d_rhs(rhs)
    arr = Float64.(rhs)
    if ndims(arr) == 1
        return reshape(arr, :, 1)
    elseif ndims(arr) == 2
        return arr
    else
        throw(ArgumentError("rhs must be a vector or a matrix"))
    end
end

function _jacobi_preconditioner(matrix::SparseMatrixCSC{Float64,Int})
    diagv = Array(diag(matrix))
    inv = zeros(Float64, length(diagv))
    for i in eachindex(diagv)
        if abs(diagv[i]) > 0
            inv[i] = 1.0 / diagv[i]
        end
    end
    return inv
end

function _solve_with_cg(
    matrix::SparseMatrixCSC{Float64,Int},
    rhs;
    tol::Float64,
    maxiter::Union{Nothing,Int},
    preconditioner=nothing,
)
    rhs_2d = _as_2d_rhs(rhs)
    out = zeros(Float64, size(rhs_2d))
    infos = Int[]
    n = size(matrix, 1)
    maxiter_eff = something(maxiter, max(10 * max(n, 1), 100))

    for j in 1:size(rhs_2d, 2)
        b = rhs_2d[:, j]
        x = zeros(Float64, n)
        r = copy(b)
        z = preconditioner === nothing ? copy(r) : preconditioner .* r
        p = copy(z)
        rz_old = dot(r, z)
        bnorm = max(norm(b), 1.0)
        info = 0

        if norm(r) <= tol * bnorm
            out[:, j] .= x
            push!(infos, 0)
            continue
        end

        converged = false
        for iter in 1:maxiter_eff
            Ap = matrix * p
            denom = dot(p, Ap)
            if !(denom > 0)
                info = -1
                break
            end

            α = rz_old / denom
            x .+= α .* p
            r .-= α .* Ap

            if norm(r) <= tol * bnorm
                info = iter
                converged = true
                break
            end

            z = preconditioner === nothing ? copy(r) : preconditioner .* r
            rz_new = dot(r, z)
            β = rz_new / rz_old
            p .= z .+ β .* p
            rz_old = rz_new
            info = iter
        end

        if info < 0
            error("conjugate-gradient solve failed")
        elseif !converged
            x .= matrix \ b
        end

        out[:, j] .= x
        push!(infos, info)
    end

    return out, infos
end

function _solve_spd_system(
    matrix::SparseMatrixCSC{Float64,Int},
    rhs;
    solver::AbstractString,
    tol::Float64,
    maxiter::Union{Nothing,Int},
    solver_options=Dict{String,Any}(),
)
    rhs_2d = _as_2d_rhs(rhs)
    solver_name = lowercase(strip(String(solver)))
    n = size(matrix, 1)

    if n == 0
        return zeros(Float64, size(rhs_2d)), DirichletSolveInfo(solver="none", fixed_vertices=0, free_vertices=0, iterations=Int[])
    end

    direct_threshold = Int(get(solver_options, "direct_threshold", 12_000))
    if solver_name in ("direct", "spsolve")
        out = Matrix{Float64}(matrix \ rhs_2d)
        return out, DirichletSolveInfo(solver="direct", fixed_vertices=0, free_vertices=n, iterations=fill(0, size(rhs_2d, 2)))
    elseif solver_name == "auto" && n <= direct_threshold
        out = Matrix{Float64}(matrix \ rhs_2d)
        return out, DirichletSolveInfo(solver="direct", fixed_vertices=0, free_vertices=n, iterations=fill(0, size(rhs_2d, 2)))
    elseif solver_name in ("amg", "pyamg")
        throw(ArgumentError("AMG is not enabled in this Julia port; use solver=\"cg\" or solver=\"auto\" instead"))
    elseif solver_name ∉ ("auto", "cg")
        throw(ArgumentError("unknown solver $(repr(solver))"))
    end

    preconditioner = _jacobi_preconditioner(matrix)
    out, infos = _solve_with_cg(matrix, rhs_2d; tol=tol, maxiter=maxiter, preconditioner=preconditioner)
    return out, DirichletSolveInfo(solver="cg+jacobi", fixed_vertices=0, free_vertices=n, iterations=infos)
end

function solve_dirichlet_laplacian(
    laplacian::SparseMatrixCSC{Float64,Int},
    fixed_vertices,
    fixed_values;
    solver::AbstractString="auto",
    tol::Float64=1.0e-9,
    maxiter::Union{Nothing,Int}=nothing,
    solver_options=Dict{String,Any}(),
)
    size(laplacian, 1) == size(laplacian, 2) || throw(ArgumentError("laplacian must be square"))
    n = size(laplacian, 1)

    fixed = Int.(collect(fixed_vertices))
    !isempty(fixed) || throw(ArgumentError("at least one fixed vertex is required"))
    minimum(fixed) >= 0 || throw(ArgumentError("fixed vertex indices must be >= 0"))
    maximum(fixed) < n || throw(ArgumentError("fixed vertex indices must lie in [0, n)"))
    length(unique(fixed)) == length(fixed) || throw(ArgumentError("fixed_vertices must be distinct"))

    values = _as_2d_rhs(fixed_values)
    size(values, 1) == length(fixed) || throw(ArgumentError("fixed_values must have one row per fixed vertex"))

    full = zeros(Float64, n, size(values, 2))
    full[fixed .+ 1, :] .= values

    free_mask = trues(n)
    free_mask[fixed .+ 1] .= false
    free = findall(free_mask) .- 1

    if isempty(free)
        info = DirichletSolveInfo(solver="none", fixed_vertices=length(fixed), free_vertices=0, iterations=fill(0, size(values, 2)))
        return full, info
    end

    free1 = Int.(free .+ 1)
    fixed1 = fixed .+ 1
    reduced = laplacian[free1, free1]
    rhs = -(laplacian[free1, fixed1] * values)

    free_values, info = _solve_spd_system(reduced, rhs; solver=solver, tol=tol, maxiter=maxiter, solver_options=solver_options)
    full[free1, :] .= free_values
    info.fixed_vertices = length(fixed)
    info.free_vertices = length(free)
    return full, info
end

function _pick_source_vertex(num_vertices::Integer, boundary_vertices; source_vertex=nothing, seed=nothing)
    boundary = Int.(collect(boundary_vertices))
    boundary_mask = falses(Int(num_vertices))
    boundary_mask[boundary .+ 1] .= true
    interior = findall(.!boundary_mask) .- 1
    isempty(interior) && throw(ArgumentError("no interior vertex is available for harmonic-measure boundary placement"))
    if source_vertex === nothing
        rng = seed === nothing ? Random.default_rng() : MersenneTwister(Int(seed))
        return rand(rng, interior)
    else
        src = Int(source_vertex)
        0 <= src < Int(num_vertices) || throw(ArgumentError("source_vertex must be an interior vertex"))
        !boundary_mask[src + 1] || throw(ArgumentError("source_vertex must be an interior vertex"))
        return src
    end
end

function harmonic_measure_boundary_positions(
    graph::WeightedGraph,
    boundary_vertices;
    source_vertex=nothing,
    seed=nothing,
    radius=nothing,
    solver::AbstractString="auto",
    tol::Float64=1.0e-9,
    maxiter::Union{Nothing,Int}=nothing,
    solver_options=Dict{String,Any}(),
)
    boundary = Int.(collect(boundary_vertices))
    !isempty(boundary) || throw(ArgumentError("boundary_vertices may not be empty"))
    length(unique(boundary)) == length(boundary) || throw(ArgumentError("boundary_vertices must be distinct"))

    n = graph.num_vertices
    src = _pick_source_vertex(n, boundary; source_vertex=source_vertex, seed=seed)

    fixed_vertices = vcat(boundary, [src])
    fixed_values = zeros(Float64, length(fixed_vertices), 1)
    fixed_values[end, 1] = 1.0

    psi, info = solve_dirichlet_laplacian(
        graph.laplacian,
        fixed_vertices,
        fixed_values;
        solver=solver,
        tol=tol,
        maxiter=maxiter,
        solver_options=solver_options,
    )
    psi = psi[:, 1]
    info.source_vertex = src

    mask = trues(n)
    mask[boundary .+ 1] .= false
    interior1 = findall(mask)
    increments = Vector{Float64}(graph.adjacency[boundary .+ 1, interior1] * psi[interior1])
    increments .= max.(increments, 0.0)
    total = sum(increments)

    radius_value = radius === nothing ? sqrt(max(n, 1)) : float(radius)
    radius_value > 0 || throw(ArgumentError("radius must be positive"))

    angles = if total <= 0.0
        range(0.0, 2π; length=length(boundary)+1)[1:end-1]
    else
        cumulative = cumsum(increments)
        start_cdf = vcat(0.0, cumulative[1:end-1])
        2π .* start_cdf ./ total
    end

    positions = zeros(Float64, length(boundary), 2)
    positions[:, 1] .= radius_value .* cos.(angles)
    positions[:, 2] .= radius_value .* sin.(angles)
    return positions, info
end

function compute_tutte_layout(
    num_vertices::Integer,
    edges,
    boundary_vertices;
    boundary_positions=nothing,
    weights=nothing,
    source_vertex=nothing,
    seed=nothing,
    radius=nothing,
    solver::AbstractString="auto",
    tol::Float64=1.0e-9,
    maxiter::Union{Nothing,Int}=nothing,
    solver_options=Dict{String,Any}(),
    return_metadata::Bool=false,
)
    graph = build_weighted_graph(num_vertices, edges; weights=weights, drop_loops=true)
    boundary = Int.(collect(boundary_vertices))
    !isempty(boundary) || throw(ArgumentError("boundary_vertices may not be empty"))
    minimum(boundary) >= 0 || throw(ArgumentError("boundary vertices must be >= 0"))
    maximum(boundary) < graph.num_vertices || throw(ArgumentError("boundary vertices must lie in [0, num_vertices)"))
    length(unique(boundary)) == length(boundary) || throw(ArgumentError("boundary vertices must be distinct"))

    metadata = Dict{String,Any}(
        "solver_requested" => String(solver),
        "num_vertices" => graph.num_vertices,
        "num_edges" => size(graph.edges, 1),
        "boundary_vertices" => collect(boundary),
    )

    if boundary_positions === nothing
        boundary_xy, boundary_info = harmonic_measure_boundary_positions(
            graph,
            boundary;
            source_vertex=source_vertex,
            seed=seed,
            radius=radius,
            solver=solver,
            tol=tol,
            maxiter=maxiter,
            solver_options=solver_options,
        )
        metadata["boundary_info"] = Dict(
            "solver" => boundary_info.solver,
            "source_vertex" => boundary_info.source_vertex,
            "iterations" => boundary_info.iterations,
        )
        metadata["boundary_positions_mode"] = "harmonic_measure"
    else
        boundary_xy = Float64.(boundary_positions)
        size(boundary_xy) == (length(boundary), 2) || throw(ArgumentError("boundary_positions must have shape (len(boundary_vertices), 2)"))
        metadata["boundary_positions_mode"] = "explicit"
    end

    full_pos, embedding_info = solve_dirichlet_laplacian(
        graph.laplacian,
        boundary,
        boundary_xy;
        solver=solver,
        tol=tol,
        maxiter=maxiter,
        solver_options=solver_options,
    )
    positions = Float64.(full_pos)
    metadata["embedding_info"] = Dict(
        "solver" => embedding_info.solver,
        "iterations" => embedding_info.iterations,
        "fixed_vertices" => embedding_info.fixed_vertices,
        "free_vertices" => embedding_info.free_vertices,
    )

    return return_metadata ? (positions, metadata) : positions
end

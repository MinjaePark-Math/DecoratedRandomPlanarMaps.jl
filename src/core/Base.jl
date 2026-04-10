# Shared interface for concrete random planar map models.

abstract type AbstractRandomPlanarMap end

"""
    num_vertices(map)

Return the number of vertices in a concrete random-map model.
"""
function num_vertices(::AbstractRandomPlanarMap)
    throw(MethodError(num_vertices, ()))
end

"""
    num_faces(map)

Return the number of combinatorial faces represented by a concrete model.
"""
function num_faces(::AbstractRandomPlanarMap)
    throw(MethodError(num_faces, ()))
end

"""
    layout_edges(map; drop_loops=true)

Return a 0-based `(E, 2)` edge array used by the layout engine.
"""
function layout_edges(::AbstractRandomPlanarMap; drop_loops::Bool=true)
    throw(MethodError(layout_edges, ()))
end

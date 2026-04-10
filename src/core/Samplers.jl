# Exact samplers for Dyck paths and related combinatorial objects.

function sample_uniform_dyck_path(m::Integer, rng::AbstractRNG=Random.default_rng())::Vector{Int8}
    m_int = Int(m)
    m_int >= 0 || throw(ArgumentError("m must be nonnegative"))
    m_int == 0 && return Int8[]

    up_left = m_int
    down_left = m_int
    steps = Vector{Int8}()
    sizehint!(steps, 2m_int)

    for _ in 1:(2m_int)
        total = (down_left + up_left) * (down_left - up_left + 1)
        down_threshold = (down_left + 1) * (down_left - up_left)
        if down_threshold > 0 && rand(rng, 0:(total - 1)) < down_threshold
            push!(steps, Int8(-1))
            down_left -= 1
        else
            push!(steps, Int8(1))
            up_left -= 1
        end
    end
    return steps
end

function sample_uniform_primitive_dyck_path(m::Integer, rng::AbstractRNG=Random.default_rng())::Vector{Int8}
    m_int = Int(m)
    m_int >= 0 || throw(ArgumentError("m must be nonnegative"))
    m_int == 0 && return Int8[]
    m_int == 1 && return Int8[1, -1]
    return vcat(Int8[1], sample_uniform_dyck_path(m_int - 1, rng), Int8[-1])
end

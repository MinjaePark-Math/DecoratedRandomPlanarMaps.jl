# Small helpers for recording wall-clock timings in programmatic pipelines.

struct StepTiming
    name::String
    seconds::Float64
end

mutable struct TimingRecorder
    steps::Vector{StepTiming}
end

TimingRecorder() = TimingRecorder(StepTiming[])

function track!(f::Function, recorder::TimingRecorder, name::AbstractString)
    t0 = time_ns()
    value = f()
    elapsed = (time_ns() - t0) / 1.0e9
    push!(recorder.steps, StepTiming(String(name), elapsed))
    return value
end

as_dict(recorder::TimingRecorder) = Dict(step.name => step.seconds for step in recorder.steps)
total_seconds(recorder::TimingRecorder) = sum(step.seconds for step in recorder.steps)

function pretty_lines(recorder::TimingRecorder)
    lines = ["- $(step.name): $(round(step.seconds; digits=3))s" for step in recorder.steps]
    push!(lines, "- total: $(round(total_seconds(recorder); digits=3))s")
    return lines
end

function write_json(recorder::TimingRecorder, path)
    out_path = abspath(String(path))
    mkpath(dirname(out_path))
    payload = Dict(
        "steps" => [Dict("name" => step.name, "seconds" => step.seconds) for step in recorder.steps],
        "total_seconds" => total_seconds(recorder),
    )
    open(out_path, "w") do io
        print(io, JSON3.write(payload; allow_inf=true))
    end
    return out_path
end

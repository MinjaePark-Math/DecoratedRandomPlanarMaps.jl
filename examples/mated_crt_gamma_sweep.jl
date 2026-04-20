using Base.Threads
using Dates
using DecoratedRandomPlanarMaps
using JSON3
using Printf

const DEFAULT_GAMMAS = [round(0.2 * i; digits=1) for i in 1:9]
const DEFAULT_OUTPUT_DIR = joinpath(@__DIR__, "out", "mated_crt_sphere_gamma_sweep")

const USAGE = """
Usage:
  julia --project=. examples/mated_crt_gamma_sweep.jl [options]

Options:
  --vertices N          Number of vertices per map (default: 200)
  --seed N              Shared seed for map generation and layout (default: 712)
  --output PATH         Output directory for the sweep export
  --gammas LIST         Comma-separated gamma values (default: 0.2,0.4,...,1.8)
  --sphere-sampler MODE Sphere sampler: exact or approx (default: exact)
  --sphere-tail N       Exact tail cutoff for approx sphere sampler (default: 64)
  --iterations N        Optional SFDP iteration count
  --K X                 Optional SFDP spring constant
  --repulsiveforce X    Optional SFDP repulsive force
  --overlap MODE        Optional SFDP overlap mode
  --web-scale X         Scale factor for web export coordinates (default: 10.0)
  --serial              Disable threaded generation even if Julia has multiple threads
  --help                Show this help message

Examples:
  julia --project=. examples/mated_crt_gamma_sweep.jl
  julia --project=. examples/mated_crt_gamma_sweep.jl --vertices 1000 --output demo/mated_crt_1000
  julia --project=. examples/mated_crt_gamma_sweep.jl --sphere-sampler approx --sphere-tail 96
  julia --project=. examples/mated_crt_gamma_sweep.jl --gammas 0.4,0.8,1.2,1.6 --iterations 500
"""

function _consume_option(args::Vector{String}, i::Int)
    arg = args[i]
    if occursin('=', arg)
        _, value = split(arg, '='; limit=2)
        return value, i
    end
    i < length(args) || error("missing value for option $(arg)")
    return args[i + 1], i + 1
end

function _parse_gamma_list(raw)::Vector{Float64}
    parts = split(strip(raw), ','; keepempty=false)
    isempty(parts) && error("gamma list must not be empty")
    gammas = Float64[]
    for part in parts
        γ = parse(Float64, strip(part))
        0.0 < γ < 2.0 || error("gamma values must lie in (0, 2); got $(γ)")
        push!(gammas, round(γ; digits=6))
    end
    sort!(unique!(gammas))
    return gammas
end

function _parse_cli(args::Vector{String})
    cfg = Dict{String,Any}(
        "vertices" => 200,
        "seed" => 712,
        "output" => DEFAULT_OUTPUT_DIR,
        "gammas" => DEFAULT_GAMMAS,
        "sphere_sampler" => "exact",
        "sphere_exact_tail_cutoff" => 64,
        "iterations" => nothing,
        "K" => nothing,
        "repulsiveforce" => nothing,
        "overlap" => nothing,
        "web_scale" => 10.0,
        "threaded" => true,
    )

    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--help" || arg == "-h"
            println(USAGE)
            exit(0)
        elseif startswith(arg, "--vertices")
            value, i = _consume_option(args, i)
            cfg["vertices"] = parse(Int, value)
        elseif startswith(arg, "--seed")
            value, i = _consume_option(args, i)
            cfg["seed"] = parse(Int, value)
        elseif startswith(arg, "--output")
            value, i = _consume_option(args, i)
            cfg["output"] = value
        elseif startswith(arg, "--gammas")
            value, i = _consume_option(args, i)
            cfg["gammas"] = _parse_gamma_list(value)
        elseif startswith(arg, "--sphere-sampler")
            value, i = _consume_option(args, i)
            cfg["sphere_sampler"] = lowercase(strip(value))
        elseif startswith(arg, "--sphere-tail")
            value, i = _consume_option(args, i)
            cfg["sphere_exact_tail_cutoff"] = parse(Int, value)
        elseif startswith(arg, "--iterations")
            value, i = _consume_option(args, i)
            cfg["iterations"] = parse(Int, value)
        elseif startswith(arg, "--K")
            value, i = _consume_option(args, i)
            cfg["K"] = parse(Float64, value)
        elseif startswith(arg, "--repulsiveforce")
            value, i = _consume_option(args, i)
            cfg["repulsiveforce"] = parse(Float64, value)
        elseif startswith(arg, "--overlap")
            value, i = _consume_option(args, i)
            cfg["overlap"] = strip(value)
        elseif startswith(arg, "--web-scale")
            value, i = _consume_option(args, i)
            cfg["web_scale"] = parse(Float64, value)
        elseif arg == "--serial"
            cfg["threaded"] = false
        else
            error("unknown option $(arg)\n\n$(USAGE)")
        end
        i += 1
    end

    Int(cfg["vertices"]) >= 4 || error("vertices must be at least 4")
    cfg["sphere_sampler"] in ("exact", "approx") || error("sphere sampler must be `exact` or `approx`")
    return cfg
end

_gamma_label(γ::Real) = @sprintf("%.1f", Float64(γ))
_gamma_slug(γ::Real) = replace(_gamma_label(γ), "." => "_")

function _merge_metadata(parts...)
    merged = Dict{String,Any}()
    for part in parts
        part === nothing && continue
        for (k, v) in pairs(part)
            merged[string(k)] = v
        end
    end
    return merged
end

function _matter_central_charge(γ::Real)
    gamma = Float64(γ)
    Q = 2.0 / gamma + gamma / 2.0
    c = 25.0 - 6.0 * Q^2
    return Q, c
end

function _layout_metadata(cfg)
    md = Dict{String,Any}(
        "engine" => "sfdp",
        "dimension" => 3,
    )
    cfg["iterations"] !== nothing && (md["iterations"] = Int(cfg["iterations"]))
    cfg["K"] !== nothing && (md["K"] = Float64(cfg["K"]))
    cfg["repulsiveforce"] !== nothing && (md["repulsiveforce"] = Float64(cfg["repulsiveforce"]))
    cfg["overlap"] !== nothing && (md["overlap"] = string(cfg["overlap"]))
    return md
end

function _build_gamma_entry(γ::Float64, cfg)
    vertices = Int(cfg["vertices"])
    seed = Int(cfg["seed"])
    web_scale = Float64(cfg["web_scale"])
    gamma_label = _gamma_label(γ)
    gamma_slug = _gamma_slug(γ)
    output_root = abspath(string(cfg["output"]))
    export_dir = joinpath(output_root, "gamma_" * gamma_slug)

    t0 = time()
    map_data = build_mated_crt_map(
        ;
        vertices=vertices,
        topology=:sphere,
        gamma=γ,
        seed=seed,
        sphere_sampler=cfg["sphere_sampler"],
        sphere_exact_tail_cutoff=cfg["sphere_exact_tail_cutoff"],
    )
    t1 = time()

    problem = prepare_layout_problem(map_data; dimension=3, options=Dict("engine" => "sfdp"))
    t2 = time()

    pos = compute_sfdp_layout(
        problem.num_vertices,
        problem.edges;
        scale=1.0,
        seed=seed,
        K=cfg["K"],
        repulsiveforce=cfg["repulsiveforce"],
        iterations=cfg["iterations"],
        overlap=cfg["overlap"],
        options=Dict("engine" => "sfdp"),
    )
    t3 = time()

    Q, c = _matter_central_charge(γ)
    metadata = _merge_metadata(
        Dict(
            "model" => "mated_crt",
            "topology" => "sphere",
            "vertices" => vertices,
            "seed" => seed,
            "gamma" => map_data.gamma,
            "gamma_prime" => map_data.gamma_prime,
            "kappa" => map_data.kappa,
            "kappa_prime" => map_data.kappa_prime,
            "correlation" => map_data.brownian_correlation,
            "matter_central_charge" => c,
            "liouville_Q" => Q,
            "gamma_label" => gamma_label,
        ),
        problem.metadata,
        _layout_metadata(cfg),
    )

    export_web_binaries(
        problem.render_map_data,
        pos .* web_scale,
        export_dir;
        edge_groups=problem.edge_groups,
        faces=problem.faces,
        triangles=problem.surface_triangles,
        metadata=metadata,
    )
    t4 = time()

    return Dict{String,Any}(
        "gamma" => γ,
        "gamma_label" => gamma_label,
        "gamma_slug" => gamma_slug,
        "viewer_path" => "gamma_$(gamma_slug)/index.html",
        "export_dir" => "gamma_$(gamma_slug)",
        "seed" => seed,
        "vertices" => vertices,
        "sphere_sampler" => string(cfg["sphere_sampler"]),
        "sphere_exact_tail_cutoff" => Int(cfg["sphere_exact_tail_cutoff"]),
        "num_faces" => num_faces(map_data),
        "layout_vertex_count" => problem.num_vertices,
        "layout_edge_count" => size(problem.edges, 1),
        "correlation" => map_data.brownian_correlation,
        "liouville_Q" => Q,
        "matter_central_charge" => c,
        "sampler" => map_data.sampler,
        "sphere_construction" => get(problem.metadata, "sphere_construction", nothing),
        "layout_graph_mode" => get(problem.metadata, "layout_graph_mode", nothing),
        "generate_seconds" => round(t1 - t0; digits=3),
        "prepare_seconds" => round(t2 - t1; digits=3),
        "layout_seconds" => round(t3 - t2; digits=3),
        "export_seconds" => round(t4 - t3; digits=3),
        "total_seconds" => round(t4 - t0; digits=3),
    )
end

function _write_manifest(output_root::AbstractString, manifest::Dict{String,Any})
    open(joinpath(output_root, "manifest.json"), "w") do io
        print(io, JSON3.write(manifest; allow_inf=true))
    end
end

function _dashboard_html(manifest::Dict{String,Any})
    manifest_json = JSON3.write(manifest; allow_inf=true)
    template = raw"""
<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>Mated-CRT Gamma Sweep</title>
  <style>
    :root {
      --bg: #f4f1eb;
      --panel: rgba(255, 252, 247, 0.95);
      --ink: #1f2a2e;
      --muted: #5d6c72;
      --accent: #9d4d2f;
      --accent-soft: #f0d8c8;
      --line: #ddcfc2;
      --shadow: 0 18px 48px rgba(67, 43, 24, 0.12);
    }
    * { box-sizing: border-box; }
    html, body {
      margin: 0;
      min-height: 100%;
      background:
        radial-gradient(circle at top left, rgba(208, 148, 96, 0.22), transparent 28%),
        radial-gradient(circle at bottom right, rgba(91, 134, 140, 0.16), transparent 24%),
        linear-gradient(180deg, #f8f5ef 0%, var(--bg) 100%);
      color: var(--ink);
      font-family: Georgia, "Times New Roman", serif;
    }
    body {
      display: grid;
      grid-template-columns: minmax(280px, 360px) minmax(0, 1fr);
      gap: 20px;
      padding: 20px;
    }
    aside, main {
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 22px;
      box-shadow: var(--shadow);
      min-height: calc(100vh - 40px);
    }
    aside {
      padding: 22px 20px;
      display: flex;
      flex-direction: column;
      gap: 18px;
    }
    main {
      padding: 16px;
      display: grid;
      grid-template-rows: auto auto minmax(0, 1fr);
      gap: 14px;
      min-width: 0;
    }
    h1 {
      margin: 0;
      font-size: 1.75rem;
      line-height: 1.05;
      letter-spacing: -0.03em;
    }
    .subtitle {
      color: var(--muted);
      font: 14px/1.5 "Trebuchet MS", "Segoe UI", sans-serif;
    }
    .formula, .note, .meta-grid, .summary-grid, .controls {
      font: 14px/1.45 "Trebuchet MS", "Segoe UI", sans-serif;
    }
    .formula, .summary-grid, .entry-meta, .note {
      background: rgba(255, 255, 255, 0.72);
      border: 1px solid var(--line);
      border-radius: 16px;
      padding: 14px 15px;
    }
    .summary-grid, .entry-meta {
      display: grid;
      grid-template-columns: repeat(2, minmax(0, 1fr));
      gap: 10px 14px;
    }
    .summary-grid div, .entry-meta div {
      min-width: 0;
    }
    .summary-grid strong, .entry-meta strong {
      display: block;
      font-size: 0.78rem;
      letter-spacing: 0.08em;
      text-transform: uppercase;
      color: var(--muted);
      margin-bottom: 4px;
    }
    #entry-list {
      display: grid;
      gap: 10px;
      overflow: auto;
      padding-right: 4px;
    }
    .entry-button {
      width: 100%;
      border: 1px solid var(--line);
      background: rgba(255, 255, 255, 0.82);
      border-radius: 16px;
      padding: 12px 14px;
      text-align: left;
      cursor: pointer;
      transition: transform 120ms ease, border-color 120ms ease, background 120ms ease;
      font: 14px/1.4 "Trebuchet MS", "Segoe UI", sans-serif;
      color: var(--ink);
    }
    .entry-button:hover {
      transform: translateY(-1px);
      border-color: var(--accent);
    }
    .entry-button.active {
      border-color: var(--accent);
      background: linear-gradient(135deg, rgba(240, 216, 200, 0.96), rgba(255, 255, 255, 0.96));
    }
    .entry-button .gamma {
      font-size: 1.02rem;
      font-weight: 700;
      letter-spacing: -0.02em;
    }
    .entry-button .charge {
      color: var(--muted);
      margin-top: 2px;
    }
    .controls {
      display: flex;
      flex-wrap: wrap;
      align-items: center;
      justify-content: space-between;
      gap: 10px 16px;
      padding: 0 4px;
    }
    .controls a {
      color: var(--accent);
      text-decoration: none;
      font-weight: 600;
    }
    .viewer-frame {
      width: 100%;
      height: 100%;
      min-height: 60vh;
      border: 1px solid var(--line);
      border-radius: 18px;
      background: #eef2f4;
    }
    .note {
      color: var(--muted);
    }
    code {
      font-family: "Cascadia Mono", Consolas, monospace;
      font-size: 0.92em;
    }
    @media (max-width: 980px) {
      body {
        grid-template-columns: 1fr;
      }
      aside, main {
        min-height: auto;
      }
      .viewer-frame {
        min-height: 70vh;
      }
    }
  </style>
</head>
<body>
  <aside>
    <div>
      <h1>Mated-CRT Sphere Sweep</h1>
      <div class="subtitle">3D SFDP exports across a gamma sweep, bundled into one browser dashboard.</div>
    </div>
    <div class="formula">
      <strong>Recorded Matter Charge</strong><br />
      <span><code>Q = 2 / γ + γ / 2</code></span><br />
      <span><code>c = 25 - 6 Q²</code></span>
    </div>
    <div class="summary-grid" id="summary-grid"></div>
    <div id="entry-list"></div>
    <div class="note">
      If the embedded viewer cannot fetch its binary files from <code>file://</code>, serve this directory over HTTP first.
    </div>
  </aside>
  <main>
    <div class="controls">
      <div class="subtitle" id="selection-label">Loading sweep…</div>
      <a id="open-link" href="#" target="_blank" rel="noreferrer">Open selected viewer</a>
    </div>
    <div class="entry-meta" id="entry-meta"></div>
    <iframe id="viewer-frame" class="viewer-frame" title="Mated-CRT gamma viewer"></iframe>
  </main>

  <script id="manifest-data" type="application/json">__MANIFEST_JSON__</script>
  <script>
    const manifest = JSON.parse(document.getElementById('manifest-data').textContent);
    const entries = manifest.entries || [];
    const entryList = document.getElementById('entry-list');
    const summaryGrid = document.getElementById('summary-grid');
    const entryMeta = document.getElementById('entry-meta');
    const viewerFrame = document.getElementById('viewer-frame');
    const selectionLabel = document.getElementById('selection-label');
    const openLink = document.getElementById('open-link');

    const fmt = (value, digits = 3) => {
      if (value === null || value === undefined || !Number.isFinite(Number(value))) return 'n/a';
      return Number(value).toFixed(digits).replace(/\\.0+$/, '').replace(/(\\.\\d*[1-9])0+$/, '$1');
    };

    summaryGrid.innerHTML = `
      <div><strong>Vertices</strong>${manifest.vertices}</div>
      <div><strong>Seed</strong>${manifest.seed}</div>
      <div><strong>Gammas</strong>${entries.length}</div>
      <div><strong>Threads</strong>${manifest.thread_count}</div>
      <div><strong>Sphere Sampler</strong>${manifest.sphere_sampler}</div>
      <div><strong>Approx Tail</strong>${manifest.sphere_exact_tail_cutoff}</div>
      <div><strong>Layout</strong>${manifest.layout_engine.toUpperCase()} ${manifest.dimension}D</div>
      <div><strong>Generated</strong>${manifest.generated_at}</div>
    `;

    const defaultIndex = Number.isInteger(manifest.default_index) ? manifest.default_index : 0;
    let activeIndex = Math.min(Math.max(defaultIndex, 0), Math.max(entries.length - 1, 0));

    function renderEntryButtons() {
      entryList.innerHTML = '';
      entries.forEach((entry, idx) => {
        const button = document.createElement('button');
        button.type = 'button';
        button.className = 'entry-button' + (idx === activeIndex ? ' active' : '');
        button.innerHTML = `
          <div class="gamma">γ = ${entry.gamma_label}</div>
          <div class="charge">c = ${fmt(entry.matter_central_charge, 4)}</div>
        `;
        button.addEventListener('click', () => {
          activeIndex = idx;
          updateSelection();
          renderEntryButtons();
        });
        entryList.appendChild(button);
      });
    }

    function updateSelection() {
      const entry = entries[activeIndex];
      if (!entry) {
        selectionLabel.textContent = 'No entries in manifest';
        entryMeta.innerHTML = '';
        viewerFrame.removeAttribute('src');
        openLink.removeAttribute('href');
        return;
      }

      selectionLabel.textContent = `γ = ${entry.gamma_label}, c = ${fmt(entry.matter_central_charge, 4)}, ρ = ${fmt(entry.correlation, 6)}`;
      entryMeta.innerHTML = `
        <div><strong>Gamma</strong>${entry.gamma_label}</div>
        <div><strong>Matter Charge</strong>${fmt(entry.matter_central_charge, 6)}</div>
        <div><strong>Liouville Q</strong>${fmt(entry.liouville_Q, 6)}</div>
        <div><strong>Correlation</strong>${fmt(entry.correlation, 6)}</div>
        <div><strong>Sampler</strong>${entry.sampler || 'n/a'}</div>
        <div><strong>Construction</strong>${entry.sphere_construction || 'n/a'}</div>
        <div><strong>Layout Graph</strong>${entry.layout_graph_mode || 'n/a'}</div>
        <div><strong>Layout Edges</strong>${entry.layout_edge_count}</div>
        <div><strong>Sampler Mode</strong>${entry.sphere_sampler || 'n/a'}</div>
        <div><strong>Approx Tail</strong>${entry.sphere_exact_tail_cutoff ?? 'n/a'}</div>
        <div><strong>Generate Time</strong>${fmt(entry.generate_seconds, 3)} s</div>
        <div><strong>Layout Time</strong>${fmt(entry.layout_seconds, 3)} s</div>
      `;
      viewerFrame.src = entry.viewer_path;
      openLink.href = entry.viewer_path;
      window.location.hash = `gamma=${encodeURIComponent(entry.gamma_label)}`;
    }

    const requested = window.location.hash.startsWith('#gamma=')
      ? decodeURIComponent(window.location.hash.slice(7))
      : null;
    if (requested) {
      const idx = entries.findIndex(entry => entry.gamma_label === requested);
      if (idx >= 0) activeIndex = idx;
    }

    renderEntryButtons();
    updateSelection();
  </script>
</body>
</html>
"""
    return replace(template, "__MANIFEST_JSON__" => manifest_json)
end

function _write_dashboard(output_root::AbstractString, manifest::Dict{String,Any})
    open(joinpath(output_root, "index.html"), "w") do io
        write(io, _dashboard_html(manifest))
    end
end

function _closest_default_index(gammas::Vector{Float64})
    isempty(gammas) && return 0
    target = sqrt(2.0)
    best_idx = 1
    best_dist = abs(gammas[1] - target)
    for idx in 2:length(gammas)
        dist = abs(gammas[idx] - target)
        if dist < best_dist
            best_dist = dist
            best_idx = idx
        end
    end
    return best_idx - 1
end

function main(args=ARGS)
    cfg = _parse_cli(collect(args))
    output_root = abspath(string(cfg["output"]))
    gammas = Float64[Float64(γ) for γ in cfg["gammas"]]
    mkpath(output_root)

    println("Writing gamma sweep to $(output_root)")
    println("   vertices=$(cfg["vertices"]) seed=$(cfg["seed"]) gammas=$(join(_gamma_label.(gammas), ", "))")
    println("   sphere_sampler=$(cfg["sphere_sampler"]) exact_tail_cutoff=$(cfg["sphere_exact_tail_cutoff"])")
    println("   julia_threads=$(nthreads()) threaded=$(cfg["threaded"] && nthreads() > 1)")

    print_lock = ReentrantLock()
    entries = Vector{Dict{String,Any}}(undef, length(gammas))
    errors = Vector{Any}(undef, length(gammas))
    fill!(errors, nothing)

    worker = function (idx::Int)
        γ = gammas[idx]
        label = _gamma_label(γ)
        lock(print_lock) do
            println("[start] gamma=$(label) thread=$(threadid())")
        end
        try
            entry = _build_gamma_entry(γ, cfg)
            entries[idx] = entry
            lock(print_lock) do
                println("[done ] gamma=$(label) total=$(entry["total_seconds"])s export=$(joinpath(output_root, entry["export_dir"]))")
            end
        catch err
            errors[idx] = (err, stacktrace(catch_backtrace()))
            lock(print_lock) do
                println("[fail ] gamma=$(label) :: $(sprint(showerror, err))")
            end
        end
        return nothing
    end

    if cfg["threaded"] && nthreads() > 1 && length(gammas) > 1
        Threads.@threads for idx in eachindex(gammas)
            worker(idx)
        end
    else
        for idx in eachindex(gammas)
            worker(idx)
        end
    end

    failed = findall(!isnothing, errors)
    if !isempty(failed)
        idx = first(failed)
        err, bt = errors[idx]
        showerror(stderr, err, bt)
        println(stderr)
        error("gamma sweep failed for $(length(failed)) value(s)")
    end

    manifest = Dict{String,Any}(
        "title" => "Mated-CRT Sphere Gamma Sweep",
        "generated_at" => string(Dates.now()),
        "vertices" => Int(cfg["vertices"]),
        "seed" => Int(cfg["seed"]),
        "dimension" => 3,
        "layout_engine" => "sfdp",
        "thread_count" => nthreads(),
        "threaded_generation" => cfg["threaded"] && nthreads() > 1,
        "sphere_sampler" => string(cfg["sphere_sampler"]),
        "sphere_exact_tail_cutoff" => Int(cfg["sphere_exact_tail_cutoff"]),
        "web_scale" => Float64(cfg["web_scale"]),
        "gammas" => gammas,
        "default_index" => _closest_default_index(gammas),
        "entries" => entries,
    )

    _write_manifest(output_root, manifest)
    _write_dashboard(output_root, manifest)

    println("Wrote dashboard:")
    println("   $(joinpath(output_root, "index.html"))")
    println("Wrote manifest:")
    println("   $(joinpath(output_root, "manifest.json"))")
end

main()

using Base.Threads
using Dates
using DecoratedRandomPlanarMaps
using JSON3
using Printf

function _gamma_from_matter_central_charge(c::Real)
    charge = Float64(c)
    charge <= 1.0 || error("matter charges must be <= 1; got $(charge)")
    Q = sqrt((25.0 - charge) / 6.0)
    return Q - sqrt(Q^2 - 4.0)
end

const DEFAULT_CHARGES = Float64[1.0, 0.5, 0.0, -1.0, -2.0, -5.0, -8.0, -12.5]
const DEFAULT_GAMMAS = let
    vals = Float64[0.5]
    append!(vals, _gamma_from_matter_central_charge.(DEFAULT_CHARGES))
    sort!(unique!(round.(vals; digits=12)))
end
const DEFAULT_OUTPUT_DIR = joinpath(@__DIR__, "out", "mated_crt_sphere_gamma_sweep")

const USAGE = """
Usage:
  julia --project=. examples/mated_crt_gamma_sweep.jl [options]

Options:
  --vertices N          Number of vertices per map (default: 10000)
  --seed N              Shared seed for map generation and layout (default: 712)
  --output PATH         Output directory for the sweep export
  --gammas LIST         Comma-separated gamma values to include
  --charges LIST        Comma-separated matter charges to include (converted to gamma)
  --sphere-sampler MODE Sphere sampler: exact or approx (default: approx)
  --sphere-tail N       Exact tail cutoff for approx sphere sampler (default: 64)
  --grid-columns N      Number of viewer columns in the dashboard grid (default: 3)
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
  julia --project=. examples/mated_crt_gamma_sweep.jl --sphere-tail 96 --grid-columns 3
  julia --project=. examples/mated_crt_gamma_sweep.jl --charges 1,0,-2,-8,-12.5 --gammas 0.5
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
        0.0 < γ <= 2.0 || error("gamma values must lie in (0, 2]; got $(γ)")
        push!(gammas, round(γ; digits=12))
    end
    sort!(unique!(gammas))
    return gammas
end

function _parse_charge_list(raw)::Vector{Float64}
    parts = split(strip(raw), ','; keepempty=false)
    isempty(parts) && error("charge list must not be empty")
    charges = Float64[]
    for part in parts
        charge = parse(Float64, strip(part))
        charge <= 1.0 || error("matter charges must be <= 1; got $(charge)")
        push!(charges, round(charge; digits=12))
    end
    sort!(unique!(charges); rev=true)
    return charges
end

function _start_custom_value_override!(cfg)
    if !Bool(cfg["custom_values"])
        cfg["gammas"] = Float64[]
        cfg["custom_values"] = true
    end
    return cfg
end

function _parse_cli(args::Vector{String})
    cfg = Dict{String,Any}(
        "vertices" => 10_000,
        "seed" => 712,
        "output" => DEFAULT_OUTPUT_DIR,
        "gammas" => DEFAULT_GAMMAS,
        "sphere_sampler" => "approx",
        "sphere_exact_tail_cutoff" => 64,
        "grid_columns" => 3,
        "custom_values" => false,
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
            _start_custom_value_override!(cfg)
            append!(cfg["gammas"], _parse_gamma_list(value))
        elseif startswith(arg, "--charges")
            value, i = _consume_option(args, i)
            _start_custom_value_override!(cfg)
            append!(cfg["gammas"], _gamma_from_matter_central_charge.(_parse_charge_list(value)))
        elseif startswith(arg, "--sphere-sampler")
            value, i = _consume_option(args, i)
            cfg["sphere_sampler"] = lowercase(strip(value))
        elseif startswith(arg, "--sphere-tail")
            value, i = _consume_option(args, i)
            cfg["sphere_exact_tail_cutoff"] = parse(Int, value)
        elseif startswith(arg, "--grid-columns")
            value, i = _consume_option(args, i)
            cfg["grid_columns"] = parse(Int, value)
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
    Int(cfg["grid_columns"]) >= 1 || error("grid-columns must be at least 1")
    cfg["sphere_sampler"] in ("exact", "approx") || error("sphere sampler must be `exact` or `approx`")
    cfg["gammas"] = sort!(unique!(round.(Float64.(cfg["gammas"]); digits=12)))
    return cfg
end

function _trim_fixed_decimal(x::Real; digits::Int=6)
    s = @sprintf("%.*f", digits, Float64(x))
    s = replace(s, r"(\.\d*?[1-9])0+$" => s"\1")
    s = replace(s, r"\.0+$" => "")
    return s
end

function _approx_same(a::Real, b::Real; atol::Real=1.0e-8)
    return abs(Float64(a) - Float64(b)) <= float(atol)
end

_gamma_label(γ::Real) = _trim_fixed_decimal(γ; digits=6)
_gamma_slug(γ::Real) = replace(replace(_gamma_label(γ), "." => "_"), "-" => "m")

function _gamma_display_label(γ::Real)
    gamma = Float64(γ)
    for (value, label) in (
        (1.0, "1"),
        (sqrt(2.0), "√2"),
        (sqrt(8.0 / 3.0), "√(8/3)"),
        (sqrt(3.0), "√3"),
        (2.0, "2"),
    )
        _approx_same(gamma, value) && return label
    end
    return @sprintf("%.2f", gamma)
end

function _charge_display_label(c::Real)
    charge = Float64(c)
    nearest_int = round(Int, charge)
    if _approx_same(charge, nearest_int)
        return string(nearest_int)
    end
    nearest_half = round(charge * 2.0) / 2.0
    if _approx_same(charge, nearest_half)
        return @sprintf("%.1f", charge)
    end
    return @sprintf("%.2f", charge)
end

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
    gamma_display_label = _gamma_display_label(γ)
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
    charge_display_label = _charge_display_label(c)
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
            "gamma_display_label" => gamma_display_label,
            "matter_central_charge_label" => charge_display_label,
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
        "gamma_display_label" => gamma_display_label,
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
        "matter_central_charge_label" => charge_display_label,
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
  <title>Liouville Quantum Gravity (LQG) Gamma Sweep</title>
  <style>
    :root {
      --bg: #f3eee4;
      --panel: rgba(255, 251, 246, 0.96);
      --ink: #20272b;
      --accent: #9b452a;
      --line: #decfbe;
      --shadow: 0 16px 36px rgba(79, 52, 31, 0.12);
      --grid-cols: 3;
    }
    * { box-sizing: border-box; }
    html, body {
      margin: 0;
      min-height: 100%;
      background:
        radial-gradient(circle at top left, rgba(204, 144, 96, 0.23), transparent 28%),
        radial-gradient(circle at bottom right, rgba(87, 131, 136, 0.14), transparent 26%),
        linear-gradient(180deg, #f8f5ef 0%, var(--bg) 100%);
      color: var(--ink);
      font-family: Georgia, "Times New Roman", serif;
    }
    body {
      padding: 18px;
    }
    .page {
      display: grid;
      grid-template-rows: auto minmax(0, 1fr);
      gap: 16px;
      max-width: 1800px;
      margin: 0 auto;
      min-height: calc(100vh - 36px);
    }
    .page.one-row {
      grid-template-rows: auto auto;
      min-height: 0;
    }
    .page.multi-row {
      grid-template-rows: auto minmax(0, 1fr);
      min-height: calc(100vh - 36px);
    }
    header, .viewer-scroller {
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 22px;
      box-shadow: var(--shadow);
      padding: 18px 20px;
    }
    header {
      display: flex;
      align-items: flex-start;
      justify-content: space-between;
      gap: 16px;
    }
    .header-copy {
      min-width: 0;
    }
    h1 {
      margin: 0;
      font-size: clamp(1.9rem, 3vw, 3rem);
      line-height: 1.02;
      letter-spacing: -0.03em;
    }
    .subtitle {
      color: #5c6970;
      font: 14px/1.5 "Trebuchet MS", "Segoe UI", sans-serif;
      margin: 0;
      margin-top: 6px;
    }
    .header-controls {
      display: flex;
      align-items: center;
      gap: 10px;
      flex-wrap: wrap;
      color: #5c6970;
      font: 13px/1.4 "Trebuchet MS", "Segoe UI", sans-serif;
    }
    .row-picker {
      display: inline-flex;
      gap: 6px;
      align-items: center;
      padding: 6px;
      border: 1px solid var(--line);
      border-radius: 999px;
      background: rgba(255, 255, 255, 0.76);
    }
    .row-label {
      padding: 0 6px 0 4px;
      font-weight: 600;
      color: #5c6970;
    }
    .row-button {
      border: 1px solid transparent;
      border-radius: 999px;
      background: transparent;
      color: var(--ink);
      padding: 7px 12px;
      cursor: pointer;
      font: inherit;
    }
    .row-button.active {
      background: rgba(155, 69, 42, 0.12);
      border-color: rgba(155, 69, 42, 0.24);
      color: var(--accent);
      font-weight: 700;
    }
    .action-button {
      border: 1px solid var(--line);
      border-radius: 999px;
      background: rgba(255, 255, 255, 0.82);
      color: var(--ink);
      padding: 9px 14px;
      cursor: pointer;
      font: inherit;
      font-weight: 600;
    }
    .viewer-scroller {
      min-height: 0;
      padding-bottom: 14px;
      scrollbar-width: auto;
      scrollbar-color: rgba(155, 69, 42, 0.82) rgba(222, 207, 190, 0.9);
    }
    .viewer-scroller::-webkit-scrollbar {
      width: 16px;
      height: 16px;
    }
    .viewer-scroller::-webkit-scrollbar-track {
      background: rgba(222, 207, 190, 0.9);
      border-radius: 999px;
    }
    .viewer-scroller::-webkit-scrollbar-thumb {
      background: linear-gradient(180deg, rgba(155, 69, 42, 0.9), rgba(122, 53, 32, 0.92));
      border-radius: 999px;
      border: 3px solid rgba(222, 207, 190, 0.92);
    }
    .viewer-scroller.one-row {
      align-self: start;
      overflow-x: scroll;
      overflow-y: hidden;
      padding-bottom: 18px;
    }
    .viewer-scroller.multi-row {
      overflow-x: hidden;
      overflow-y: scroll;
      padding-right: 18px;
    }
    .viewer-grid {
      display: grid;
      gap: 16px;
      min-height: 0;
    }
    .viewer-grid.one-row {
      grid-auto-flow: column;
      grid-auto-columns: minmax(180px, 15.5vw);
      width: max-content;
    }
    .viewer-grid.multi-row {
      width: 100%;
    }
    .viewer-card {
      min-width: 0;
      display: grid;
      grid-template-rows: auto minmax(220px, 1fr);
      gap: 10px;
      background: rgba(255, 255, 255, 0.7);
      border: 1px solid var(--line);
      border-radius: 18px;
      padding: 14px;
    }
    .card-head {
      display: flex;
      align-items: baseline;
      justify-content: space-between;
      gap: 12px;
      min-width: 0;
    }
    .card-head .gamma {
      font-size: 1.12rem;
      font-weight: 700;
      letter-spacing: -0.02em;
      color: var(--ink);
    }
    .card-head .charge {
      flex-shrink: 0;
      color: var(--accent);
      font-weight: 700;
    }
    .viewer-frame {
      border: 1px solid var(--line);
      border-radius: 16px;
      width: 100%;
      height: clamp(220px, 26vh, 340px);
      background: #eef2f4;
    }
    @media (max-width: 1080px) {
      .page.multi-row {
        min-height: calc(100vh - 36px);
      }
      .viewer-grid.one-row {
        grid-auto-columns: minmax(170px, 35vw);
      }
    }
    @media (max-width: 760px) {
      body {
        padding: 12px;
      }
      .page.multi-row {
        min-height: calc(100vh - 24px);
      }
      header {
        flex-direction: column;
        align-items: stretch;
      }
      .header-controls {
        justify-content: flex-start;
      }
      .viewer-grid.multi-row {
        grid-template-columns: 1fr !important;
      }
      .viewer-grid.one-row {
        grid-auto-columns: minmax(150px, 44vw);
      }
      .viewer-frame {
        height: 54vh;
      }
    }
  </style>
</head>
<body>
  <div class="page" id="page">
    <header>
      <div class="header-copy">
        <h1>Liouville Quantum Gravity (LQG) Surfaces for Various γ</h1>
        <div class="subtitle">Mated-CRT sphere approximations across the gamma sweep.</div>
      </div>
      <div class="header-controls">
        <div class="row-picker" id="row-picker">
          <span class="row-label">Rows</span>
          <button type="button" class="row-button" data-rows="1">1</button>
          <button type="button" class="row-button" data-rows="2">2</button>
          <button type="button" class="row-button" data-rows="3">3</button>
        </div>
        <button type="button" class="action-button" id="motion-button">Pause Motion</button>
        <button type="button" class="action-button" id="reset-button">Reset Views</button>
      </div>
    </header>

    <section class="viewer-scroller" id="viewer-scroller">
      <div class="viewer-grid" id="viewer-grid"></div>
    </section>
  </div>

  <script id="manifest-data" type="application/json">__MANIFEST_JSON__</script>
  <script>
    const manifest = JSON.parse(document.getElementById('manifest-data').textContent);
    const entries = manifest.entries || [];
    const page = document.getElementById('page');
    const viewerScroller = document.getElementById('viewer-scroller');
    const viewerGrid = document.getElementById('viewer-grid');
    const rowButtons = Array.from(document.querySelectorAll('.row-button'));
    const motionButton = document.getElementById('motion-button');
    const resetButton = document.getElementById('reset-button');
    const viewerFrames = [];
    let motionEnabled = true;

    const fmt = (value, digits = 3) => {
      if (value === null || value === undefined || !Number.isFinite(Number(value))) return 'n/a';
      return Number(value).toFixed(digits).replace(/\\.0+$/, '').replace(/(\\.\\d*[1-9])0+$/, '$1');
    };

    function broadcastViewerCommand(action, value = null) {
      const message = { kind: 'randmaps-embed-control', action };
      if (value !== null) {
        message.value = value;
      }
      viewerFrames.forEach((frame) => {
        try {
          frame.contentWindow?.postMessage(message, '*');
        } catch (_) {
        }
      });
    }

    function syncMotionButton() {
      motionButton.textContent = motionEnabled ? 'Pause Motion' : 'Play Motion';
    }

    function clampRows(value) {
      const parsed = Number(value);
      return parsed === 1 || parsed === 2 || parsed === 3 ? parsed : 3;
    }

    function initialRows() {
      return 2;
    }

    function applyRowLayout(rows) {
      const rowCount = clampRows(rows);
      const cols = Math.max(1, Math.ceil(Math.max(entries.length, 1) / rowCount));

      page.classList.toggle('one-row', rowCount === 1);
      page.classList.toggle('multi-row', rowCount !== 1);
      viewerScroller.classList.toggle('one-row', rowCount === 1);
      viewerScroller.classList.toggle('multi-row', rowCount !== 1);
      viewerGrid.classList.toggle('one-row', rowCount === 1);
      viewerGrid.classList.toggle('multi-row', rowCount !== 1);

      if (rowCount === 1) {
        viewerGrid.style.gridTemplateColumns = '';
      } else {
        viewerGrid.style.gridTemplateColumns = `repeat(${cols}, minmax(0, 1fr))`;
      }

      rowButtons.forEach((button) => {
        button.classList.toggle('active', Number(button.dataset.rows) === rowCount);
      });
    }

    entries.forEach((entry) => {
      const card = document.createElement('article');
      card.className = 'viewer-card';
      card.innerHTML = `
        <div class="card-head">
          <div class="gamma">γ = ${entry.gamma_display_label || entry.gamma_label}</div>
          <div class="charge">c = ${entry.matter_central_charge_label || fmt(entry.matter_central_charge, 2)}</div>
        </div>
        <iframe class="viewer-frame" title="γ = ${entry.gamma_display_label || entry.gamma_label}" src="${entry.viewer_path}?embed=1&autorotate=1" loading="lazy"></iframe>
      `;
      const frame = card.querySelector('iframe');
      frame.addEventListener('load', () => {
        try {
          frame.contentWindow?.postMessage(
            { kind: 'randmaps-embed-control', action: 'setAutoRotate', value: motionEnabled },
            '*',
          );
        } catch (_) {
        }
      });
      viewerFrames.push(frame);
      viewerGrid.appendChild(card);
    });

    const selectedRows = initialRows();
    applyRowLayout(selectedRows);
    rowButtons.forEach((button) => {
      button.addEventListener('click', () => applyRowLayout(button.dataset.rows));
    });
    syncMotionButton();
    motionButton.addEventListener('click', () => {
      motionEnabled = !motionEnabled;
      syncMotionButton();
      broadcastViewerCommand('setAutoRotate', motionEnabled);
    });
    resetButton.addEventListener('click', () => {
      broadcastViewerCommand('resetView');
    });
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
        "title" => "Liouville Quantum Gravity (LQG) Gamma Sweep",
        "generated_at" => string(Dates.now()),
        "vertices" => Int(cfg["vertices"]),
        "seed" => Int(cfg["seed"]),
        "dimension" => 3,
        "layout_engine" => "sfdp",
        "thread_count" => nthreads(),
        "threaded_generation" => cfg["threaded"] && nthreads() > 1,
        "sphere_sampler" => string(cfg["sphere_sampler"]),
        "sphere_exact_tail_cutoff" => Int(cfg["sphere_exact_tail_cutoff"]),
        "grid_columns" => Int(cfg["grid_columns"]),
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

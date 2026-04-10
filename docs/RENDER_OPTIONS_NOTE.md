# Render options and defaults

## Entry points

### Web export

```julia
export_web_binaries(
    map_data,
    pos,
    output_dir;
    edge_groups=nothing,
    faces=nothing,
    triangles=nothing,
    metadata=nothing,
    write_index_html=true,
)
```

### Makie 2D

```julia
render_makie_2d(
    map_data,
    pos;
    edge_groups=nothing,
    faces=nothing,
    triangles=nothing,
    metadata=nothing,
    title="2D map",
)
```

### Makie 3D

```julia
render_makie_3d(
    map_data,
    pos;
    edge_groups=nothing,
    faces=nothing,
    triangles=nothing,
    metadata=nothing,
    title="3D map",
)
```

## Shared concepts

All renderers use:
- vertex positions `pos` with shape `(N, 2)` or `(N, 3)`
- triangulated surface faces from `triangles` or `faces`
- edge groups from `edge_groups` or from `map_data`
- an optional FK-only exploration overlay

Web, Makie, and SVG previews ignore non-finite vertices when framing the view and skip any face / edge / exploration segment that touches a non-finite endpoint. This is especially important for FK gasket-style 2D layouts, where some helper code paths may still leave NaNs on removed exterior vertices.

## Exploration overlay

Exploration is shown only when the input looks FK-like:
- `map_data isa FKMap`, or
- `metadata["model"] == "fk"` / `"spanning_tree"` / `"hc"`, or
- `edge_groups` looks like FK colors (`generic/green` plus some of `red`, `blue`, `purple`, `orange`)

Exploration is reconstructed from:
- the current triangle set
- the current FK generic-edge set
- midpoint-to-midpoint segments inside triangles that have exactly two generic edges

Non-finite vertices are skipped, so NaNs outside a gasket do not create bad segments.

## `web_meta.json` fields

The web exporter writes:
- `model`
- `dimension`
- `num_vertices`
- `num_edges`
- `edge_groups`
- `files.vertices`
- `files.faces`
- `files.exploration`

Each `edge_groups` entry contains:
- `name`
- `file`
- `color`
- `num_edges`
- `kind`
  - `index_pairs` for usual edge groups
  - `positions` for explicit line endpoint buffers like exploration
- `default_visible`

## Visibility defaults

### Web
- `faces`: off in 2D, on in 3D
- regular edge groups: on
- `exploration`: off

### Makie
- `faces`: on
- Schnyder `generic` in 2D: off
- all other regular edge groups: on
- FK `exploration`: off

## Color conventions

### Uniform
- `generic`: gray

### Schnyder
- `generic`: light gray and thin in 2D
- `orange`: orange
- `navy`: dark blue
- `green`: dark green and bold
- `outer`: dark boundary color

### FK / spanning-tree
- `green` is normalized to `generic`
- `generic`: light gray
- `red`: red
- `blue`: blue
- `purple`: purple
- `orange`: orange
- `exploration`: green

## Likely failure modes

If exploration still does not show:
1. `web_meta.json` does not contain an `edge_groups` entry named `exploration`
2. the input is not being recognized as FK-like
3. the current triangle set has no triangles with exactly two generic edges
4. all candidate midpoint segments touch non-finite vertices and are filtered out


## Web summary and Makie titles

Web exports expose both `parameters` and `parameters_summary` in `web_meta.json`. The bundled viewer shows the parameter summary in the HUD. Makie renderers append the same parameter summary to the title when metadata is available; for default titles this happens automatically.

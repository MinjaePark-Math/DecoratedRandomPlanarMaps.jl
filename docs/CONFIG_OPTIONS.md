# Config options

`run_pipeline` accepts either a nested `Dict` or a YAML file loaded with `load_config`.

## Top-level sections

- `model`
- `layout`
- `output`

## `model`

Common keys:

- `type`: `uniform`, `schnyder`, `fk`, `spanning_tree`
- `faces`: number of faces
- `seed`: RNG seed

FK-only keys:

- `p`: FK parameter in `[0, 1]`
- `q`: alternate FK parameter, converted to `p`
- `sampling_method`: `auto`, `exact_rejection`, or `approx`
- `pool_size`: optional proposal pool size for the approximate sampler; when omitted, a smooth size- and `p`-dependent default is used
- `mh_steps`: optional Metropolis–Hastings step count for the approximate sampler; when omitted, the default is `4 * pool_size`
- `exact_pilot_samples`: short pilot sample count used by `sampling_method: auto`
- `exact_min_acceptance`: minimum estimated acceptance rate required before `auto` chooses exact rejection
- `max_exact_tries`: safety cap for `exact_rejection`

## `layout`

Common keys:

- `dimension`: `2` or `3`
- `seed`: layout RNG seed
- `boundary_scale`: boundary radius / scale for fixed 2D boundary placement
- `options`: additional model-specific options

2D keys:

- `engine`: `tutte` or `circle_packing`
- `solver`: `auto`, `direct`, `spsolve`, or `cg`
- `tol`: iterative tolerance for `cg` or convergence tolerance for circle packing
- `maxiter`: optional max iteration count
- `radius`: optional boundary radius override
- `source_vertex`: optional harmonic-measure source vertex
- `relaxation`: optional circle-packing radius relaxation factor
- `initial_radius`: optional circle-packing starting radius
- `min_radius`: optional lower clamp for circle-packing radii

3D keys:

- `engine`: currently `sfdp`
- `normalize_scale`: target radius after normalization
- `K`: optional SFDP ideal edge-length parameter
- `repulsiveforce`: optional SFDP repulsion-strength parameter
- `iterations`: optional iteration count for the built-in force-layout fallback
- `overlap`: optional Graphviz overlap mode

`layout.options` currently supports:

- `hc_boundary_mode`: `h_gasket`, `c_gasket`, or `face` (`face` uses a triangular outer face with explicit equilateral boundary positions in 2D)
- `outer_face_index`: explicit outer face choice for uniform 2D

## `output`

- `export_web`: folder for Three.js export
- `export_stl`: path for STL export (3D only)
- `preview_svg`: path for a 2D SVG preview
- `preview_web`: folder for a 3D preview web export
- `show`: whether to open the default preview viewer
- `timings_json`: optional JSON file for timing output
- `web_scale`: export scaling for web
- `stl_scale`: export scaling for STL
- `show_scale`: preview scaling for the viewer, and for `preview_svg` when set

## Minimal examples

### Uniform 3D

```yaml
model:
  type: uniform
  faces: 200
  seed: 7

layout:
  dimension: 3
  engine: sfdp
  normalize_scale: 1.0

output:
  export_web: ./exports/uniform_3d
  export_stl: ./exports/uniform_3d.stl
```

### FK 2D h-gasket

```yaml
model:
  type: fk
  faces: 400
  p: 0.25
  seed: 11

layout:
  dimension: 2
  engine: tutte
  solver: auto
  options:
    hc_boundary_mode: h_gasket

output:
  export_web: ./exports/fk_h_gasket
  preview_svg: ./exports/fk_h_gasket/preview.svg
```

### Schnyder 2D circle packing

```yaml
model:
  type: schnyder
  faces: 200
  seed: 7

layout:
  dimension: 2
  engine: circle_packing
  maxiter: 100

output:
  preview_svg: ./exports/schnyder_circle_packing.svg
```


### FK sampling note

`sampling_method: exact_rejection` is exact, but it can become very slow when `faces` is large or `p` is too close to `1`. The default `sampling_method: auto` estimates the rejection acceptance rate first and falls back to the approximate sampler when exact rejection looks too expensive.


### Approximate FK sampler note

When `sampling_method: approx` is used and `pool_size` / `mh_steps` are omitted, the package now chooses them automatically from the map size and `p`:

- base scale `≈ sqrt(faces)`
- smoothly increases with `p`
- reaches about `2 * sqrt(faces)` by `p = 0.5`
- `mh_steps = 4 * pool_size`

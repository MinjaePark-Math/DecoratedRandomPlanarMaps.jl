# Examples

Run from the package root with:

```bash
julia --project=. examples/uniform_examples.jl
```

Swap in any of:

- `examples/uniform_examples.jl`
- `examples/schnyder_examples.jl`
- `examples/fk_examples.jl`
- `examples/spanning_tree_examples.jl`
- `examples/meandric.jl`
- `examples/meandric_arc_diagrams.jl`

Each file includes all relevant export blocks:

- SVG preview
- web export
- STL export for 3D
- optional Makie blocks

Comment out the blocks you do not want.

Each example also prints a simple timing breakdown for:

- map generation
- 2D/3D layout preparation
- Tutte or SFDP layout
- each non-interactive export step

Makie is left commented out and is intentionally not included in the timing breakdown because the interactive window holds the event loop open.

# BioEMMA

BioEMMA is an early-stage Python library for building Escher-compatible
metabolic maps from KEGG pathway layouts and genome-scale metabolic models.

The current main workflow is:

1. Parse a KEGG KGML/XML pathway with `KeggMap`.
2. Convert KEGG compounds and reactions to BiGG/SEED identifiers using bundled
   MetaNetX-derived mapping tables.
3. Build an Escher JSON map with `EscherMapper`.
4. Optionally save a reproducible workflow output directory with the Escher map,
   the reconstructed KEGG map, flux data, summaries, and merged maps.

The project is currently in alpha. The public API may still change while the
package structure is being prepared for PyPI.

## Installation

For local development:

```bash
pip install -e .
```

Runtime dependencies can also be installed from:

```bash
pip install -r requirements.txt
```

## Basic Usage

The workflow API is the recommended user-facing entry point. It accepts a COBRA
model path (or an in-memory `cobra.Model`) and either a KEGG pathway identifier
or a local KGML file.

```python
from bioemma.workflow import build_outputs


result = build_outputs(
    model="path/to/model.xml",
    pathway="rn00010",
    output_dir="out",
    database="BIGG",
    run_fba=True,
)

escher_map = result.escher_map
kegg_reconstruction = result.kegg_reconstruction
```

`escher_map` is a Python object compatible with the Escher JSON map structure,
and `kegg_reconstruction` is a normalized analytical representation of the KEGG
layout and mapped identifiers. When `save_kegg_map=True`, BioEMMA also writes
`kegg_escher_map.json`: a pure KEGG-layout Escher map before model filtering or
secondary metabolite addition.

With `output_dir`, BioEMMA writes:

```text
out/rn00010/
  escher_map.json
  kegg_escher_map.json     # when save_kegg_map=True or --save-kegg-map
  kegg_source_reconstruction.json
  summary.json
  fluxes.json              # when fluxes are provided or run_fba=True
  escher_map.html          # when save_html=True
  escher_map_with_fluxes.html  # when flux data and HTML output are requested
```

HTML output requires the `escher` package. BioEMMA does not export PNG files
directly; open the HTML output in Escher and use Escher's built-in PNG export
when a raster image is needed.

Visualization layout settings can be tuned with `VisualizationOptions`:

```python
from bioemma.workflow import build_outputs
from bioemma.visualization import VisualizationOptions


result = build_outputs(
    model="path/to/model.xml",
    pathway="rn00010",
    output_dir="out",
    visualization_options=VisualizationOptions(
        scaling_factor=4,
        axis_epsilon=2,
        markers_dist=10,
        metabolite_label_shift=(10, 10),
        reaction_label_shift=(10, 10),
        canvas_margin_x=160,
        canvas_margin_y=160,
        axis_offset=20,
    ),
)
```

The defaults are conservative starting values for KEGG layouts: coordinates are
scaled up for Escher readability, aligned reaction lanes keep a small tolerance,
and secondary metabolites get enough spacing after scaling.

## Command Line Usage

Build one map from a KEGG pathway identifier:

```bash
bioemma build --model path/to/model.xml --pathway rn00010 --output-dir out
```

Build one map from a local KGML file:

```bash
bioemma build --model path/to/model.xml --kgml path/to/rn00010.xml --output-dir out
```

Build multiple maps and merge them:

```bash
bioemma build --model path/to/model.xml --pathway rn00010 rn00020 --output-dir out
```

The same works with local KGML files:

```bash
bioemma build --model path/to/model.xml --kgml path/to/rn00010.xml path/to/rn00020.xml --output-dir out
```

For multiple inputs, BioEMMA writes each individual map into its own subfolder
and writes a merged Escher map at:

```text
out/merged_escher_map.json
```

Use `--no-merge` to skip the merged map.

The legacy single-file JSON output is still available:

```bash
bioemma build --model path/to/model.xml --kgml path/to/rn00010.xml --output escher_map.json
```

`summary.json` includes `map_stats`, a stage-by-stage count of total elements,
nodes, reactions, and segments added or removed while the map is built. To print
the same reduction statistics in the CLI, add `--map-stats`:

```bash
bioemma build --model path/to/model.xml --kgml path/to/rn00010.xml --output-dir out --map-stats
```

To save the unfiltered KEGG Escher map next to the normal model-derived map,
add `--save-kegg-map`:

```bash
bioemma build --model path/to/model.xml --kgml path/to/rn00010.xml --output-dir out --save-kegg-map
```

The same visualization settings are available in the CLI, for example:

```bash
bioemma build --model path/to/model.xml --kgml path/to/rn00010.xml --output-dir out --scaling-factor 4 --canvas-margin-x 160 --canvas-margin-y 160
```

If cobrapy cannot access its default cache directory on Windows, set a local
cache directory before running tests or CLI commands:

```cmd
set BIOEMMA_COBRA_CACHE_DIR=%CD%\.cobra-cache
```

## Included Mapping Data

BioEMMA currently bundles two compact runtime mapping files:

- `metabolite_mapping.tsv`
- `reaction_mapping.tsv`

These files are derived from MetaNetX cross-reference tables and are used to
map KEGG identifiers to BiGG and SEED identifiers. The large raw MetaNetX
download cache is not intended to be included in the Python package.

See `NOTICE.md` for third-party data attribution and usage notes.

## License

BioEMMA's source code is distributed under the MIT License. Bundled mapping
data are derived from third-party database resources and may be subject to
their own license terms. See `LICENSE` and `NOTICE.md`.

By default, the workflow keeps the KEGG reactions and compounds that can be
matched to the COBRA model. To preserve KEGG-only elements that are not present
in the model, pass `include_kegg_only=True` in Python or use
`--include-kegg-only` in the CLI.

## Development Notes

The package code lives in:

```text
src/bioemma/
```

The current core modules are:

- `bioemma.maps.KeggMap`
- `bioemma.mapper_base.EscherMapper`
- `bioemma.metanetx_mapper.MetaNetXMapper`
- `bioemma.merger.EscherMerger`
- `bioemma.workflow.build_outputs`
- `bioemma.workflow.build_many_outputs`

The script for regenerating mapping tables is kept separately in:

```text
scripts/prepare_db_mapping.py
```

Run the test suite from a source checkout with:

```cmd
set PYTHONPATH=%CD%\src
set BIOEMMA_COBRA_CACHE_DIR=%CD%\.pytest-cobra-cache
python -m pytest -q
```

## Publishing

Before publishing, bump the version in `pyproject.toml`, run tests, and build
fresh distribution artifacts:

```cmd
python -m pip install --upgrade build twine
rmdir /s /q dist
python -m build
python -m twine check dist/*
```

Upload to TestPyPI first:

```cmd
python -m twine upload --repository testpypi dist/*
```

Install from TestPyPI in a clean environment and smoke-test the CLI. Then upload
the same checked artifacts to PyPI:

```cmd
python -m twine upload dist/*
```

## Status

BioEMMA is not yet a stable release. Before publishing to PyPI, the package
still needs a final check of bundled data, license compatibility, and user-facing
visualization dependencies.


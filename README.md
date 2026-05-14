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
layout and mapped identifiers.

With `output_dir`, BioEMMA writes:

```text
out/rn00010/
  map.json
  kegg_reconstruction.json
  summary.json
  fluxes.json              # when fluxes are provided or run_fba=True
  map.html                 # when save_html=True or save_png=True
  map.png                  # when save_png=True
  map_with_fluxes.html     # when flux data and HTML/PNG output are requested
  map_with_fluxes.png      # when flux data and PNG output are requested
```

HTML output requires the `escher` package. PNG output additionally requires
`playwright` and a browser installed for Playwright. These visualization
dependencies are not installed automatically by BioEMMA.

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
out/merged_map.json
```

Use `--no-merge` to skip the merged map.

The legacy single-file JSON output is still available:

```bash
bioemma build --model path/to/model.xml --kgml path/to/rn00010.xml --output map.json
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


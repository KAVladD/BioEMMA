# BioEMMA

BioEMMA is an early-stage Python library for building Escher-compatible
metabolic maps from KEGG pathway layouts and genome-scale metabolic models.

The current main workflow is:

1. Parse a KEGG KGML/XML pathway with `KeggMap`.
2. Convert KEGG compounds and reactions to BiGG/SEED identifiers using bundled
   MetaNetX-derived mapping tables.
3. Build an Escher JSON map with `EscherMapper`.

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

```python
import cobra

from bioemma import KeggMap, EscherMapper


kegg_map = KeggMap()

with open("path/to/rn00010.xml", "r", encoding="utf-8") as file:
    kegg_map.read_from_file(file)

metabolites = kegg_map.get_metabolites()
reactions = kegg_map.get_reactions()

model = cobra.io.read_sbml_model("path/to/model.xml")

mapper = EscherMapper(
    metabolites=metabolites,
    reactions=reactions,
    database="BIGG",
)

escher_map = mapper.build_map(model)
```

`escher_map` is a Python object compatible with the Escher JSON map structure.
It can be serialized with `json.dump`.

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

The script for regenerating mapping tables is kept separately in:

```text
scripts/prepare_db_mapping.py
```

## Status

BioEMMA is not yet a stable release. Before publishing to PyPI, the package
still needs packaging metadata, tests, and a final check of bundled data and
license compatibility.


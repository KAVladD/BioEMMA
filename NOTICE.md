# Third-Party Data Notice

BioEMMA includes compact runtime mapping tables:

- `src/bioemma/resources/metabolite_mapping.tsv`
- `src/bioemma/resources/reaction_mapping.tsv`

These files are derived from MetaNetX/MNXref cross-reference data. MetaNetX
states that, except where otherwise noted, data available from its site are
licensed under the Creative Commons Attribution 4.0 International License
(CC BY 4.0). MetaNetX also notes that MNXref uses information sourced from
external resources and that licensing agreements for those resources are
specified in the downloadable files.

Users are responsible for complying with the licenses and terms of the
underlying data sources referenced by the mapping tables, including MetaNetX,
KEGG, BiGG, SEED, and other external databases.

The BioEMMA command-line interface can download KGML files from the KEGG REST
API when `--pathway` is used. KEGG states that the KEGG API at `rest.kegg.jp`
is available only for academic use by academic users, and asks users to limit
API calls to no more than three requests per second.

References:

- MetaNetX/MNXref namespace: https://beta.metanetx.org/mnxdoc/mnxref.html
- KEGG API: https://www.genome.jp/kegg/rest/
- KEGG commercial licensing: https://kegg.net/en/licensing.html

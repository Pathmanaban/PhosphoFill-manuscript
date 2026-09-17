# Canonical proteome analysis table

`proteome_results_complete_canonical.tsv.gz` is the consolidated
PhosphoFill output table used for the proteome-scale analyses.

- 306,900 rows: 102,300 phosphosites × three ranked poses.
- Uncompressed size: 241,267,715 bytes.
- Compressed size: 23,678,898 bytes.
- SHA-256 of the uncompressed TSV:
  `c8b7681b52ebcd7503406a5a2a13dae17098a0ada518018dd8ea5a803234af99`.

The gzip file was verified by streaming decompression and comparing its SHA-256
with the source TSV. Pandas reads it directly; no manual decompression is
required:

```python
import pandas as pd
results = pd.read_csv(
    "proteome_results_complete_canonical.tsv.gz",
    sep="\t",
)
```

The table contains numeric and categorical outputs for all three independently
minimised ranks. Per-protein coordinate models, rich JSON/TSV reports, and HTML
visualisations are not duplicated in this repository; those storage-heavy
artifacts are distributed through Scop3P.


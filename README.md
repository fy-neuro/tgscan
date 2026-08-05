# tgscan

**Transgene hitchhiker gene screener** — identify and validate genes captured by BAC transgene constructs that cause RNA-seq count artifacts in transgenic mouse models.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)

## Overview

When a BAC transgene construct is built, it sometimes captures neighboring gene exons from the genome. After multi-copy integration and transcription, these "hitchhiker" sequences are miscounted as endogenous gene expression, causing RNA-seq count artifacts.

`tgscan` provides:
- **8 GEO matrix format parsers** (xlsx, csv.gz, tsv.gz, txt.gz, xls.gz, RAW.tar simple/Kallisto, h5ad)
- **Two-stage validation**: Pearson correlation + percentile, then cis-enrichment hypergeometric test
- **CLI** for single/batch verification + catalog search
- **Catalog** of 14 confirmed + 4 candidate hitchhikers (as of 2026-08)

## Installation

```bash
pip install -e .
# Optional: h5ad support
pip install -e ".[h5ad]"
# Optional: development
pip install -e ".[dev]"
```

## Quick start

### Command line

```bash
# Verify a single candidate
tgscan verify \
  --matrix GSE130842.xlsx \
  --driver Nfil3 \
  --candidate Auh \
  --gtf mm39.gtf

# Batch validate
tgscan batch --jobs candidates.tsv --geo-dir ./geo_data --output results.tsv

# Run cis-enrichment on stage1 results
tgscan cis --input stage1.tsv --gtf mm39.gtf

# Search known hitchhikers
tgscan catalog list
tgscan catalog search --driver Nfil3
```

### Python API

```python
from tgscan import verify

result = verify(
    matrix='GSE130842.xlsx',
    driver='Nfil3',
    candidate='Auh',
    gtf='mm39.gtf'
)
print(result.verdict)    # "CONFIRMED"
print(result.r)          # 0.965
print(result.cis_p)      # 6.8e-06
```

## Status definitions

| Status | Meaning |
|--------|---------|
| ✅ Confirmed | Percentile >95% + cis-enrichment p < 1e-3 |
| ⚠️ Candidate | Percentile >95% + cis-enrichment p in [1e-3, 1e-2] |
| Stage 1 passed | Percentile >95% but cis not tested |
| Negative | Verified, no signal in clean background |
| Insufficient | Data format/design/quality prevents assessment |
| Predicted only | Candidate predicted from BAC span but not yet validated |

See [Methods](https://github.com/fy-neuro/tgscan/blob/main/docs/methods.md) for cis-enrichment test details.

## Supported formats

| Format | Extension | Notes |
|--------|-----------|-------|
| Excel | `.xlsx` | Skips DE-result sheets (LRT/WALD) |
| Excel gzipped | `.xlsx.gz` | |
| Legacy Excel | `.xls.gz` | Via calamine (xlrd 2.0 dropped support) |
| CSV / TSV / TXT | `.csv`, `.tsv`, `.txt` (+`.gz`) | Auto-separator detection |
| 10x RAW.tar | `_RAW.tar` | Simple per-sample + Kallisto transcript-level |
| h5ad | `.h5ad` (+`.gz`) | scanpy pseudo-bulk by cell type (optional dep) |

## Catalog

The catalog of known hitchhikers is bundled with the package. See:
- [Online catalog](http://47.97.243.155/)
- `data/known_hitchhikers.tsv` in this repo

## License

MIT

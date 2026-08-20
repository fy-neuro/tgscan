# tgscan

**Transgene hitchhiker gene screener** — identify and validate genes captured by BAC transgene constructs that cause RNA-seq count artifacts in transgenic mouse models.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)

## Overview

When a BAC transgene construct is built, it sometimes captures neighboring gene exons from the genome. After multi-copy integration and transcription, these "hitchhiker" sequences are miscounted as endogenous gene expression, causing RNA-seq count artifacts.

`tgscan` provides:
- **9 GEO matrix format parsers** (xlsx, csv.gz, tsv.gz, txt.gz, xls.gz, RAW.tar simple/Kallisto/Cufflinks, h5ad pseudo-bulk)
- **Three-stage validation (v0.3)**: Stage 0 design gate (FACS-sorted / single-cell / TRAP-IP / z-score / miRNA-only / over-filtered matrices blocked, curated per-GSE design registry), Stage 1 Pearson + genome-wide percentile (with background-dirty guard), Stage 2 cis-1Mb hypergeometric enrichment (the confirmation gold standard, p<1e-3)
- **Design gate**: miRNA-only / over-filtered / z-score matrices flagged automatically; v0.3 adds FACS-fraction (GFP+/GFP-), single-cell (cells-as-samples) and TRAP/RiboTag detection + `data/known_designs.tsv` registry (ground truth: GSE83356/GSE115934/GSE127845)
- **Construct gate**: empirically excluded alleles (promoter cassettes that cannot capture neighbors, e.g. Vil1-cre 12.4kb) are auto-skipped in batch runs
- **Bundled catalog**: 21 confirmed + 6 candidate hitchhikers (2026-08-17) with evidence levels
- **CLI + Python API** for single/batch verification

## Installation

```bash
pip install .                # from a checkout
pip install -e ".[dev]"      # development install (pytest)
pip install ".[h5ad]"        # + scanpy for h5ad pseudo-bulk
```

> **Note (exFAT/NAS users)**: building directly from a source tree on exFAT can
> stall (setuptools writes a `build/` dir; slow + case-insensitive FS). Copy the
> tree to a local ext4/tmpfs dir first, or install from the git URL — pip then
> builds in an ext4 temp dir automatically.

## Quick start

### 0. Selftest (do this first)

```bash
tgscan selftest
```

**Evidence card (v0.3 module)** — transparent per-gene dossier instead of a one-word verdict:

```bash
tgscan card -d Lfng -c Ttyh3 -s evidence_store.tsv --ssot known_hitchhikers.tsv
# shows every dataset (r/n/pct/background/status), pooled r + CI + I²,
# Fisher-combined cis p, structural capture, and automatic flags
# (single-dataset / small-n / r-spread / sign-flip / capture-unverified / lineage-gap)
# labels: replicated | needs-review | L1-channel. Formats: --format text|tsv|json
```
Build the store with `Task/Task25_validation_queue/build_evidence_store.py`.

Runs the full pipeline (parser → stage 1 → cis) on a bundled synthetic
locus — a co-amplified hitchhiker and a negative control. No downloads needed.
All checks must print `PASS`.

### 1. Verify one candidate

```bash
tgscan verify --matrix GSE130842.xlsx --driver Nfil3 --candidate Auh --gtf mm39.gtf
```

Output: Stage-1 r / percentile / background mean, Stage-2 cis p (Top-10/50/100)
and a verdict: **CONFIRMED** (cis p<1e-3) / **CANDIDATE** (1e-3..1e-2) / WEAK.

### 2. Batch verify

```bash
tgscan batch --jobs jobs.tsv --gtf mm39.gtf --geo-dir /path/to/geo_files -o results.tsv
```

`jobs.tsv` columns: `gse, driver, candidate, dist_kb, allele` (optional `matrix`
for explicit paths). Resumable; alleles on the exclusion list are skipped with
status `excluded_construct`.

### 3. Catalog

```bash
tgscan catalog list --status confirmed
tgscan catalog search --driver Aldh1l1
tgscan catalog stats
```

### Python API

```python
from tgscan import verify
r = verify("GSE130842.xlsx", "Nfil3", "Auh", "mm39.gtf")
print(r.verdict, r.r, r.cis_p)   # CONFIRMED 0.965 6.8e-06
```

## How to read the numbers

- **Percentile is dataset-relative**: it is the rank of the driver–candidate
  Pearson r within *that* matrix's genome-wide r distribution. In a clean
  background (mean r ≈ 0) even a modest r can rank >99%. Confirmation therefore
  always requires the independent cis test — never percentile alone.
- **BACKGROUND_TOO_HIGH** (mean r > 0.3): the design (cell-type gradient,
  infection response, sort-marker dominance) swamps the transgene-load signal —
  the dataset is unusable, not the candidate negative.
- **Promoter-cassette trap**: a 12-kb promoter construct cannot capture
  neighbors; cell-type gradients can still produce driver–neighbor correlations
  that pass both stages. The construct gate + `excluded_constructs.tsv`
  guard against this (see Task25 audit, 2026-08-17).

## Data & provenance

Catalog SSOT: `src/tgscan/data/known_hitchhikers.tsv` (26 rows; Ensembl IDs
verified against mm39). Live browser: [tgscan catalog](http://47.97.243.155/).

## License

MIT

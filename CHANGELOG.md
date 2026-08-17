# Changelog

## 0.2.0 (2026-08-17)

Focus: field-test fixes from the Task24/25 screening campaigns + first-run
usability. 15/15 tests pass; `tgscan selftest` green on clean venv install.

### Fixed
- **Packaging**: catalog TSV now ships *inside* the package
  (`src/tgscan/data/`) and loads via `importlib.resources`. v0.1.0 pointed
  outside site-packages after install — every `tgscan catalog` command
  crashed. Also `.gitignore` added; stale `build/` dirs caused wheel-build
  failures on exFAT checkouts.
- **Cufflinks duplicate value columns** (GSE71463 style): repeated `FPKM`
  headers deduped by pandas to `FPKM.1/.2` are now merged (mean), previously
  produced wrong-column parses or "no common genes" errors.
- **Cufflinks tracking files inside RAW.tar**: header-aware column selection
  (`gene_short_name` + first `FPKM`) instead of positional cols.

### Added
- **`tgscan selftest`**: end-to-end installation check on bundled synthetic
  data (synthetic GTF + count matrix, no downloads). Includes regression
  coverage for numeric-leading MGI symbols (`0610005C13Rik`) — the root cause
  of the NO_VALID_R failures in the 2026-08-17 web pipeline.
- **Construct gate**: `data/excluded_constructs.tsv` (empirical blacklist from
  PI decisions + Task25 audits); `verify --batch` skips excluded alleles with
  status `excluded_construct`. Python: `catalog.construct_gate()` /
  `construct_warning()`. Note: `*` in MGI symbols is a literal mutant marker
  (`Tg(Prnp*)CDah`), never a wildcard.
- **Design guard**: z-score normalized matrices flagged
  (`zscore_normalized`) via negative-value fraction.
- Tests for all of the above (7 new; 15 total).

### Changed
- Catalog SSOT updated to 26 rows: 21 confirmed + 6 candidate
  (new: Ttyh3/Lfng, Usp47/Dkk3; Arsi stats completed via full-depth
  re-analysis).
- README rewritten: selftest-first quickstart, number-reading guide,
  exFAT build note.

## 0.1.0 (2026-08-05)

Initial MVP: 8 GEO parsers, two-stage validation (percentile + cis),
CLI (verify/batch/cis/catalog), 26-gene catalog, real-data validation
report (5/5).

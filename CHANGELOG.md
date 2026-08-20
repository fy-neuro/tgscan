# Changelog

## 0.4.0 (2026-08-20)

Focus: **Bayesian reporting channels** in the evidence card (designs A/C/D of
`BAYESIAN_DESIGN.md`), developed on branch `bayes-v2` off the sealed v0.3.1
baseline and merged after PI review. The calibrated scorer deliberately
**keeps v1 bare features**: a 3-variant LOO comparison on the same 141-pair
adjudication set gave v1 0.977 > v3 (posterior z) 0.912 > v2 (P_theta)
0.882 — the design gate has already cleared the "high-r-but-fake" cases, so
posterior transforms add no discriminative gain there. Bayesian channels and
the gate are complementary, not replacements.

### Added
- **`tgscan/bayes.py`**: DerSimonian-Laird τ² + normal-conjugate posterior
  (no MCMC needed at k=2-4). Card gains `tau2 / pooled_r_re / r_ci95_re /
  P_theta / bayes_band` (design A) and `sign_consistency_p` — a Beta(2,8)
  posterior turning the binary sign-flip flag into a continuous suspicion
  score (design D; Cerkl/Ttc12 = 0.500 vs Ankk1 = 0.969 in the pilot).
- **`bg_sd / z_abs` absolute background anchor** (design C, simplified):
  full empirical-Bayes shrinkage was quantified (<1e-3 per-dataset effect at
  ~20k background genes) and dropped; the honest simplification is noted in
  code comments.
- `BAYESIAN_DESIGN.md`: full design doc (math, adaptation, pilot table,
  decision points for PI) + `bayesian_pilot.py` companion for recomputation.

### Changed
- **card**: cis combination now pools only design-clean rows; registry gains
  GSE137572 (single-cell, BLOCKED) and GSE182815 (filter, BLOCKED).
- **design**: paired plus/minus word-fraction detection (Venus_plus/
  Venus_minus); single-arm genotype labels stay unblocked per the pairing rule.

### Data
- `known_hitchhikers.tsv`: 29 rows = 21 confirmed + 8 candidate — Cerkl
  retracted to candidate (sole dataset GSE137572 is a 10x per-cell matrix);
  Ttc12 added (L2c).

### Tests
- pytest 23 + selftest green; bayes unit tests use the pilot digits as gold
  standard (digit-identical regression).

## 0.3.0 (2026-08-20)

Focus: **design gate as Stage 0** + evidence transparency. Ground truth from
the Ttyh3/Cerkl design review: FACS-sorted / single-cell / TRAP-IP designs
manufacture spurious correlations (GSE83356: cross-group r=-0.77 with
within-GFP+ r=+0.91). 22 pytest + 14 selftest checks green; real-data
regression vs the SSOT digit-identical.

### Added
- **Stage 0 design gate** (`tgscan/design.py`): matrix-intrinsic detection of
  FACS-sorted fraction pairs (GFP+/GFP-), single-cell signatures (cells as
  samples: >24 columns with >80% zeros), TRAP/RiboTag/IP-enriched sample
  names — on top of the existing miRNA-only / over-filtered / z-score /
  n<3 checks. Wired into `verify()`/`verify_batch()` (previously imported
  but never called). `--gse` consults a curated per-GSE design registry
  (`data/known_designs.tsv`, seeded with GSE83356/115934/127845 + 2 more);
  `--no-gate` preserves legacy behaviour.
- **`tgscan design`** subcommand: audit a matrix + GSE against the gate
  without running correlations.
- **`tgscan card`**: evidence dossier per (driver, gene) — every dataset row,
  Fisher-z pooling, I2, sign-flip flag, cis combination, structure capture;
  design-unclean store rows are excluded from the usable pool.
- **`tgscan score`**: transparent logistic probability (all coefficients
  disclosed, LOO-CV reported) — ranking aid, not doctrine.
- `VerifyResult.design_issues` / `design_verdict`; batch output gains a
  `design_issues` column.

### Fixed
- **Composite-ID header parsing** (`_parse_simple`): files whose first column
  is comma-composite (`id_gene,gene_name,gene_type<TAB>samples...`) made the
  separator sniffer pick "," and collapse all sample columns into one field —
  `tgscan verify` failed with "too few samples (0)" (GSE94145). Headers with
  tabs now parse as TSV and the composite column expands into real columns.
- `screen_design` registry-BLOCKED verdict could be silently downgraded to
  PASS (missing branch); caught with a clean 15k-gene control matrix.

### Data
- `known_hitchhikers.tsv`: 29 rows = 22 confirmed + 7 candidate (Ankk1
  upgraded on two-dataset evidence; Ttc12 candidate with disclosed sign-flip;
  Drd2 BAC RP23-161H15 capture verified from the original publication).

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

"""Example 1: Verify a single candidate hitchhiker.

This is the simplest way to use tgscan. Given a GEO supplementary file,
a driver gene, and a candidate gene, tgscan will:
  1. Parse the matrix (auto-detect format)
  2. Compute Pearson r + percentile (Stage 1)
  3. Run cis-enrichment test (Stage 2) if Stage 1 passes
  4. Return a verdict (CONFIRMED / CANDIDATE / WEAK_OR_FALSE_POSITIVE)
"""
from tgscan import verify

# Replace these paths with your own
result = verify(
    matrix_path='GSE130842_Count_table_Delacher_et_al_2019.xlsx',  # GEO supp file
    driver='Nfil3',           # driver gene symbol
    candidate='Auh',          # candidate hitchhiker
    gtf_path='mm39.gtf',      # GTF for mouse genome
)

print(f"Driver:    {result.driver}")
print(f"Candidate: {result.gene}")
print(f"Matrix:    {result.matrix_format}")
if result.error:
    print(f"ERROR:     {result.error}")
else:
    print(f"\nStage 1:")
    print(f"  r = {result.r:.4f}")
    print(f"  Percentile = {result.percentile:.2f}%")
    print(f"  Status: {result.analysis.status}")
    if result.cis:
        print(f"\nStage 2 (cis-enrichment):")
        print(f"  Top-10 p: {result.cis.cis_top10_p:.2e}  (fold {result.cis.fold_top10:.1f}x)")
        print(f"  Top cis hits: {', '.join(result.cis.top_cis_genes)}")
    print(f"\n>>> VERDICT: {result.verdict}")

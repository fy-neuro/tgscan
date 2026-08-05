"""Example 2: Batch validate candidates from a TSV.

The TSV should have columns: gse, driver, candidate, dist_kb, allele
Optionally: matrix (full path to file), if you want to skip auto-detection.

tgscan will:
  - Parse each matrix once (shared across candidates of the same GSE)
  - Compute Stage 1 Pearson + percentile for each candidate
  - Save incrementally to output TSV (resumable)
"""
from tgscan import verify_batch

output = verify_batch(
    jobs_tsv='candidates.tsv',         # see format above
    gtf_path='mm39.gtf',
    output_tsv='results.tsv',
    geo_dir='/path/to/geo_files',      # directory with GEO supp files
    resume=True,                        # skip GSEs already in output
)

print(f"\nResults saved to: {output}")

# Read and summarize
import pandas as pd
df = pd.read_csv(output, sep='\t')
print(df['status'].value_counts())

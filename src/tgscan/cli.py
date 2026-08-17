"""tgscan CLI — typer-based command-line interface."""
from __future__ import annotations
from pathlib import Path
from typing import Optional
import typer
import json
from . import __version__
from .runner import verify, verify_batch
from .parsers import parse_matrix, detect_format
from .gtf import get_gtf_index
from .cis import cis_test
from .parsers import numeric_sample_cols
from . import catalog as cat_mod

app = typer.Typer(
    name='tgscan',
    help='Transgene hitchhiker gene screener',
    no_args_is_help=True,
    add_completion=False,
)


@app.command(name='verify')
def verify_cmd(
    matrix: Path = typer.Option(..., '--matrix', '-m', help='Path to GEO supplementary file'),
    driver: str = typer.Option(..., '--driver', '-d', help='Driver gene symbol'),
    candidate: str = typer.Option(..., '--candidate', '-c', help='Candidate hitchhiker gene symbol'),
    gtf: Path = typer.Option(..., '--gtf', '-g', help='Path to GTF file'),
    no_cis: bool = typer.Option(False, '--no-cis', help='Skip cis-enrichment test (stage 1 only)'),
    json_out: bool = typer.Option(False, '--json', help='Output as JSON'),
):
    """Verify a single candidate hitchhiker."""
    result = verify(
        matrix_path=str(matrix),
        driver=driver,
        candidate=candidate,
        gtf_path=str(gtf),
        run_cis=not no_cis,
    )
    if json_out:
        typer.echo(json.dumps(result.to_dict(), indent=2, default=str))
    else:
        typer.echo(f"\n=== Verification Result ===")
        typer.echo(f"Driver:    {result.driver}")
        typer.echo(f"Candidate: {result.gene}")
        typer.echo(f"Matrix:    {result.matrix_format}")
        if result.error:
            typer.echo(f"ERROR:     {result.error}")
        else:
            typer.echo(f"\n--- Stage 1 (Pearson + percentile) ---")
            if result.analysis:
                a = result.analysis
                typer.echo(f"  r = {a.r:.4f}  (p = {a.p_value:.2e})")
                typer.echo(f"  Percentile: {a.percentile:.2f}%  (genome-wide mean_r = {a.mean_r:.4f})")
                typer.echo(f"  Status: {a.status}")
            if result.cis:
                typer.echo(f"\n--- Stage 2 (cis-enrichment) ---")
                c = result.cis
                typer.echo(f"  Top-10 p:  {c.cis_top10_p:.2e}  (fold {c.fold_top10:.1f}x)")
                typer.echo(f"  Top-50 p:  {c.cis_top50_p:.2e}  (fold {c.fold_top50:.1f}x)")
                typer.echo(f"  Top-100 p: {c.cis_top100_p:.2e}  (fold {c.fold_top100:.1f}x)")
                typer.echo(f"  Top cis genes: {', '.join(c.top_cis_genes)}")
            typer.echo(f"\n>>> VERDICT: {result.verdict}")


@app.command()
def batch(
    jobs: Path = typer.Option(..., '--jobs', '-j', help='TSV with columns: gse, driver, candidate, dist_kb, allele'),
    gtf: Path = typer.Option(..., '--gtf', '-g', help='Path to GTF file'),
    output: Path = typer.Option('results.tsv', '--output', '-o', help='Output TSV path'),
    geo_dir: Optional[Path] = typer.Option(None, '--geo-dir', help='Directory containing GEO files'),
    no_resume: bool = typer.Option(False, '--no-resume', help='Disable resume (re-run all)'),
):
    """Batch verify candidates from a TSV."""
    verify_batch(
        jobs_tsv=str(jobs),
        gtf_path=str(gtf),
        output_tsv=str(output),
        geo_dir=str(geo_dir) if geo_dir else None,
        resume=not no_resume,
    )


@app.command(name='cis')
def cis_cmd(
    matrix: Path = typer.Option(..., '--matrix', '-m', help='Path to GEO supplementary file'),
    driver: str = typer.Option(..., '--driver', '-d', help='Driver gene symbol'),
    gtf: Path = typer.Option(..., '--gtf', '-g', help='Path to GTF file'),
):
    """Run cis-enrichment test on a single driver (find all cis hits)."""
    gtf_idx = get_gtf_index(str(gtf))
    df, fmt = parse_matrix(str(matrix))
    id_col = df.columns[0]
    sample_cols = numeric_sample_cols(df, id_col)
    drv_eid = gtf_idx.symbol_to_ensembl(driver)
    if drv_eid is None:
        typer.echo(f"Driver '{driver}' not in GTF", err=True)
        raise typer.Exit(1)
    result = cis_test(df, drv_eid, sample_cols, gtf_idx, id_col=id_col)
    typer.echo(f"\n=== cis-enrichment for {driver} ===")
    typer.echo(f"Top-10 p:  {result.cis_top10_p:.2e}  (fold {result.fold_top10:.1f}x)")
    typer.echo(f"Top-50 p:  {result.cis_top50_p:.2e}  (fold {result.fold_top50:.1f}x)")
    typer.echo(f"Top-100 p: {result.cis_top100_p:.2e}  (fold {result.fold_top100:.1f}x)")
    typer.echo(f"Top cis genes: {', '.join(result.top_cis_genes)}")
    typer.echo(f"\n>>> VERDICT: {result.verdict}")


catalog_app = typer.Typer(name='catalog', help='Known hitchhiker catalog operations')
app.add_typer(catalog_app, name='catalog')


@catalog_app.command(name='list')
def catalog_list(
    status: Optional[str] = typer.Option(None, '--status', '-s',
                                          help='Filter by status: confirmed / candidate'),
):
    """List known hitchhikers."""
    if status == 'confirmed':
        df = cat_mod.list_confirmed()
    elif status == 'candidate':
        df = cat_mod.list_candidates()
    else:
        df = cat_mod.list_all()
    if len(df) == 0:
        typer.echo("(empty)")
        return
    typer.echo(df[['gene', 'driver', 'gene_biotype', 'status', 'r', 'percentile',
                    'cis_top10_p']].to_string(index=False))


@catalog_app.command(name='search')
def catalog_search(
    driver: Optional[str] = typer.Option(None, '--driver', '-d'),
    gene: Optional[str] = typer.Option(None, '--gene', '-g'),
):
    """Search catalog."""
    df = cat_mod.search(driver=driver, gene=gene)
    if len(df) == 0:
        typer.echo("(no matches)")
        return
    typer.echo(df[['gene', 'driver', 'status', 'r', 'percentile']].to_string(index=False))


@catalog_app.command(name='stats')
def catalog_stats():
    """Show catalog statistics."""
    s = cat_mod.stats()
    typer.echo(f"Total entries:        {s['total']}")
    typer.echo(f"Confirmed:            {s['confirmed']}")
    typer.echo(f"Candidate:            {s['candidate']}")
    typer.echo(f"Unique driver genes:  {s['unique_drivers']}")
    typer.echo(f"Unique hitchhikers:   {s['unique_hitchhikers']}")


@app.command()
def version():
    """Print version."""
    typer.echo(__version__)


@app.command(name='selftest')
def selftest_cmd():
    """End-to-end installation check with bundled synthetic data (no downloads)."""
    from .selftest import run_selftest
    typer.echo("Running tgscan selftest (synthetic GTF + matrix, ~5s)...")
    ok = run_selftest(verbose=True)
    if not ok:
        raise typer.Exit(1)


if __name__ == '__main__':
    app()

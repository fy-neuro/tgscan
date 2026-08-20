"""tgscan CLI — typer-based command-line interface."""
from __future__ import annotations
from pathlib import Path
from typing import Optional
import typer
import json
import csv
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


@app.command(name='score')
def score_cmd(
    driver: str = typer.Option(..., '--driver', '-d', help='Driver gene symbol'),
    candidate: str = typer.Option(..., '--candidate', '-c', help='Candidate gene symbol'),
    store: Path = typer.Option(..., '--store', '-s', help='evidence_store.tsv'),
    ssot: Path = typer.Option(None, '--ssot', help='SSOT TSV (training positives + structure)'),
):
    """Provisional v0.3 probability score (transparent logistic, calibrated on verdict set).
    Weights are disclosed; treat as ranking aid, not doctrine."""
    from . import score as score_mod, card as card_mod
    rows = card_mod.load_store(store)
    ssot_rows = list(csv.DictReader(open(ssot), delimiter='	')) if ssot else []
    X, y, _ = score_mod.build_training(rows, ssot_rows)
    m = score_mod.fit_logistic(X, y)
    ssot_row = next((r for r in ssot_rows
                     if r.get('driver') == driver and r.get('gene') == candidate), None)
    c, p, contrib = score_mod.score_pair(driver, candidate, rows, ssot_row, m)
    typer.echo(card_mod.render_text(c))
    typer.echo('')
    if p is None:
        typer.echo('score: n/a (no usable evidence rows — L1 channel?)')
        return
    typer.echo(f'P(hitchhiker) [provisional v0.3] = {p:.3f}')
    typer.echo('  contributions: ' + ', '.join(f'{k} {v:+.2f}' for k, v in
                                               sorted(contrib.items(), key=lambda t: -abs(t[1]))))


@app.command(name='card')
def card_cmd(
    driver: str = typer.Option(..., '--driver', '-d', help='Driver gene symbol'),
    candidate: str = typer.Option(..., '--candidate', '-c', help='Candidate hitchhiker gene symbol'),
    store: Path = typer.Option(..., '--store', '-s', help='evidence_store.tsv (per-dataset evidence rows)'),
    ssot: Optional[Path] = typer.Option(None, '--ssot', help='Optional SSOT TSV (adds status/structure)'),
    fmt: str = typer.Option('text', '--format', '-f', help='text | tsv | json'),
):
    """Generate an evidence card: all datasets, pooled stats, flags (no single-word verdicts)."""
    from . import card as card_mod
    rows = card_mod.load_store(store)
    ssot_row = None
    if ssot and Path(ssot).exists():
        with open(ssot) as fh:
            for r in csv.DictReader(fh, delimiter='\t'):
                if r.get('driver') == driver and r.get('gene') == candidate:
                    ssot_row = r
                    break
    c = card_mod.make_card(driver, candidate, rows, ssot_row)
    if fmt == 'json':
        typer.echo(card_mod.render_json(c))
    elif fmt == 'tsv':
        typer.echo('\t'.join(card_mod.CARD_TSV_FIELDS))
        typer.echo(card_mod.render_tsv_row(c))
    else:
        typer.echo(card_mod.render_text(c))


@app.command(name='verify')
def verify_cmd(
    matrix: Path = typer.Option(..., '--matrix', '-m', help='Path to GEO supplementary file'),
    driver: str = typer.Option(..., '--driver', '-d', help='Driver gene symbol'),
    candidate: str = typer.Option(..., '--candidate', '-c', help='Candidate hitchhiker gene symbol'),
    gtf: Path = typer.Option(..., '--gtf', '-g', help='Path to GTF file'),
    gse: Optional[str] = typer.Option(None, '--gse', help='GSE accession (enables curated design-registry lookup)'),
    no_cis: bool = typer.Option(False, '--no-cis', help='Skip cis-enrichment test (stage 1 only)'),
    no_gate: bool = typer.Option(False, '--no-gate', help='Disable v0.3 design gate (legacy behaviour)'),
    json_out: bool = typer.Option(False, '--json', help='Output as JSON'),
):
    """Verify a single candidate hitchhiker (v0.3: design gate first)."""
    result = verify(
        matrix_path=str(matrix),
        driver=driver,
        candidate=candidate,
        gtf_path=str(gtf),
        run_cis=not no_cis,
        gse=gse,
        gate=not no_gate,
    )
    if json_out:
        typer.echo(json.dumps(result.to_dict(), indent=2, default=str))
        return
    typer.echo(f"\n=== Verification Result ===")
    typer.echo(f"Driver:    {result.driver}")
    typer.echo(f"Candidate: {result.gene}")
    typer.echo(f"Matrix:    {result.matrix_format}")
    if result.design_verdict:
        tag = {'BLOCK': 'BLOCKED', 'WEAK': 'WEAK', 'PASS': 'PASS'}[result.design_verdict]
        typer.echo(f"Design:    {tag}" +
                   (f" ({'; '.join(result.design_issues)})" if result.design_issues else ""))
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


@app.command(name='design')
def design_cmd(
    matrix: Path = typer.Option(..., '--matrix', '-m', help='Path to GEO supplementary file'),
    gse: Optional[str] = typer.Option(None, '--gse', help='GSE accession (curated registry lookup)'),
):
    """Audit a matrix/GSE against the v0.3 design gate (no correlation run)."""
    from .design import screen_design, lookup_registry
    df, fmt = parse_matrix(str(matrix))
    id_col = df.columns[0]
    sample_cols = numeric_sample_cols(df, id_col)
    typer.echo(f"format: {fmt}  genes: {len(df)}  samples: {len(sample_cols)}")
    if sample_cols:
        typer.echo("sample cols: " + ", ".join(str(c) for c in sample_cols[:12]) +
                   (' ...' if len(sample_cols) > 12 else ''))
    res = screen_design(df, sample_cols, gse)
    reg = res['registry']
    if reg:
        typer.echo(f"registry [{res['source']}]: {reg.get('design')} -> {reg.get('verdict')}"
                   f"  ({reg.get('evidence')})")
    typer.echo(f"issues: {res['issues'] or '(none)'}")
    typer.echo(f">>> DESIGN GATE: {res['verdict']}")


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

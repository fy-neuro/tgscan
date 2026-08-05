"""High-level API: verify single candidate + batch verification."""
from __future__ import annotations
import os, signal, time
from typing import Optional, List, Dict, Any
import pandas as pd
from .models import VerifyResult, AnalysisResult, CisResult
from .parsers import parse_matrix, numeric_sample_cols
from .analysis import analyze
from .cis import cis_test
from .gtf import GtfIndex, get_gtf_index
from .design import detect_design_issues


class TimeoutError(Exception):
    pass


def _timeout_handler(signum, frame):
    raise TimeoutError(f"Operation timed out")


def verify(matrix_path: str, driver: str, candidate: str, gtf_path: str,
           driver_eid: Optional[str] = None, candidate_eid: Optional[str] = None,
           driver_entrez: Optional[str] = None, candidate_entrez: Optional[str] = None,
           run_cis: bool = True, timeout_sec: int = 120) -> VerifyResult:
    """Verify a single (driver, candidate) pair.

    Args:
        matrix_path: path to GEO supplementary file
        driver: driver gene symbol (e.g. 'Nfil3')
        candidate: candidate hitchhiker gene symbol (e.g. 'Auh')
        gtf_path: path to GTF file (mm39.gtf or similar)
        driver_eid/candidate_eid: optional Ensembl IDs (auto-fetched if not provided)
        driver_entrez/candidate_entrez: optional Entrez IDs (for Entrez-format matrices)
        run_cis: whether to run Stage 2 cis-enrichment test (default True)
        timeout_sec: timeout in seconds for the full pipeline

    Returns:
        VerifyResult with analysis (stage 1) and cis (stage 2, if run).

    Example:
        >>> from tgscan import verify
        >>> r = verify('GSE130842.xlsx', 'Nfil3', 'Auh', 'mm39.gtf')
        >>> print(r.verdict, r.r, r.cis_p)
        CONFIRMED 0.965 6.8e-06
    """
    # Set timeout
    old_handler = signal.signal(signal.SIGALRM, _timeout_handler)
    signal.alarm(timeout_sec)
    try:
        return _verify_inner(matrix_path, driver, candidate, gtf_path,
                             driver_eid, candidate_eid, driver_entrez, candidate_entrez,
                             run_cis)
    except TimeoutError as e:
        result = VerifyResult(gene=candidate, driver=driver, matrix_format='?',
                              error=f'timeout ({timeout_sec}s)')
        return result
    except Exception as e:
        result = VerifyResult(gene=candidate, driver=driver, matrix_format='?',
                              error=f'{type(e).__name__}: {e}')
        return result
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old_handler)


def _verify_inner(matrix_path, driver, candidate, gtf_path,
                  driver_eid, candidate_eid, driver_entrez, candidate_entrez, run_cis):
    gtf = get_gtf_index(gtf_path)
    # resolve Ensembl IDs
    if driver_eid is None:
        driver_eid = gtf.symbol_to_ensembl(driver)
    if candidate_eid is None:
        candidate_eid = gtf.symbol_to_ensembl(candidate)

    # parse matrix
    df, fmt = parse_matrix(matrix_path)

    id_col = df.columns[0]
    sample_cols = numeric_sample_cols(df, id_col)
    if len(sample_cols) < 3:
        raise ValueError(f"too few samples ({len(sample_cols)})")

    # find driver
    drv_row, _ = gtf.find_in_matrix(df, driver, driver_eid, driver_entrez)
    if drv_row is None:
        raise ValueError(f"driver '{driver}' not in matrix")

    # find candidate
    cand_row, _ = gtf.find_in_matrix(df, candidate, candidate_eid, candidate_entrez)
    if cand_row is None:
        raise ValueError(f"candidate '{candidate}' not in matrix")

    # stage 1
    analysis = analyze(df, drv_row, cand_row, sample_cols)

    # stage 2 (cis)
    cis_res = None
    if run_cis and driver_eid is not None and analysis.status in ('HIGH_CONFIDENCE', 'MODERATE'):
        try:
            cis_res = cis_test(df, driver_eid, sample_cols, gtf, id_col=id_col)
        except Exception as e:
            # cis can fail (e.g., driver not in GTF). Keep analysis result.
            pass

    return VerifyResult(gene=candidate, driver=driver, matrix_format=fmt,
                        analysis=analysis, cis=cis_res)


def verify_batch(jobs_tsv: str, gtf_path: str, output_tsv: str,
                 geo_dir: Optional[str] = None, resume: bool = True,
                 timeout_sec: int = 120) -> str:
    """Batch verify candidates from a TSV file.

    Expected TSV columns: gse, driver, candidate, matrix (optional), dist_kb (optional), allele (optional)
    If matrix column is missing, requires geo_dir + gse to locate files.

    Args:
        jobs_tsv: input TSV
        gtf_path: GTF file
        output_tsv: output TSV
        geo_dir: directory containing GEO batch files (required if matrix not in TSV)
        resume: skip GSE already in output_tsv (incremental save)
        timeout_sec: per-job timeout

    Returns:
        Path to output_tsv.
    """
    gtf = get_gtf_index(gtf_path)
    jobs = pd.read_csv(jobs_tsv, sep='\t')

    # resume support
    existing_gse = set()
    if resume and os.path.exists(output_tsv):
        try:
            ex = pd.read_csv(output_tsv, sep='\t')
            existing_gse = set(ex.get('gse', []).unique())
            print(f"Resume: skipping {len(existing_gse)} GSEs already in {output_tsv}")
        except Exception:
            pass

    # group jobs by GSE for shared parsing
    if 'gse' not in jobs.columns:
        raise ValueError("jobs_tsv must have 'gse' column")
    by_gse: Dict[str, List[dict]] = {}
    for _, row in jobs.iterrows():
        by_gse.setdefault(str(row['gse']), []).append(row.to_dict())

    results = []
    n_total = len(by_gse)
    for i, (gse, gse_jobs) in enumerate(by_gse.items(), 1):
        if gse in existing_gse:
            continue
        print(f"[{i}/{n_total}] {gse} ({len(gse_jobs)} candidates)", flush=True)
        # locate matrix file
        if 'matrix' in gse_jobs[0] and pd.notna(gse_jobs[0].get('matrix')):
            matrix_path = gse_jobs[0]['matrix']
        elif geo_dir is not None:
            matrix_path = _find_geo_file(geo_dir, gse)
            if matrix_path is None:
                for job in gse_jobs:
                    results.append(_error_result(job, gse, f'no file for {gse} in {geo_dir}'))
                continue
        else:
            raise ValueError("Either 'matrix' column or geo_dir must be provided")

        try:
            df, fmt = parse_matrix(matrix_path)
        except Exception as e:
            for job in gse_jobs:
                results.append(_error_result(job, gse, f'parse: {type(e).__name__}: {e}'))
            continue

        id_col = df.columns[0]
        sample_cols = numeric_sample_cols(df, id_col)

        # group by driver
        by_driver: Dict[str, List[dict]] = {}
        for job in gse_jobs:
            by_driver.setdefault(job['driver'], []).append(job)

        for driver, driver_jobs in by_driver.items():
            drv_eid = gtf.symbol_to_ensembl(driver)
            drv_row, _ = gtf.find_in_matrix(df, driver, drv_eid)
            if drv_row is None:
                for job in driver_jobs:
                    results.append(_error_result(job, gse, f"driver '{driver}' not in matrix"))
                continue
            # stage 1 for each candidate
            for job in driver_jobs:
                cand = job['candidate']
                cand_eid = gtf.symbol_to_ensembl(cand)
                cand_row, _ = gtf.find_in_matrix(df, cand, cand_eid)
                if cand_row is None:
                    results.append(_result_dict(job, gse, fmt, len(sample_cols),
                                                status=f'candidate_not_found'))
                    continue
                try:
                    an = analyze(df, drv_row, cand_row, sample_cols)
                except Exception as e:
                    results.append(_error_result(job, gse, f'analyze: {e}'))
                    continue
                results.append(_result_dict(job, gse, fmt, len(sample_cols),
                                            analysis=an))

        # incremental save every 5 GSE
        if i % 5 == 0:
            _save_results(results, output_tsv)
            print(f"  [checkpoint] saved {len(results)} results")

    _save_results(results, output_tsv)
    print(f"\nDONE. {len(results)} results -> {output_tsv}")
    return output_tsv


def _find_geo_file(geo_dir: str, gse: str) -> Optional[str]:
    """Find the supplementary data file for a GSE in geo_dir."""
    if not os.path.isdir(geo_dir):
        return None
    candidates = []
    for f in os.listdir(geo_dir):
        if f.startswith(gse + '_') and not f.endswith('supp.txt'):
            candidates.append(os.path.join(geo_dir, f))
    if not candidates:
        return None
    # prefer non-RAW.tar; if all .tar, take the first
    non_tar = [c for c in candidates if not c.endswith('.tar')]
    return (non_tar or candidates)[0]


def _result_dict(job, gse, fmt, n_samples, analysis: Optional[AnalysisResult] = None,
                  status: Optional[str] = None) -> dict:
    d = {
        'gse': gse,
        'allele': job.get('allele', ''),
        'driver': job['driver'],
        'candidate': job['candidate'],
        'dist_kb': job.get('dist_kb', ''),
        'format': fmt,
        'n_samples': n_samples,
        'status': status or (analysis.status if analysis else 'unknown'),
    }
    if analysis:
        d.update({
            'r': analysis.r, 'p': analysis.p_value,
            'pct': analysis.percentile, 'mean_r': analysis.mean_r,
            'n_genes': analysis.n_genes,
        })
    return d


def _error_result(job, gse, error: str) -> dict:
    return {
        'gse': gse, 'allele': job.get('allele', ''),
        'driver': job['driver'], 'candidate': job['candidate'],
        'dist_kb': job.get('dist_kb', ''), 'format': '',
        'n_samples': 0, 'status': 'error', 'error': error,
    }


def _save_results(results: list, output_tsv: str):
    if not results:
        return
    pd.DataFrame(results).to_csv(output_tsv, sep='\t', index=False)

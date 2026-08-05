"""Matrix parsers for 8 GEO supplementary file formats."""
from __future__ import annotations
import os, gzip, shutil, re, tarfile
from typing import Tuple, Optional
import pandas as pd


# Column names that should NOT be treated as samples
META_COLS = {
    'Chr','Start','End','Strand','Length','chr','start','end','strand','length',
    'gene_biotype','gene_type','Chromosome','Feature','Source','Type','description',
    'UniProt','Pfam','GO','KEGG','NR','NT','COG','Biotype','Position',
    'entrez_id','EntrezID','RefSeq','refseq','transcript_id','gene_id','Geneid',
    'gene_name','Gene','Symbol','name','Accession','Chromosome','accession','annotation',
    'gene_chr','gene_start','gene_end','gene_strand','gene_length',
    'gene_description','tf_family',
    'GC','gc','width','transcriptCount','BaseMean','baseMean','log2FC','lfcSE','pvalue','padj',
    'LocusTag','locus_tag',
}


def _is_meta_col(c, id_col: str) -> bool:
    """True if column should be excluded from samples."""
    if c == id_col or c in META_COLS:
        return True
    if isinstance(c, str):
        if c.startswith('gene_') and c not in ('gene_id', 'gene_name'):
            return True
        if c.startswith('Unnamed'):
            return True
        # edgeR / DESeq2 result columns
        cl = c.lower()
        for pat in ('logfc','logcpm','pvalue','fdr','padj','lr ','lrt ','wald ',
                    'avg.','avgp.','ave.','stat ','lfcse',
                    'log2ratio','log2foldchange','q value','q-value','qvalue',
                    'log2(','lfcase','ctrl_vs_','_vs_.log','comparison fpkm','control fpkm'):
            if pat in cl:
                return True
        for tok in ('.logFC', '.PValue', '.FDR', '.pvalue', '.fdr', '.padj',
                    '.log2FoldChange', '.baseMean', '.lfcSE', '.stat'):
            if c.endswith(tok) or c.lower().endswith(tok.lower()):
                return True
    return False


def numeric_sample_cols(df: pd.DataFrame, id_col: str):
    """Return list of column names that look like sample columns."""
    out = []
    for c in df.columns:
        if _is_meta_col(c, id_col):
            continue
        try:
            df[c].astype(float)
            out.append(c)
        except (ValueError, TypeError):
            continue
    return out


# ---------- format detection ----------
def detect_format(path: str) -> str:
    """Detect format from file extension."""
    p = path.lower()
    if p.endswith('.tar'):
        return 'raw_tar'
    if p.endswith('.xlsx'):
        return 'xlsx'
    if p.endswith('.xlsx.gz'):
        return 'xlsx_gz'
    if p.endswith('.xls.gz'):
        return 'xls_gz'
    if p.endswith('.h5ad.gz') or p.endswith('.h5ad'):
        return 'h5ad'
    if p.endswith('.csv.gz') or p.endswith('.csv'):
        return 'csv'
    if p.endswith('.tsv.gz') or p.endswith('.tsv'):
        return 'tsv'
    if p.endswith('.txt.gz') or p.endswith('.txt'):
        return 'txt'
    return 'unknown'


# ---------- parsers ----------
def _parse_simple(path: str) -> pd.DataFrame:
    """csv/tsv/txt.gz with auto-detect separator + composite ID normalization."""
    df = pd.read_csv(path, sep=None, engine='python', comment='#', on_bad_lines='skip')
    if len(df.columns) > 0:
        col = df.columns[0]
        s = df[col].astype(str)
        # pacbio composite: ENSMUST_X_ENSMUSG_Y → extract gene ID
        if s.head(50).str.contains('ENSMUST.*ENSMUSG').any():
            s = s.str.extract(r'(ENSMUSG\d+)', expand=False).fillna(s)
        # strip version suffixes
        s = s.str.replace(r'(ENSMUSG\d+)\.\d+', r'\1', regex=True)
        s = s.str.replace(r'(ENSMUSG\d+)-\d+', r'\1', regex=True)
        df[col] = s
    return df


def _parse_xlsx(path: str) -> pd.DataFrame:
    """Single/multi-sheet xlsx. Skip DE-result sheets."""
    xl = pd.ExcelFile(path, engine='openpyxl')
    de_patterns = ('LRT', 'WALD', 'DESeq', 'edgeR', 'log2FC', 'INSTRUCTIONS', 'README')
    data_sheets = [s for s in xl.sheet_names if not any(p in str(s) for p in de_patterns)]
    if not data_sheets:
        raise ValueError(f"DE-only xlsx (sheets: {xl.sheet_names})")
    # pick sheet with most columns (likely count matrix)
    best, best_ncol = None, 0
    for s in data_sheets:
        df = pd.read_excel(path, sheet_name=s, engine='openpyxl', nrows=5)
        if len(df.columns) > best_ncol:
            best, best_ncol = s, len(df.columns)
    return pd.read_excel(path, sheet_name=best, engine='openpyxl')


def _parse_xlsx_gz(path: str) -> pd.DataFrame:
    raw = path.replace('.gz', '')
    if not os.path.exists(raw):
        with gzip.open(path, 'rb') as fi, open(raw, 'wb') as fo:
            shutil.copyfileobj(fi, fo)
    return _parse_xlsx(raw)


def _parse_xls_gz(path: str) -> pd.DataFrame:
    """Legacy .xls.gz — try TSV first (.xls extension is often just text), then calamine."""
    raw = path.replace('.gz', '')
    if not os.path.exists(raw):
        with gzip.open(path, 'rb') as fi, open(raw, 'wb') as fo:
            shutil.copyfileobj(fi, fo)
    try:
        return pd.read_csv(raw, sep='\t', comment='#', on_bad_lines='skip')
    except Exception:
        return pd.read_excel(raw, engine='calamine')


def _parse_h5ad(path: str) -> pd.DataFrame:
    """scanpy h5ad — pseudo-bulk by cell type, skip z-score normalized data."""
    try:
        import scanpy as sc
        import numpy as np
    except ImportError as e:
        raise ImportError("h5ad parsing requires scanpy: pip install 'tgscan[h5ad]'") from e
    if path.endswith('.gz'):
        raw = path.replace('.gz', '')
        if not os.path.exists(raw):
            with gzip.open(path, 'rb') as fi, open(raw, 'wb') as fo:
                shutil.copyfileobj(fi, fo)
        path = raw
    adata = sc.read_h5ad(path)
    X = adata.X
    if hasattr(X, 'toarray'):
        X = X.toarray()
    else:
        X = np.asarray(X)
    # detect z-score normalization (per-row std≈1, sum≈0)
    sample_std = np.array([X[i].std() for i in range(min(100, X.shape[0]))])
    if 0.95 < np.median(sample_std) < 1.05:
        raise ValueError("z-score normalized scRNA-seq (per-row std≈1); pseudo-bulk mean-z unreliable, skipping")
    obs = adata.obs
    for cand_col in ['level2', 'level1', 'cell_type', 'cluster', 'louvain', 'leiden', 'batch', 'timepoint']:
        if cand_col in obs.columns:
            group_col = cand_col
            break
    else:
        group_col = obs.columns[0]
    groups = obs[group_col].astype(str).unique().tolist()
    pb = {}
    for g in groups:
        mask = (obs[group_col].astype(str) == g).values
        if mask.sum() < 30:
            continue
        pb[g] = X[mask].sum(axis=0)
    df = pd.DataFrame(pb, index=adata.var_names).reset_index()
    df.columns = ['Symbol'] + list(df.columns[1:])
    if 'Accession' in adata.var.columns:
        df['_Ensembl'] = adata.var['Accession'].values
    return df


def _parse_raw_tar(path: str, extract_dir: Optional[str] = None) -> pd.DataFrame:
    """RAW.tar — auto-detects simple (gene, count) and Kallisto transcript-level formats.

    Raises ValueError for 10x cellranger 3-file format (barcodes/features/matrix.mtx).
    """
    if extract_dir is None:
        extract_dir = path + '_extracted'
    os.makedirs(extract_dir, exist_ok=True)
    parseable_exts = ('.txt.gz', '.txt', '.csv.gz', '.csv', '.tsv.gz', '.tsv')
    tarfiles = [f for f in os.listdir(extract_dir) if f.endswith(parseable_exts)]
    if not tarfiles:
        with tarfile.open(path) as tar:
            tar.extractall(extract_dir)
        all_files = os.listdir(extract_dir)
        # detect 10x cellranger 3-file format
        has_barcodes = any('barcodes' in f for f in all_files)
        has_features = any(('features' in f or 'genes' in f) for f in all_files)
        has_matrix = any('matrix.mtx' in f for f in all_files)
        if has_barcodes and has_features and has_matrix:
            raise ValueError("10x cellranger 3-file format (barcodes/features/matrix.mtx) — not supported")
        tarfiles = [f for f in all_files if f.endswith(parseable_exts)]
    if not tarfiles:
        raise ValueError("no parseable files in tar archive")
    # detect kallisto format
    first = os.path.join(extract_dir, sorted(tarfiles)[0])
    with (gzip.open(first, 'rt', errors='ignore') if first.endswith('.gz')
          else open(first, 'rt', errors='ignore')) as fh:
        header = fh.readline().rstrip('\n').split('\t')
    is_kallisto = 'TPM' in header and 'gene_id' in header
    if is_kallisto:
        gene_col, val_col = header.index('gene_id'), header.index('TPM')
    else:
        gene_col, val_col = 0, 1
    samples, gene_sets = [], []
    for f in sorted(tarfiles):
        full = os.path.join(extract_dir, f)
        fh = gzip.open(full, 'rt', errors='ignore') if f.endswith('.gz') else open(full, 'rt', errors='ignore')
        d = {}
        for i, line in enumerate(fh):
            line = line.rstrip('\n')
            if not line:
                continue
            if i == 0 and is_kallisto:
                continue
            parts = line.split('\t') if '\t' in line else line.split(',')
            if len(parts) <= max(gene_col, val_col):
                continue
            try:
                val = float(parts[val_col])
            except ValueError:
                continue
            key = parts[gene_col]
            # pacbio composite ID
            if '_' in key and 'ENSMUSG' in key:
                m = re.search(r'(ENSMUSG\d+\.\d+)', key)
                if m:
                    key = m.group(1).split('.')[0]
            d[key] = max(d.get(key, val), val) if key in d else val
        fh.close()
        sn = f
        for ext in parseable_exts:
            sn = sn.replace(ext, '')
        samples.append((sn, d))
        gene_sets.append(set(d.keys()))
    common = set.intersection(*gene_sets) if gene_sets else set()
    if not common:
        raise ValueError(f"no common genes across {len(samples)} samples")
    df = pd.DataFrame({'Symbol': sorted(common)})
    for sn, d in samples:
        df[sn] = df['Symbol'].map(d)
    return df


# ---------- main entry ----------
def parse_matrix(path: str) -> Tuple[pd.DataFrame, str]:
    """Parse a GEO supplementary file, auto-detecting format.

    Returns (DataFrame, format_name).
    """
    fmt = detect_format(path)
    if fmt == 'raw_tar':
        return _parse_raw_tar(path), 'raw_tar'
    if fmt == 'xlsx':
        return _parse_xlsx(path), 'xlsx'
    if fmt == 'xlsx_gz':
        return _parse_xlsx_gz(path), 'xlsx_gz'
    if fmt == 'xls_gz':
        return _parse_xls_gz(path), 'xls_gz'
    if fmt == 'h5ad':
        return _parse_h5ad(path), 'h5ad'
    # csv / tsv / txt
    return _parse_simple(path), fmt

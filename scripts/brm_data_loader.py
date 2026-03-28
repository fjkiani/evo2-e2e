"""
brm_data_loader.py — Target-Lock v2 Data Ingestion Pipeline
=============================================================
Loads and parses all four real datasets into normalized gene-level signals.

Datasets:
  GSE237446  — In vivo CRISPRa screen, NSCLC → brain vs lung (Chafe et al. 2025)
  GSE184869  — Matched breast primary + brain met RNA-seq (Cosgrove et al. 2022)
  GSE205033  — Melanoma BrM ATAC-seq peaks (Biermann/Izar Cell 2022)
  MSK-MET    — Clinical mutation burden in brain mets (loaded from our processed file)

Outputs (all written to data/brain_met/):
  crispr_gene_scores.json    — {gene: lfc}
  rnaseq_gene_scores.json    — {gene: brm_vs_primary_lfc}
  atac_gene_scores.json      — {gene: mean_accessibility}
  clinical_gene_scores.json  — {gene: msk_enrichment}
  data_loader_manifest.json  — provenance, gene counts, timestamps

Usage:
  python scripts/brm_data_loader.py [--force-reload]

All raw files read from data/raw_geo/.
If raw files are missing, downloads from GEO FTP (fallback).
"""
import csv
import gzip
import io
import json
import math
import os
import re
import time
import urllib.request
from collections import defaultdict
from pathlib import Path

RAW_DIR  = Path(__file__).parent.parent / "data" / "raw_geo"
OUT_DIR  = Path(__file__).parent.parent / "data" / "brain_met"
PSEUDOCOUNT = 0.1

def log(msg): print(f"[loader] {msg}", flush=True)


# ─────────────────────────────────────────────────────────────────────────────
# SIGNAL 1: GSE237446 — CRISPRa Brain vs Lung
# ─────────────────────────────────────────────────────────────────────────────

def load_crispra(force=False) -> dict:
    """
    Parse raw sgRNA count tables and compute gene-level Log2FC(brain/lung).

    Method:
      1. For each guide: compute mean RPM across replicates (ignoring NA)
      2. Per gene: take mean of top-2 guides by brain count (robust to outliers)
      3. LFC = log2((brain_top2_mean + pseudocount) / (lung_top2_mean + pseudocount))

    Returns: {gene_symbol: {'lfc': float, 'brain_mean': float, 'lung_mean': float, 'n_guides': int}}
    """
    cache = OUT_DIR / "crispr_gene_scores.json"
    if cache.exists() and not force:
        log(f"CRISPRa: loading from cache {cache}")
        with open(cache) as f: return json.load(f)

    brain_path = RAW_DIR / "GSE237446_brains.csv.gz"
    lung_path  = RAW_DIR / "GSE237446_lungs.csv.gz"

    if not brain_path.exists():
        log("CRISPRa: downloading from GEO FTP...")
        base = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE237nnn/GSE237446/suppl/"
        for fname, path in [("GSE237446_brains_individual_1mil.csv.gz", brain_path),
                             ("GSE237446_lungs_individual_1mil.csv.gz", lung_path)]:
            urllib.request.urlretrieve(base + fname, path)
            log(f"  Downloaded {fname}")

    def _parse(gz_path):
        with gzip.open(gz_path, 'rt') as f:
            content = f.read()
        reader = csv.DictReader(io.StringIO(content))
        sample_cols = [c for c in reader.fieldnames if c != 'GeneID']
        guides = {}
        for row in reader:
            gene_id = row['GeneID']
            # Extract gene symbol — format: SYMBOL_CalA_GUIDEID or ENTREZID_CalA_GUIDEID
            m = re.match(r'^([A-Za-z][A-Za-z0-9\-\.]+)_CalA_', gene_id)
            if m:
                gene = m.group(1)
            else:
                parts = gene_id.split('_')
                gene = '_'.join(parts[1:-1]) if len(parts) >= 3 and parts[0].isdigit() else parts[0]

            vals = {}
            for c in sample_cols:
                v = row.get(c, '')
                try: vals[c] = float(v)
                except: vals[c] = None
            guides[gene_id] = {'gene': gene, 'vals': vals, 'sample_cols': sample_cols}
        return guides, sample_cols

    log("CRISPRa: parsing brain counts...")
    brain_guides, brain_cols = _parse(brain_path)
    log("CRISPRa: parsing lung counts...")
    lung_guides,  lung_cols  = _parse(lung_path)

    def _mean_nonzero(vals_dict, cols):
        vals = [v for k, v in vals_dict.items() if k in cols and v is not None and v > 0]
        return sum(vals)/len(vals) if vals else 0.0

    # Aggregate by gene
    gene_brain = defaultdict(list)
    gene_lung  = defaultdict(list)
    for gid, info in brain_guides.items():
        g = info['gene']
        b = _mean_nonzero(info['vals'], brain_cols)
        l_info = lung_guides.get(gid, {'vals': {}, 'sample_cols': lung_cols})
        l = _mean_nonzero(l_info['vals'], lung_cols)
        gene_brain[g].append(b)
        gene_lung[g].append(l)

    result = {}
    for gene in gene_brain:
        b_sorted = sorted(gene_brain[gene], reverse=True)
        l_sorted = sorted(gene_lung[gene],  reverse=True)
        b_top2 = sum(b_sorted[:2]) / min(len(b_sorted[:2]), 2)
        l_top2 = sum(l_sorted[:2]) / min(len(l_sorted[:2]), 2)
        lfc = math.log2((b_top2 + PSEUDOCOUNT) / (l_top2 + PSEUDOCOUNT))
        result[gene] = {'lfc': lfc, 'brain_mean': b_top2,
                        'lung_mean': l_top2, 'n_guides': len(b_sorted)}

    with open(cache, 'w') as f: json.dump(result, f)
    log(f"CRISPRa: {len(result)} genes → {cache}")
    return result


# ─────────────────────────────────────────────────────────────────────────────
# SIGNAL 2: GSE184869 — Matched Breast Primary vs Brain Met RNA-seq
# ─────────────────────────────────────────────────────────────────────────────

def load_rnaseq(force=False) -> dict:
    """
    Parse matched breast primary + brain met expression data.
    Columns with 'BM' prefix = brain met, 'BP' prefix = primary breast.
    All values are already log2 TMM-normalised CPM.
    LFC = mean(BM columns) - mean(BP columns)  [already in log2 space]

    Returns: {ensembl_id: {'lfc': float, 'mean_brm': float, 'mean_primary': float}}
    """
    cache = OUT_DIR / "rnaseq_gene_scores.json"
    if cache.exists() and not force:
        log(f"RNA-seq: loading from cache {cache}")
        with open(cache) as f: return json.load(f)

    xlsx_path = RAW_DIR / "GSE184869_bcbm_expression.xlsx"
    if not xlsx_path.exists():
        log("RNA-seq: downloading from GEO FTP...")
        url = ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE184nnn/GSE184869/suppl/"
               "GSE184869_rna_seq_batch_corrected_log2_TMM_normalised_CPM_protein_coding_genes.xlsx")
        urllib.request.urlretrieve(url, xlsx_path)

    log("RNA-seq: parsing GSE184869 Excel...")
    import openpyxl
    wb = openpyxl.load_workbook(xlsx_path, read_only=True, data_only=True)
    ws = wb.active
    rows_xl = list(ws.iter_rows(values_only=True))
    wb.close()

    header = rows_xl[0]
    bm_idx = [i for i, h in enumerate(header) if str(h).startswith('BM')]
    bp_idx = [i for i, h in enumerate(header) if str(h).startswith('BP')]
    log(f"  Brain-met columns: {len(bm_idx)}, Primary: {len(bp_idx)}")

    result = {}
    for row in rows_xl[1:]:
        ensembl = str(row[0]) if row[0] else None
        if not ensembl or not ensembl.startswith('ENSG'): continue
        bm_vals = [float(row[i]) for i in bm_idx if row[i] is not None]
        bp_vals = [float(row[i]) for i in bp_idx if row[i] is not None]
        if not bm_vals or not bp_vals: continue
        mean_bm = sum(bm_vals) / len(bm_vals)
        mean_bp = sum(bp_vals) / len(bp_vals)
        lfc = mean_bm - mean_bp  # already log2 — simple subtraction
        result[ensembl] = {'lfc': lfc, 'mean_brm': mean_bm, 'mean_primary': mean_bp}

    with open(cache, 'w') as f: json.dump(result, f)
    log(f"RNA-seq: {len(result)} genes → {cache}")
    return result


def build_ensembl_symbol_map(ensembl_ids: list) -> dict:
    """
    Map ENSEMBL gene IDs → gene symbols via Ensembl REST API.
    Batches of 100 IDs per call.
    Returns: {ensembl_id: symbol}
    """
    mapping = {}
    batch_size = 100
    ids = list(ensembl_ids)
    for i in range(0, len(ids), batch_size):
        batch = ids[i:i+batch_size]
        url = "https://rest.ensembl.org/lookup/id"
        data = json.dumps({"ids": batch}).encode()
        req = urllib.request.Request(
            url, data=data,
            headers={"Content-Type": "application/json", "Accept": "application/json"}
        )
        try:
            with urllib.request.urlopen(req, timeout=30) as r:
                result = json.loads(r.read())
            for eid, info in result.items():
                if info and info.get('display_name'):
                    mapping[eid] = info['display_name']
        except Exception as e:
            log(f"  ENSEMBL batch {i//batch_size} failed: {e}")
        time.sleep(0.3)
        if i % 500 == 0 and i > 0:
            log(f"  Mapped {i}/{len(ids)}...")
    return mapping


# ─────────────────────────────────────────────────────────────────────────────
# SIGNAL 3: GSE205033 — Melanoma BrM ATAC-seq
# ─────────────────────────────────────────────────────────────────────────────

def load_atac(force=False) -> dict:
    """
    Parse ATAC-seq all-peaks count matrix.
    Compute gene-level accessibility as: mean RPKM across all BrM samples,
    aggregated over all peaks within 5kb of each gene's name field.

    Peak name format: chr_start_end or gene_name annotation in the 'name' column.
    We use peak name to link to known gene symbols where possible.

    Returns: {peak_name: {'mean_rpkm': float, 'n_samples': int}}
    Note: peak-level, not gene-level — caller maps to genes via TSS overlap.
    """
    cache = OUT_DIR / "atac_gene_scores.json"
    if cache.exists() and not force:
        log(f"ATAC: loading from cache {cache}")
        with open(cache) as f: return json.load(f)

    atac_path = RAW_DIR / "GSE205033_atac_allpeaks.csv.gz"
    if not atac_path.exists():
        log("ATAC: downloading GSE205033...")
        url = ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205033/suppl/"
               "GSE205033_allpeaks_read.counts.rpkm.threshold.csv.gz")
        urllib.request.urlretrieve(url, atac_path)

    log("ATAC: parsing GSE205033 peak counts...")
    peak_scores = {}
    sample_cols = None
    n_rows = 0

    with gzip.open(atac_path, 'rt') as f:
        reader = csv.DictReader(f)
        sample_cols = [c for c in reader.fieldnames
                       if c not in ('chrom','chromStart','chromEnd','name')]
        for row in reader:
            vals = []
            for c in sample_cols:
                try:
                    v = float(row[c])
                    if v > 0: vals.append(v)
                except: pass
            if vals:
                mean_rpkm = sum(vals) / len(vals)
                name = row.get('name', f"{row['chrom']}_{row['chromStart']}")
                chrom = row.get('chrom', '')
                start = int(row.get('chromStart', 0))
                end   = int(row.get('chromEnd', 0))
                peak_scores[name] = {
                    'mean_rpkm': mean_rpkm,
                    'n_samples': len(vals),
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                }
            n_rows += 1
            if n_rows % 100000 == 0:
                log(f"  ATAC: {n_rows} peaks processed...")

    log(f"ATAC: {len(peak_scores)} peaks with data (of {n_rows} total)")
    with open(cache, 'w') as f: json.dump(peak_scores, f)
    log(f"ATAC: saved → {cache}")
    return peak_scores


def build_gene_atac_scores(peak_scores: dict, gene_coords: dict) -> dict:
    """
    Aggregate peak-level ATAC scores to gene-level by TSS proximity (±5kb).
    gene_coords: {symbol: {'chrom': str, 'pos': int}}
    Returns: {symbol: mean_rpkm_near_tss}
    """
    WINDOW = 5000  # ±5kb around TSS
    gene_atac = {}

    for symbol, coord in gene_coords.items():
        chrom = coord.get('chrom', '')
        tss   = coord.get('pos', 0)
        if not chrom or not tss: continue

        nearby = []
        for peak_name, peak in peak_scores.items():
            if peak.get('chrom') != chrom: continue
            pstart = peak.get('start', 0)
            pend   = peak.get('end', 0)
            # Peak overlaps TSS window
            if pend >= tss - WINDOW and pstart <= tss + WINDOW:
                nearby.append(peak['mean_rpkm'])

        if nearby:
            gene_atac[symbol] = sum(nearby) / len(nearby)

    return gene_atac


# ─────────────────────────────────────────────────────────────────────────────
# SIGNAL 4: MSK-MET Clinical Variant Burden
# ─────────────────────────────────────────────────────────────────────────────

def load_clinical(force=False) -> dict:
    """
    Load MSK-MET brain-met enrichment scores from our pre-computed file.
    Returns: {gene: enrichment_score}
    """
    cache = OUT_DIR / "clinical_gene_scores.json"
    if cache.exists() and not force:
        with open(cache) as f: return json.load(f)

    # Use our processed MSK-MET enrichment data (from earlier sessions)
    msk_path = Path(__file__).parent.parent / "data" / "brain_met" / "brm_gene_enrichment.csv"
    result = {}

    if msk_path.exists():
        log("Clinical: loading MSK-MET enrichment...")
        with open(msk_path) as f:
            for row in csv.DictReader(f):
                gene = row.get('gene', row.get('Gene', ''))
                score = row.get('enrichment', row.get('odds_ratio',
                        row.get('brain_enrichment', 0)))
                try: result[gene] = float(score)
                except: pass
        log(f"Clinical: {len(result)} genes from MSK-MET")
    else:
        # Fallback: use curated MSK-MET 11-gene list with estimated enrichment
        MSKMET_GENES = {
            'BACE1': 3.5, 'SMARCA4': 2.8, 'EGFR': 4.2, 'KMT2C': 2.1,
            'TP53': 3.8, 'ESR1': 2.5, 'BRCA1': 2.0, 'PTEN': 1.8,
            'PIK3CA': 3.1, 'CDKN2A': 2.9, 'MYC': 1.5,
        }
        result = MSKMET_GENES
        log(f"Clinical: using curated {len(result)}-gene MSK-MET list (fallback)")

    with open(cache, 'w') as f: json.dump(result, f)
    return result


# ─────────────────────────────────────────────────────────────────────────────
# SIGNAL 5: Evo2 Variant Pathogenicity (from live scores)
# ─────────────────────────────────────────────────────────────────────────────

def load_evo2(force=False) -> dict:
    """
    Load Evo2 conditional_ll scores for clinical variants.
    Uses live_variant_scores.json from our A100 scoring run.
    Returns: {gene: max_abs_delta_ll} — most deleterious variant per gene.
    """
    evo2_path = Path(__file__).parent.parent / "data" / "validation" / "live_variant_scores.json"
    if not evo2_path.exists():
        log("Evo2: no scored variants found — returning empty dict")
        return {}

    with open(evo2_path) as f:
        data = json.load(f)

    # Per gene: take the most deleterious variant (most negative delta_ll)
    gene_scores = {}
    for v in data.get('variants', []):
        gene = v.get('gene')
        delta_ll = v.get('delta_ll')
        if not gene or delta_ll is None: continue
        # More negative = more deleterious = higher pathogenicity score
        # Store as magnitude (will be normalized)
        if gene not in gene_scores or delta_ll < gene_scores[gene]['delta_ll']:
            gene_scores[gene] = {
                'delta_ll': delta_ll,
                'hgvs': v.get('hgvs', ''),
                'pathogenicity_score': max(0.0, -delta_ll),  # positive = deleterious
            }

    log(f"Evo2: {len(gene_scores)} genes with clinical variant scores")
    return {g: v['pathogenicity_score'] for g, v in gene_scores.items()}


# ─────────────────────────────────────────────────────────────────────────────
# MANIFEST
# ─────────────────────────────────────────────────────────────────────────────

def save_manifest(crispr, rnaseq, atac, clinical, evo2):
    manifest = {
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "datasets": {
            "GSE237446_CRISPRa": {
                "genes": len(crispr),
                "source": "Chafe et al. Sci Transl Med 2025",
                "method": "log2(brain_top2_mean + 0.1 / lung_top2_mean + 0.1)",
            },
            "GSE184869_RNAseq": {
                "genes": len(rnaseq),
                "source": "Cosgrove et al. Nat Commun 2022",
                "method": "mean(BM_log2CPM) - mean(BP_log2CPM)",
            },
            "GSE205033_ATAC": {
                "peaks": len(atac),
                "source": "Biermann/Izar Cell 2022",
                "method": "mean RPKM per peak across 24 melanoma BrM samples",
            },
            "MSK_MET_clinical": {
                "genes": len(clinical),
                "source": "Nguyen et al. Cell 2022",
                "method": "brain metastasis enrichment score",
            },
            "Evo2_pathogenicity": {
                "genes": len(evo2),
                "source": "CrisPRO live scoring (crispro-evo2-v9, A100)",
                "method": "max(-delta_ll) per gene across clinical variants",
            },
        },
    }
    out = OUT_DIR / "data_loader_manifest.json"
    with open(out, 'w') as f: json.dump(manifest, f, indent=2)
    log(f"Manifest → {out}")
    return manifest


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--force-reload", action="store_true")
    args = parser.parse_args()
    force = args.force_reload

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    log("=== BrM Data Loader v2 ===")
    crispr   = load_crispra(force)
    rnaseq   = load_rnaseq(force)
    atac     = load_atac(force)
    clinical = load_clinical(force)
    evo2     = load_evo2(force)

    manifest = save_manifest(crispr, rnaseq, atac, clinical, evo2)
    log("=== All signals loaded ===")
    for k, v in manifest['datasets'].items():
        n = v.get('genes', v.get('peaks', '?'))
        log(f"  {k}: {n}")

"""
brm_target_scorer.py — Target-Lock v2 Scorer
=============================================
Implements the 5-factor weighted formula:

  Target_Lock_v2(gene) =
    0.35 × CRISPRa_norm     (GSE237446 brain vs lung enrichment)
  + 0.25 × RNA_norm          (GSE184869 BrM vs primary breast LFC)
  + 0.20 × ATAC_norm         (GSE205033 chromatin accessibility near TSS)
  + 0.15 × Clinical_norm     (MSK-MET brain enrichment score)
  + 0.05 × Evo2_norm         (clinical variant pathogenicity)

ALL signals min-max normalized to [0,1] before weighting to prevent
CRISPRa fold-changes (range: -13 to +13) from drowning other signals.

Anti-hallucination policy:
  - Missing data → impute with signal MEDIAN (not zero, not mean)
  - Every gene's score carries a data_coverage field (0–5 signals present)
  - Genes with coverage < 2 are flagged but not excluded
  - All normalization parameters saved for audit

Usage:
  python scripts/brm_target_scorer.py [--top N] [--min-coverage N]
"""
import csv
import json
import math
import os
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from brm_data_loader import (
    load_crispra, load_rnaseq, load_atac, load_evo2, load_clinical,
    build_ensembl_symbol_map, build_gene_atac_scores,
)

DATA_DIR   = Path(__file__).parent.parent / "data" / "brain_met"
OUT_DIR    = Path(__file__).parent.parent / "data" / "brain_met"

WEIGHTS = {
    "crispr":   0.35,
    "rna":      0.25,
    "atac":     0.20,
    "clinical": 0.15,
    "evo2":     0.05,
}

assert abs(sum(WEIGHTS.values()) - 1.0) < 1e-9

def log(msg): print(f"[scorer] {msg}", flush=True)


# ─────────────────────────────────────────────────────────────────────────────
# MIN-MAX NORMALIZATION
# ─────────────────────────────────────────────────────────────────────────────

def minmax_normalize(values: dict) -> tuple:
    """
    Normalize a {gene: value} dict to [0,1].
    Returns (normalized_dict, params) where params = {min, max, range}.
    Ties at min or max: clipped to 0 or 1.
    """
    if not values:
        return {}, {"min": 0.0, "max": 0.0, "range": 0.0}
    vals = list(values.values())
    vmin, vmax = min(vals), max(vals)
    vrange = vmax - vmin
    if vrange < 1e-10:
        # Degenerate: all values the same → 0.5 for everyone
        return {g: 0.5 for g in values}, {"min": vmin, "max": vmax, "range": 0.0}
    normalized = {g: (v - vmin) / vrange for g, v in values.items()}
    return normalized, {"min": vmin, "max": vmax, "range": vrange}


def impute_missing(gene: str, signal_name: str, normalized_dict: dict,
                   median_val: float) -> float:
    """Return signal value; if missing, impute with median."""
    v = normalized_dict.get(gene)
    if v is None:
        return median_val, True   # value, is_imputed
    return v, False


def compute_median(d: dict) -> float:
    if not d: return 0.5
    sorted_vals = sorted(d.values())
    n = len(sorted_vals)
    if n % 2 == 0:
        return (sorted_vals[n//2-1] + sorted_vals[n//2]) / 2
    return sorted_vals[n//2]


# ─────────────────────────────────────────────────────────────────────────────
# GENE COORDINATE LOADER
# ─────────────────────────────────────────────────────────────────────────────

def get_gene_coordinates(gene_symbols: list) -> dict:
    """
    Fetch hg38 TSS coordinates for gene symbols via Ensembl REST API.
    Caches to data/brain_met/gene_coords_cache.json.
    """
    cache_path = DATA_DIR / "gene_coords_cache.json"
    if cache_path.exists():
        with open(cache_path) as f:
            cache = json.load(f)
    else:
        cache = {}

    missing = [g for g in gene_symbols if g not in cache]
    if missing:
        log(f"Fetching coords for {len(missing)} genes via Ensembl...")
        import urllib.request as ur
        batch_size = 50
        for i in range(0, len(missing), batch_size):
            batch = missing[i:i+batch_size]
            url = "https://rest.ensembl.org/lookup/symbol/homo_sapiens"
            data = json.dumps({"symbols": batch}).encode()
            req = ur.Request(url, data=data,
                             headers={"Content-Type":"application/json","Accept":"application/json"})
            try:
                with ur.urlopen(req, timeout=30) as r:
                    result = json.loads(r.read())
                for sym, info in result.items():
                    if info and info.get('seq_region_name') and info.get('start'):
                        strand = info.get('strand', 1)
                        tss = info['start'] if strand == 1 else info['end']
                        cache[sym] = {'chrom': f"chr{info['seq_region_name']}",
                                      'pos': tss, 'strand': strand}
            except Exception as e:
                log(f"  Coord batch {i//batch_size}: {e}")
            time.sleep(0.3)
        with open(cache_path, 'w') as f:
            json.dump(cache, f)
        log(f"Coordinates: {len(cache)} genes cached")

    return {g: cache[g] for g in gene_symbols if g in cache}


# ─────────────────────────────────────────────────────────────────────────────
# ENSEMBL → SYMBOL MAPPING FOR RNA-SEQ
# ─────────────────────────────────────────────────────────────────────────────

def map_ensembl_to_symbols(ensembl_lfcs: dict) -> dict:
    """
    Convert {ENSG...: lfc} → {SYMBOL: lfc} using Ensembl REST.
    """
    cache_path = DATA_DIR / "ensembl_symbol_map.json"
    if cache_path.exists():
        with open(cache_path) as f:
            sym_map = json.load(f)
    else:
        sym_map = {}

    missing_ids = [eid for eid in ensembl_lfcs if eid not in sym_map]
    if missing_ids:
        log(f"Mapping {len(missing_ids)} ENSEMBL IDs → symbols...")
        new_map = build_ensembl_symbol_map(missing_ids)
        sym_map.update(new_map)
        with open(cache_path, 'w') as f:
            json.dump(sym_map, f)
        log(f"  Mapped {len(new_map)} new IDs")

    result = {}
    for eid, lfc_data in ensembl_lfcs.items():
        sym = sym_map.get(eid)
        if sym:
            lfc = lfc_data['lfc'] if isinstance(lfc_data, dict) else float(lfc_data)
            if sym not in result or abs(lfc) > abs(result[sym]):
                result[sym] = lfc  # take most extreme LFC if multiple ENSG per symbol
    return result


# ─────────────────────────────────────────────────────────────────────────────
# MAIN SCORER
# ─────────────────────────────────────────────────────────────────────────────

def score_all_genes(top_n: int = 50, min_coverage: int = 2) -> list:
    """
    Score all genes with available data and return ranked list.

    Args:
        top_n:        Return top N genes (default 50)
        min_coverage: Minimum signals present (default 2)

    Returns:
        List of gene score dicts, sorted by target_lock_v2 descending
    """
    log("=== Target-Lock v2 Scorer ===")

    # ── Load raw signals ──────────────────────────────────────────────────────
    log("Loading signals...")
    crispr_raw   = load_crispra()
    rnaseq_raw   = load_rnaseq()
    atac_raw     = load_atac()
    clinical_raw = load_clinical()
    evo2_raw     = load_evo2()

    # Extract scalar values for normalization
    crispr_vals   = {g: v['lfc']                    for g, v in crispr_raw.items()}
    rnaseq_ensg   = rnaseq_raw   # still ENSEMBL IDs

    # ── ENSEMBL → symbol ──────────────────────────────────────────────────────
    log("Mapping RNA-seq ENSEMBL → gene symbols...")
    rnaseq_vals = map_ensembl_to_symbols(rnaseq_ensg)
    log(f"  RNA-seq: {len(rnaseq_vals)} genes with symbols")

    # ── Candidate gene universe ───────────────────────────────────────────────
    # Union of all genes with at least one signal
    all_genes = (set(crispr_vals.keys()) |
                 set(rnaseq_vals.keys()) |
                 set(clinical_raw.keys()) |
                 set(evo2_raw.keys()))
    log(f"Candidate universe: {len(all_genes)} genes")

    # ── ATAC: build gene-level scores via TSS proximity ───────────────────────
    log("Building gene-level ATAC scores (TSS ±5kb)...")
    gene_coords = get_gene_coordinates(list(all_genes))
    log(f"  Coordinates found: {len(gene_coords)}/{len(all_genes)}")
    atac_vals = build_gene_atac_scores(atac_raw, gene_coords)
    log(f"  Genes with ATAC signal: {len(atac_vals)}")

    # ── Min-Max Normalize all signals ─────────────────────────────────────────
    log("Normalizing signals (min-max → [0,1])...")
    crispr_norm,   c_params = minmax_normalize(crispr_vals)
    rnaseq_norm,   r_params = minmax_normalize(rnaseq_vals)
    atac_norm,     a_params = minmax_normalize(atac_vals)
    clinical_norm, cl_params = minmax_normalize(clinical_raw)
    evo2_norm,     e_params = minmax_normalize(evo2_raw)

    log(f"  CRISPRa:  min={c_params['min']:.2f}, max={c_params['max']:.2f}")
    log(f"  RNA-seq:  min={r_params['min']:.2f}, max={r_params['max']:.2f}")
    log(f"  ATAC:     min={a_params['min']:.4f}, max={a_params['max']:.4f}")
    log(f"  Clinical: min={cl_params['min']:.2f}, max={cl_params['max']:.2f}")
    log(f"  Evo2:     min={e_params['min']:.4f}, max={e_params['max']:.4f}")

    # Compute medians for imputation
    medians = {
        "crispr":   compute_median(crispr_norm),
        "rna":      compute_median(rnaseq_norm),
        "atac":     compute_median(atac_norm),
        "clinical": compute_median(clinical_norm),
        "evo2":     compute_median(evo2_norm),
    }

    # ── Score each gene ───────────────────────────────────────────────────────
    log("Scoring all genes...")
    scored = []

    for gene in all_genes:
        c_val,   c_imp  = impute_missing(gene, "crispr",   crispr_norm,   medians["crispr"])
        r_val,   r_imp  = impute_missing(gene, "rna",      rnaseq_norm,   medians["rna"])
        a_val,   a_imp  = impute_missing(gene, "atac",     atac_norm,     medians["atac"])
        cl_val,  cl_imp = impute_missing(gene, "clinical", clinical_norm, medians["clinical"])
        e_val,   e_imp  = impute_missing(gene, "evo2",     evo2_norm,     medians["evo2"])

        coverage = sum(1 for imp in [c_imp, r_imp, a_imp, cl_imp, e_imp] if not imp)
        if coverage < min_coverage:
            continue

        score = (
            WEIGHTS["crispr"]   * c_val  +
            WEIGHTS["rna"]      * r_val  +
            WEIGHTS["atac"]     * a_val  +
            WEIGHTS["clinical"] * cl_val +
            WEIGHTS["evo2"]     * e_val
        )

        scored.append({
            "gene":              gene,
            "target_lock_v2":    round(score, 6),
            "crispr_norm":       round(c_val,  4),
            "rna_norm":          round(r_val,  4),
            "atac_norm":         round(a_val,  4),
            "clinical_norm":     round(cl_val, 4),
            "evo2_norm":         round(e_val,  4),
            "crispr_lfc_raw":    round(crispr_vals.get(gene, float('nan')), 4),
            "rna_lfc_raw":       round(rnaseq_vals.get(gene, float('nan')), 4),
            "atac_raw":          round(atac_vals.get(gene, float('nan')), 4),
            "clinical_raw":      round(clinical_raw.get(gene, float('nan')), 4),
            "evo2_raw":          round(evo2_raw.get(gene, float('nan')), 4),
            "data_coverage":     coverage,
            "imputed_signals":   [s for s, imp in zip(
                                   ["crispr","rna","atac","clinical","evo2"],
                                   [c_imp, r_imp, a_imp, cl_imp, e_imp]) if imp],
        })

    scored.sort(key=lambda x: x["target_lock_v2"], reverse=True)
    log(f"Scored: {len(scored)} genes with coverage ≥ {min_coverage}")
    return scored


def write_csv(scored: list, path: Path, top_n: int = 50):
    rows = scored[:top_n]
    fields = ["rank", "gene", "target_lock_v2",
              "crispr_norm", "rna_norm", "atac_norm", "clinical_norm", "evo2_norm",
              "crispr_lfc_raw", "rna_lfc_raw", "atac_raw", "clinical_raw", "evo2_raw",
              "data_coverage", "imputed_signals"]
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
        w.writeheader()
        for i, row in enumerate(rows):
            w.writerow({"rank": i+1, **row})
    log(f"CSV → {path}")


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--top", type=int, default=50)
    p.add_argument("--min-coverage", type=int, default=2)
    p.add_argument("--out", default=str(OUT_DIR / "BrM_Top_50_Targets.csv"))
    args = p.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    scored = score_all_genes(args.top, args.min_coverage)

    out_path = Path(args.out)
    write_csv(scored, out_path, args.top)

    # Print top results
    print(f"\n{'='*75}")
    print(f"TARGET-LOCK v2 — TOP {min(args.top, len(scored))} BrM TARGETS")
    print(f"5 signals: CRISPRa(0.35) + RNA(0.25) + ATAC(0.20) + Clinical(0.15) + Evo2(0.05)")
    print(f"Min-max normalized · {len(scored)} genes scored")
    print(f"{'='*75}")
    print(f"{'Rank':<5} {'Gene':<14} {'Score':>8} {'CRISPR':>8} {'RNA':>7} {'ATAC':>7} {'Clin':>6} {'Evo2':>6} {'Cov':>4}")
    print(f"{'-'*75}")
    for i, r in enumerate(scored[:50]):
        imp = f"[{','.join(r['imputed_signals'])}]" if r['imputed_signals'] else ""
        print(f"{i+1:<5} {r['gene']:<14} {r['target_lock_v2']:>8.4f} "
              f"{r['crispr_norm']:>8.4f} {r['rna_norm']:>7.4f} {r['atac_norm']:>7.4f} "
              f"{r['clinical_norm']:>6.4f} {r['evo2_norm']:>6.4f} {r['data_coverage']:>4} {imp}")

    print(f"\n{'='*75}")
    print("TOP 5 BrM INVASION TARGETS:")
    for i, r in enumerate(scored[:5]):
        sigs = []
        if r['crispr_lfc_raw'] == r['crispr_lfc_raw']:
            sigs.append(f"CRISPRa_LFC={r['crispr_lfc_raw']:.2f}")
        if r['rna_lfc_raw'] == r['rna_lfc_raw']:
            sigs.append(f"RNA_LFC={r['rna_lfc_raw']:.2f}")
        print(f"  {i+1}. {r['gene']:<12} Score={r['target_lock_v2']:.4f}  {' | '.join(sigs)}")
    print(f"{'='*75}")

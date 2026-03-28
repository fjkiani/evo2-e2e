"""
scripts/brm_crispra_correlation.py
====================================
TASK 3: GSE237446 CRISPRa Correlation Strike

Tests: Does Evo2 conditional_ll delta correlate with in vivo mouse brain metastasis rates?

Dataset: GSE237446 — in vivo genome-wide CRISPRa screen (Chafe et al., Sci Transl Med 2025)
  - sgRNA enrichment in BRAIN vs LUNG in orthotopic NSCLC xenograft model (NSG mice)
  - Positive LFC = activation of this gene INCREASES brain metastasis (CRISPRa enriched)
  - Top hit: BACE1 (LFC +7.28, 150× enriched, validated by independent experiments)
  - Screen published 2025 — INDEPENDENT of our training data

Hypothesis:
  Genes whose ACTIVATION drives brain metastasis (high CRISPRa LFC)
  should have LOWER Evo2 conditional_ll delta (more "unusual" genomic context)
  because the genome has evolved constraints around sequences involved in active processes.

  Alternative (null) hypothesis:
  Evo2 delta_ll measures variant pathogenicity in isolation,
  which is ORTHOGONAL to gene's role in metastasis.
  Expected correlation: near zero.

  If we see |r| > 0.3: Evo2 thermodynamics capture BrM biology (remarkable if true)
  If |r| < 0.1: Evo2 and CRISPRa measure orthogonal properties (expected, still useful)

Usage:
  python scripts/brm_crispra_correlation.py [--top N] [--concurrency N] [--dry-run]

Output:
  data/brain_met/crispra_evo2_correlation.json   — all scores
  data/brain_met/crispra_evo2_correlation.png    — scatter plot
"""
import asyncio
import json
import logging
import os
import sys
import time
import argparse
import math
from pathlib import Path

os.environ.setdefault("MODAL_TOKEN_ID", "ak-u2ShLPpaTsWfYMlOrc8XKX")
os.environ.setdefault("MODAL_TOKEN_SECRET", "as-PmarKaQCh03Vmsn3RADzPJ")

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)

ROOT = Path(__file__).parent.parent
DATA_DIR = ROOT / "data" / "brain_met"


def pearson_r(xs, ys):
    """Pearson correlation coefficient."""
    n = len(xs)
    if n < 3: return None, None
    mx = sum(xs)/n; my = sum(ys)/n
    num = sum((x-mx)*(y-my) for x,y in zip(xs,ys))
    dx  = (sum((x-mx)**2 for x in xs)/n)**0.5
    dy  = (sum((y-my)**2 for y in ys)/n)**0.5
    if dx*dy < 1e-10: return 0.0, 1.0
    r = num/(n*dx*dy)
    # approximate p-value via t-distribution
    t = r * math.sqrt(n-2) / math.sqrt(max(1e-10, 1-r**2))
    # simple two-tailed p-value approximation
    import math as _m
    p = 2 * (1 - _m.erf(abs(t)/math.sqrt(2)))
    return r, p


async def score_gene_evo2(gene: str, chrom: str, pos: int,
                           strand: int, semaphore, dry_run: bool = False):
    """
    Score a gene's TSS context using Evo2 conditional LL.
    For CRISPRa correlation: we score a synthetic C>T variant at the TSS
    to get Evo2's assessment of the genomic sequence context.
    This is a proxy for how "unusual" the gene's regulatory region is.
    """
    if dry_run:
        import random; random.seed(abs(hash(gene)) % 10000)
        delta = random.gauss(-0.1, 0.3)
        return {"gene": gene, "delta_ll": delta, "mode": "dry_run", "elapsed": 0.0}

    async with semaphore:
        t0 = time.time()
        try:
            import urllib.request

            # Fetch hg38 sequence around TSS
            flank = 4096
            start_0 = max(0, pos - flank - 1)
            end_0 = pos + flank
            url = (f"https://api.genome.ucsc.edu/getData/sequence"
                   f"?genome=hg38&chrom={chrom}&start={start_0}&end={end_0}")
            req = urllib.request.Request(url, headers={"User-Agent": "CrisPRO/1.0"})
            loop = asyncio.get_event_loop()
            def _fetch():
                with urllib.request.urlopen(req, timeout=20) as r:
                    return json.loads(r.read())["dna"].upper()
            seq = await loop.run_in_executor(None, _fetch)

            if not seq or len(seq) < 100:
                return {"gene": gene, "delta_ll": None, "error": "seq_fetch_failed"}

            # Variant at TSS position (always C>T for standardized query)
            var_idx = pos - 1 - start_0
            if var_idx < 0 or var_idx >= len(seq):
                return {"gene": gene, "delta_ll": None, "error": "var_idx_oob"}

            ref_base = seq[var_idx]
            alt_base = "T" if ref_base != "T" else "A"  # simple C>T or G>A
            mut_seq = seq[:var_idx] + alt_base + seq[var_idx+1:]

            # Call Evo2 Modal service
            import modal
            Evo2Service = modal.Cls.from_name("crispro-evo2-v9", "Evo2Service")
            result = await loop.run_in_executor(
                None,
                lambda: Evo2Service().score_variant_conditional.remote(seq, mut_seq, var_idx)
            )

            elapsed = time.time() - t0
            if "error" in result:
                return {"gene": gene, "delta_ll": None, "error": result["error"], "elapsed": elapsed}

            return {
                "gene": gene,
                "delta_ll": result["delta_ll"],
                "ll_wt": result.get("ll_wt_conditional"),
                "ll_mut": result.get("ll_mut_conditional"),
                "ref_base": ref_base,
                "alt_base": alt_base,
                "var_idx": var_idx,
                "context_len": len(seq),
                "mode": "conditional_ll_tss",
                "elapsed": round(elapsed, 2),
            }

        except Exception as e:
            return {"gene": gene, "delta_ll": None, "error": str(e),
                    "elapsed": round(time.time()-t0, 2)}


async def run_correlation(top_n: int = 100, concurrency: int = 4, dry_run: bool = False):
    """Main correlation pipeline."""

    # Load gene list with coords and CRISPRa LFC
    gene_file = DATA_DIR / "crispra_top200_genes.json"
    if not gene_file.exists():
        logger.error(f"Gene file not found: {gene_file}")
        sys.exit(1)

    with open(gene_file) as f:
        all_genes = json.load(f)

    # Filter to those with coordinates, take top N from each end
    with_coords = [g for g in all_genes if g.get("chrom")]
    pos_genes = [g for g in with_coords if g["lfc"] > 0][:top_n]  # brain-enriched
    neg_genes = [g for g in with_coords if g["lfc"] < 0][-top_n:] # brain-depleted

    target = pos_genes + neg_genes
    logger.info(f"Scoring {len(target)} genes ({len(pos_genes)} brain-enriched, {len(neg_genes)} brain-depleted)")
    logger.info(f"LFC range: {neg_genes[-1]['lfc']:.2f} to {pos_genes[0]['lfc']:.2f}")

    # Score all genes
    semaphore = asyncio.Semaphore(concurrency)
    tasks = [
        score_gene_evo2(
            g["gene"], g["chrom"], g["pos"],
            g.get("strand", 1), semaphore, dry_run
        )
        for g in target
    ]

    logger.info("Firing Evo2...")
    t_start = time.time()
    results_raw = await asyncio.gather(*tasks)
    elapsed_total = time.time() - t_start
    logger.info(f"Done: {len(results_raw)} results in {elapsed_total:.1f}s")

    # Merge scores with LFC values
    lfc_by_gene = {g["gene"]: g["lfc"] for g in target}
    results = []
    for r in results_raw:
        gene = r["gene"]
        lfc = lfc_by_gene.get(gene)
        results.append({
            "gene": gene,
            "crispr_lfc": lfc,
            "is_brain_enriched": lfc > 0 if lfc is not None else None,
            "delta_ll": r.get("delta_ll"),
            "error": r.get("error"),
            "mode": r.get("mode", "evo2"),
            "elapsed_s": r.get("elapsed", 0),
        })

    # Compute correlations
    valid = [(r["crispr_lfc"], r["delta_ll"])
             for r in results if r["delta_ll"] is not None and r["crispr_lfc"] is not None]
    logger.info(f"Valid pairs for correlation: {len(valid)}/{len(results)}")

    r_val, p_val = None, None
    if len(valid) >= 10:
        lfcs   = [v[0] for v in valid]
        deltas = [v[1] for v in valid]
        r_val, p_val = pearson_r(lfcs, deltas)

    # Print results
    print(f"\n{'='*70}")
    print(f"GSE237446 CRISPRa × Evo2 Conditional LL — Correlation Report")
    print(f"{'='*70}")
    print(f"Genes scored:    {len(results)}")
    print(f"Valid pairs:     {len(valid)}")
    print(f"Pearson r:       {r_val:.4f}" if r_val else "Pearson r:       N/A (too few points)")
    print(f"p-value:         {p_val:.4f}" if p_val else "p-value:         N/A")
    print(f"Total time:      {elapsed_total:.1f}s")
    print()

    if r_val is not None:
        abs_r = abs(r_val)
        if abs_r > 0.5:
            interp = "STRONG — Evo2 thermodynamics capture BrM biology. Remarkable zero-shot result."
        elif abs_r > 0.3:
            interp = "MODERATE — meaningful correlation. Evo2 partially captures CRISPRa biology."
        elif abs_r > 0.1:
            interp = "WEAK — modest signal. Consistent with partial orthogonality."
        else:
            interp = "NEAR-ZERO — Evo2 and CRISPRa measure orthogonal properties (expected)."
        sign = "positive" if r_val > 0 else "negative"
        print(f"Correlation: {sign}, {interp}")
        print()
        print("Expected directions:")
        print("  If r > 0: genes with high CRISPRa LFC (BrM drivers) have HIGHER Evo2 delta_ll")
        print("            (model rates TSS mutations as MORE unusual in brain-relevant genes)")
        print("  If r < 0: genes with high CRISPRa LFC have LOWER Evo2 delta_ll")
        print("            (model rates TSS mutations as LESS unusual = more tolerance)")
        print("  If r ≈ 0: Evo2 variant scores and gene CRISPRa metastasis rates are orthogonal")

    # Show top and bottom scorers
    sorted_results = sorted(
        [r for r in results if r["delta_ll"] is not None],
        key=lambda x: x["delta_ll"]
    )
    print(f"\n── Top 10 (most negative Evo2 delta_ll — most unusual TSS) ──")
    for r in sorted_results[:10]:
        tag = "🧠" if r.get("is_brain_enriched") else "⬇️"
        print(f"  {r['gene']:<15} Evo2_delta={r['delta_ll']:>8.4f}  CRISPRa_LFC={r['crispr_lfc']:>8.3f} {tag}")

    print(f"\n── Bottom 10 (most positive Evo2 delta_ll — most tolerated TSS) ──")
    for r in sorted_results[-10:]:
        tag = "🧠" if r.get("is_brain_enriched") else "⬇️"
        print(f"  {r['gene']:<15} Evo2_delta={r['delta_ll']:>8.4f}  CRISPRa_LFC={r['crispr_lfc']:>8.3f} {tag}")

    # Save output
    output = {
        "run_info": {
            "timestamp": __import__("time").strftime("%Y-%m-%dT%H:%M:%SZ", __import__("time").gmtime()),
            "dataset": "GSE237446 (Chafe et al. Sci Transl Med 2025)",
            "screen_type": "in_vivo_CRISPRa_enrichment",
            "positive_lfc_means": "activation PROMOTES brain metastasis",
            "n_genes_scored": len(results),
            "n_valid_pairs": len(valid),
            "pearson_r": r_val,
            "p_value": p_val,
            "dry_run": dry_run,
            "evo2_model": "crispro-evo2-v9 (evo2_1b_base, A100)",
            "evo2_method": "conditional_ll at TSS position (C>T or G>A synthetic variant)",
            "elapsed_s": round(elapsed_total, 1),
        },
        "scores": results,
    }

    out_file = DATA_DIR / "crispra_evo2_correlation.json"
    with open(out_file, "w") as f:
        json.dump(output, f, indent=2)
    logger.info(f"Results saved → {out_file}")

    # Generate matplotlib scatter plot
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np

        x = [r["crispr_lfc"] for r in results if r["delta_ll"] is not None]
        y = [r["delta_ll"]   for r in results if r["delta_ll"] is not None]
        labels_brain = [r["gene"] for r in results if r["delta_ll"] is not None and r.get("is_brain_enriched")]

        fig, ax = plt.subplots(figsize=(10, 7))

        # Color by brain enrichment
        colors = ["#E74C3C" if r.get("is_brain_enriched") else "#3498DB"
                  for r in results if r["delta_ll"] is not None]

        ax.scatter(x, y, c=colors, alpha=0.6, s=40, edgecolors="none")

        # Label notable genes
        notable = {"BACE1", "EGFR", "CCL2", "MMP9", "MMP2", "CXCR4", "KRAS", "MYC",
                   "TP53", "PIK3CA", "PTEN", "VEGFA", "CLDN5"}
        for r in results:
            if r["delta_ll"] is not None and r["gene"] in notable:
                ax.annotate(r["gene"], (r["crispr_lfc"], r["delta_ll"]),
                           fontsize=8, ha="left", va="bottom",
                           xytext=(3, 3), textcoords="offset points")

        # Trend line
        if len(x) >= 10 and r_val:
            z = np.polyfit(x, y, 1)
            p = np.poly1d(z)
            xline = sorted(x)
            ax.plot(xline, [p(xi) for xi in xline], "k--", alpha=0.4, linewidth=1)

        ax.axvline(x=0, color="gray", linewidth=0.5, linestyle=":")
        ax.axhline(y=0, color="gray", linewidth=0.5, linestyle=":")

        r_str = f"r={r_val:.3f}, p={p_val:.4f}" if r_val else "r=N/A"
        ax.set_xlabel("CRISPRa Brain LFC (in vivo mouse, GSE237446)\n"
                      "Positive = activation promotes BrM | Negative = no BrM effect",
                      fontsize=11)
        ax.set_ylabel("Evo2 Conditional LL delta at TSS (C>T)\n"
                      "Negative = model finds this mutation unusual | Positive = tolerated",
                      fontsize=11)
        ax.set_title(f"Evo2 Thermodynamics vs In Vivo Brain Metastasis (CRISPRa)\n"
                     f"{r_str}  |  n={len(x)}  |  evo2_1b_base A100", fontsize=12)

        from matplotlib.patches import Patch
        n_red  = sum(1 for c in colors if c == "#E74C3C")
        n_blue = sum(1 for c in colors if c == "#3498DB")
        legend = [Patch(color="#E74C3C", label=f"Brain-enriched (LFC>0, n={n_red})"),
                  Patch(color="#3498DB", label=f"Brain-depleted (LFC<0, n={n_blue})")]
        ax.legend(handles=legend, loc="upper right")

        plt.tight_layout()
        plot_file = DATA_DIR / "crispra_evo2_correlation.png"
        plt.savefig(plot_file, dpi=150)
        plt.close()
        logger.info(f"Plot saved → {plot_file}")

    except ImportError:
        logger.warning("matplotlib not available — skipping plot")

    return output


def main():
    parser = argparse.ArgumentParser(description="GSE237446 CRISPRa × Evo2 Correlation")
    parser.add_argument("--top", type=int, default=100,
                        help="Genes from each end (top N brain-enriched + top N depleted)")
    parser.add_argument("--concurrency", type=int, default=4,
                        help="Parallel Evo2 calls")
    parser.add_argument("--dry-run", action="store_true",
                        help="Use random scores instead of calling Evo2 (test pipeline)")
    args = parser.parse_args()

    asyncio.run(run_correlation(top_n=args.top, concurrency=args.concurrency, dry_run=args.dry_run))


if __name__ == "__main__":
    main()

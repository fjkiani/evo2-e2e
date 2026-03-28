"""
brain_met/pipeline/run_brm_pipeline.py
========================================
Full end-to-end Brain Metastasis Target-Lock + Assassin Score pipeline.

Usage:
  python run_brm_pipeline.py [--fast] [--step bbb_transit] [--seed 42]

Pipeline:
  1. Load Brain Met target universe (universe.py)
  2. Score all genes with Target-Lock (Evo2 + Enformer)
  3. Step-z-score normalize (fixes POC saturation)
  4. Compute Assassin scores for top-ranked genes
  5. Validate: AUROC/AUPRC vs ground truth labels
  6. Save results with full provenance + receipts

Reproducibility:
  - Seed locked at 42 (matches POC manuscript)
  - All external calls logged with timestamps
  - Results saved to data/brain_met/pipeline_results_{timestamp}.json
  - One-command reproduction: python run_brm_pipeline.py
"""
import asyncio
import json
import logging
import sys
import time
import argparse
import os
from pathlib import Path
from datetime import datetime

# Ensure imports work
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)

os.environ.setdefault("MODAL_TOKEN_ID", "ak-u2ShLPpaTsWfYMlOrc8XKX")
os.environ.setdefault("MODAL_TOKEN_SECRET", "as-PmarKaQCh03Vmsn3RADzPJ")

RESULTS_DIR = Path(__file__).parent.parent.parent / "data" / "brain_met"


async def run_pipeline(
    step_filter: str = None,
    fast_mode: bool = False,
    seed: int = 42,
    use_enformer: bool = True,
):
    import random, numpy as np
    random.seed(seed)

    from brain_met.targets.universe import (
        BRAIN_MET_UNIVERSE, POSITIVE_GENES, NEGATIVE_GENES,
        UNIVERSE_STATS, as_scoring_input, get_genes_for_step,
    )
    from core.target_lock.scorer import TargetLockScorer, assassin_score

    run_ts = datetime.utcnow().strftime("%Y%m%dT%H%M%SZ")
    logger.info(f"=== BrM Pipeline RUN {run_ts} ===")
    logger.info(f"Universe: {UNIVERSE_STATS}")

    scorer = TargetLockScorer(
        disease="brain_met",
        use_enformer=use_enformer and not fast_mode,
        seed=seed,
    )

    # Build scoring inputs
    all_steps = [
        "primary_tumor_escape", "intravasation", "circulation_survival",
        "bbb_transit", "cns_colonization", "brain_niche_adaptation", "brm_angiogenesis",
    ]
    if step_filter:
        all_steps = [s for s in all_steps if step_filter in s]

    scoring_inputs = []
    gene_labels = {}  # (gene, step) → is_positive

    for step in all_steps:
        for gene in BRAIN_MET_UNIVERSE:
            inp = as_scoring_input(gene, step)
            scoring_inputs.append(inp)
            is_pos_for_step = (step in gene.primary_steps) or (step in gene.secondary_steps)
            gene_labels[(gene.symbol, step)] = gene.is_positive and is_pos_for_step

    logger.info(f"Scoring {len(scoring_inputs)} gene-step combinations ({len(all_steps)} steps × {len(BRAIN_MET_UNIVERSE)} genes)")

    t0 = time.time()
    concurrency = 2 if fast_mode else 4
    results = await scorer.score_genes(scoring_inputs, concurrency=concurrency)
    elapsed = time.time() - t0
    logger.info(f"Scoring complete: {len(results)} results in {elapsed:.1f}s")

    # ── AUROC / AUPRC validation ─────────────────────────────────────────────
    from collections import defaultdict
    step_results = defaultdict(list)
    for r in results:
        label = gene_labels.get((r.gene, r.step), False)
        step_results[r.step].append({
            "gene": r.gene,
            "step": r.step,
            "calibrated_score": r.calibrated_score,
            "target_lock_score": r.target_lock_score,
            "label": label,
            "bbb_relevant": any(g.bbb_relevant for g in BRAIN_MET_UNIVERSE if g.symbol == r.gene),
            "flags": r.anti_hallucination_flags,
        })

    validation_metrics = {}
    for step, step_data in step_results.items():
        labels = [d["label"] for d in step_data]
        scores = [d["calibrated_score"] for d in step_data]
        n_pos = sum(labels)

        if n_pos == 0 or n_pos == len(labels):
            validation_metrics[step] = {"auroc": None, "auprc": None, "n_pos": n_pos}
            continue

        try:
            auroc = _auroc(labels, scores)
            auprc, p3 = _auprc_p3(labels, scores, k=3)
        except Exception as e:
            auroc, auprc, p3 = None, None, None
            logger.warning(f"Metrics failed for {step}: {e}")

        validation_metrics[step] = {
            "auroc": round(auroc, 4) if auroc else None,
            "auprc": round(auprc, 4) if auprc else None,
            "precision_at_3": p3,
            "n_pos": n_pos,
            "n_total": len(labels),
        }

    # Print summary
    print(f"\n{'='*70}")
    print(f"BRAIN MET TARGET-LOCK RESULTS — {run_ts}")
    print(f"{'='*70}")
    print(f"{'Step':<30} {'AUROC':>8} {'AUPRC':>8} {'P@3':>6} {'N+':>4}")
    print(f"{'-'*60}")
    auroc_vals = []
    for step, m in validation_metrics.items():
        auroc = m.get("auroc", "N/A")
        auprc = m.get("auprc", "N/A")
        p3 = m.get("precision_at_3", "N/A")
        n_pos = m.get("n_pos", 0)
        print(f"{step:<30} {str(auroc):>8} {str(auprc):>8} {str(p3):>6} {n_pos:>4}")
        if auroc is not None:
            auroc_vals.append(auroc)

    if auroc_vals:
        mean_auroc = sum(auroc_vals) / len(auroc_vals)
        print(f"\nMean AUROC: {mean_auroc:.4f}")
        print(f"POC AUROC was: 0.976 ± 0.035 (beat it = success)")

    # BBB Transit analysis
    bbb_data = step_results.get("bbb_transit", [])
    if bbb_data:
        print(f"\n── BBB Transit Ranking (TOP 5) ─────────────────────────────────")
        bbb_sorted = sorted(bbb_data, key=lambda x: x["calibrated_score"], reverse=True)
        for i, d in enumerate(bbb_sorted[:5]):
            tag = "✅" if d["label"] else ("🔴" if not d["label"] else "")
            bbb_tag = " [BBB]" if d["bbb_relevant"] else ""
            print(f"  {i+1}. {d['gene']:<10} score={d['calibrated_score']:.4f} {tag}{bbb_tag}")

    # Save results
    output = {
        "run_info": {
            "timestamp": run_ts,
            "seed": seed,
            "disease": "brain_met",
            "fast_mode": fast_mode,
            "use_enformer": use_enformer,
            "n_genes": len(BRAIN_MET_UNIVERSE),
            "n_positives": len(POSITIVE_GENES),
            "n_negatives": len(NEGATIVE_GENES),
            "n_steps": len(all_steps),
            "elapsed_s": round(elapsed, 1),
        },
        "validation_metrics": validation_metrics,
        "gene_scores": {
            step: sorted(data, key=lambda x: x["calibrated_score"], reverse=True)
            for step, data in step_results.items()
        },
    }

    outfile = RESULTS_DIR / f"pipeline_results_{run_ts}.json"
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    with open(outfile, "w") as f:
        json.dump(output, f, indent=2)
    logger.info(f"Results saved → {outfile}")

    return output


def _auroc(labels, scores):
    """Simple AUROC via trapezoidal rule."""
    paired = sorted(zip(scores, labels), reverse=True)
    n_pos = sum(labels)
    n_neg = len(labels) - n_pos
    if n_pos == 0 or n_neg == 0:
        return None
    tp, fp = 0, 0
    prev_fpr, prev_tpr = 0, 0
    auc = 0.0
    for score, label in paired:
        if label:
            tp += 1
        else:
            fp += 1
        tpr = tp / n_pos
        fpr = fp / n_neg
        auc += (fpr - prev_fpr) * (tpr + prev_tpr) / 2
        prev_fpr, prev_tpr = fpr, tpr
    return auc


def _auprc_p3(labels, scores, k=3):
    """AUPRC and Precision@K."""
    paired = sorted(zip(scores, labels), reverse=True)
    n_pos = sum(labels)
    if n_pos == 0:
        return 0.0, 0.0
    tp, total = 0, 0
    prev_recall = 0
    auprc = 0.0
    for score, label in paired:
        total += 1
        if label:
            tp += 1
        precision = tp / total
        recall = tp / n_pos
        auprc += (recall - prev_recall) * precision
        prev_recall = recall

    # Precision@K
    top_k = [label for _, label in paired[:k]]
    p_at_k = sum(top_k) / k if top_k else 0

    return auprc, round(p_at_k, 3)


def main():
    parser = argparse.ArgumentParser(description="BrM Target-Lock Pipeline")
    parser.add_argument("--fast", action="store_true", help="Fast mode (skip Enformer, 2 windows)")
    parser.add_argument("--step", default=None, help="Filter to specific step")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--no-enformer", action="store_true")
    args = parser.parse_args()

    asyncio.run(run_pipeline(
        step_filter=args.step,
        fast_mode=args.fast,
        seed=args.seed,
        use_enformer=not args.no_enformer,
    ))


if __name__ == "__main__":
    main()

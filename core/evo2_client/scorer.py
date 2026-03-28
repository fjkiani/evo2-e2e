"""
crispro_evo2_scorer.py — Production Evo2 scorer wired to our Modal service.

AUDIT vs crispro-backend-v2/api/services/sequence_scorers/evo2_scorer.py:
  ✅ Kept: Multi-window ensemble logic (flanks, best-delta selection)
  ✅ Kept: Forward/reverse symmetry averaging
  ✅ Kept: SeqScore dataclass, percentile mapping, curated PGx fallbacks
  ✅ Kept: BRCA1/2 truncating variant 0.80 floor
  🔧 Fixed: No longer calls dead HTTP endpoints
  🔧 Fixed: Calls OUR Modal service score_variant_conditional() instead
  🔧 Fixed: Uses conditional_ll (not mean_ll) — needs recalibrated percentile mapping
  🔧 Fixed: percentile_like() was calibrated for mean_ll values (0.001-0.1 range)
             Recalibrated for conditional_ll values (-2.0 to 0.5 range)

Key insight from benchmark audit:
  - evo2_7b mean_ll at 1024bp gives min_delta -0.019 → AUROC 0.971
  - evo2_1b conditional_ll at 8192bp gives delta -0.6 to +0.4 for BrM variants
  - We need to use conditional_ll delta magnitude for the percentile mapping
"""
import asyncio
import logging
import urllib.request
import json
import os
import time
from dataclasses import dataclass, field
from typing import Dict, Any, List, Optional, Tuple

logger = logging.getLogger(__name__)

MODAL_TOKEN_ID     = "ak-u2ShLPpaTsWfYMlOrc8XKX"
MODAL_TOKEN_SECRET = "as-PmarKaQCh03Vmsn3RADzPJ"
MODAL_APP_NAME     = "crispro-evo2-v9"
MODAL_CLS_NAME     = "Evo2Service"


@dataclass
class SeqScore:
    """Canonical sequence score result — matches crispro-backend-v2 schema."""
    variant: Dict[str, Any]
    sequence_disruption: float          # abs(delta_ll) — raw magnitude
    delta_ll: Optional[float] = None    # raw conditional LL delta (negative = deleterious)
    min_delta: Optional[float] = None   # best delta across windows
    exon_delta: Optional[float] = None  # exon-context delta (if computed)
    calibrated_seq_percentile: float = 0.0  # percentile [0,1]
    impact_level: str = "no_impact"
    scoring_mode: str = "unknown"
    best_model: Optional[str] = None
    best_window_bp: Optional[int] = None
    scoring_strategy: Dict[str, Any] = field(default_factory=dict)
    forward_reverse_meta: Optional[Dict[str, Any]] = None
    anti_hallucination_flags: List[str] = field(default_factory=list)


def percentile_like_conditional_ll(delta: float) -> float:
    """
    Map conditional_ll delta to percentile [0,1].
    Calibrated from our 12 BrM variant run (2026-03-28):
      delta < -0.5  → top 10%  (pathogenic)
      delta -0.5 to -0.2 → 60-80%
      delta -0.2 to -0.05 → 40-60%
      delta near 0 → 20-40%
      delta > 0 → bottom 20% (tolerated/neutral)

    Note: This will be recalibrated when we run the full ClinVar 160-variant benchmark.
    """
    d = abs(delta) if delta < 0 else 0.0  # only negative deltas matter for pathogenicity
    if delta > 0.3:     return 0.05   # model actually favors alt — likely benign
    if delta > 0.1:     return 0.15
    if delta > 0.0:     return 0.25
    if d < 0.05:        return 0.30
    if d < 0.10:        return 0.45
    if d < 0.20:        return 0.60
    if d < 0.35:        return 0.72
    if d < 0.50:        return 0.82
    if d < 0.70:        return 0.90
    return 0.95


def classify_impact_level(delta: float) -> str:
    """Classify based on conditional LL delta (negative = deleterious)."""
    if delta < -1.5:  return "catastrophic"
    if delta < -0.8:  return "high"
    if delta < -0.3:  return "moderate"
    if delta < -0.05: return "low"
    return "no_impact"


COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")

def complement(b: str) -> str:
    return b.translate(COMPLEMENT)


async def fetch_hg38_sequence(chrom: str, start_0: int, end_0: int) -> Optional[str]:
    """Fetch hg38 sequence via UCSC REST API (0-based coords)."""
    try:
        url = (f"https://api.genome.ucsc.edu/getData/sequence"
               f"?genome=hg38&chrom={chrom}&start={start_0}&end={end_0}")
        req = urllib.request.Request(url, headers={"User-Agent": "CrisPRO/1.0"})
        loop = asyncio.get_event_loop()
        def _fetch():
            with urllib.request.urlopen(req, timeout=30) as r:
                return json.loads(r.read())["dna"].upper()
        return await loop.run_in_executor(None, _fetch)
    except Exception as e:
        logger.warning(f"UCSC fetch failed {chrom}:{start_0}-{end_0}: {e}")
        return None


async def score_single_window(
    chrom: str, pos: int, ref: str, alt: str, flank: int
) -> Optional[Tuple[float, int]]:
    """
    Score a single variant in a given window using conditional LL.
    Returns (delta_ll, context_len) or None on failure.
    """
    start = max(0, pos - flank - 1)
    end = pos + flank
    seq = await fetch_hg38_sequence(chrom, start, end)
    if not seq:
        return None

    var_idx = pos - 1 - start
    if var_idx < 0 or var_idx >= len(seq):
        return None

    actual = seq[var_idx]
    if actual == ref:
        eff_alt = alt
    elif actual == complement(ref):
        eff_alt = complement(alt)
    else:
        eff_alt = alt  # best effort on coord mismatch

    wt_seq  = seq
    mut_seq = seq[:var_idx] + eff_alt + seq[var_idx+1:]

    try:
        import modal
        os.environ.setdefault("MODAL_TOKEN_ID", MODAL_TOKEN_ID)
        os.environ.setdefault("MODAL_TOKEN_SECRET", MODAL_TOKEN_SECRET)
        Evo2Service = modal.Cls.from_name(MODAL_APP_NAME, MODAL_CLS_NAME)
        loop = asyncio.get_event_loop()
        result = await loop.run_in_executor(
            None,
            lambda: Evo2Service().score_variant_conditional.remote(wt_seq, mut_seq, var_idx)
        )
        if "error" in result:
            logger.warning(f"Evo2 error at flank={flank}: {result['error']}")
            return None
        return (result["delta_ll"], len(seq))
    except Exception as e:
        logger.warning(f"Evo2 call failed at flank={flank}: {e}")
        return None


async def score_variant(
    variant: Dict[str, Any],
    window_flanks: List[int] = None,
    symmetry: bool = True,
) -> SeqScore:
    """
    Score a single variant with multi-window ensemble + symmetry averaging.

    window_flanks: list of half-window sizes. Default [500, 2048, 4096, 8192]
    symmetry: if True, average forward + reverse (ref→alt and alt→ref)

    Anti-hallucination guarantees:
    - All curated fallbacks are labeled with scoring_mode='curated_fallback_*'
    - Real Evo2 calls labeled scoring_mode='evo2_conditional_ensemble'
    - If all Evo2 calls fail, returns scoring_mode='failed' with disruption=0
    """
    if window_flanks is None:
        window_flanks = [500, 2048, 4096, 8192]  # start small, go big

    chrom = str(variant.get("chrom", ""))
    pos   = int(variant.get("pos", 0))
    ref   = str(variant.get("ref", "")).upper()
    alt   = str(variant.get("alt", "")).upper()
    gene  = str(variant.get("gene", "")).upper()
    hgvs_p = str(variant.get("hgvs_p", "")).upper()
    consequence = str(variant.get("consequence", "")).lower()
    flags: List[str] = []

    # ── CURATED FALLBACKS (same as backend-v2, labeled honestly) ───────────────
    valid_bases = {"A", "C", "G", "T"}
    if ref not in valid_bases or alt not in valid_bases or not chrom or not pos:
        # Cannot call Evo2 — use curated prior
        disruption, impact = _curated_prior(gene, hgvs_p, consequence)
        return SeqScore(
            variant=variant,
            sequence_disruption=disruption,
            calibrated_seq_percentile=percentile_like_conditional_ll(-disruption),
            impact_level=impact,
            scoring_mode="curated_fallback_missing_alleles",
            anti_hallucination_flags=[f"Missing valid ref/alt: ref={ref} alt={alt}"],
        )

    # ── MULTI-WINDOW ENSEMBLE ──────────────────────────────────────────────────
    tasks = [score_single_window(chrom, pos, ref, alt, f) for f in window_flanks]
    if symmetry:
        rev_tasks = [score_single_window(chrom, pos, alt, ref, f) for f in window_flanks]
    else:
        rev_tasks = []

    all_results = await asyncio.gather(*tasks, *rev_tasks, return_exceptions=True)
    forward_results = all_results[:len(window_flanks)]
    reverse_results = all_results[len(window_flanks):]

    valid_deltas = []
    best_window = None
    windows_tested = []

    for i, (f_res, flank) in enumerate(zip(forward_results, window_flanks)):
        r_res = reverse_results[i] if symmetry and i < len(reverse_results) else None

        f_delta = None
        if isinstance(f_res, tuple):
            f_delta, ctx_len = f_res
        elif isinstance(f_res, Exception):
            flags.append(f"window flank={flank} failed: {f_res}")

        r_delta = None
        if symmetry and isinstance(r_res, tuple):
            r_delta, _ = r_res

        if f_delta is not None:
            # Average with reverse if available, else use forward only
            avg = (f_delta + (r_delta if r_delta is not None else 0.0)) / (2 if r_delta is not None else 1)
            valid_deltas.append((avg, flank))
            windows_tested.append({"flank": flank, "delta": avg, "fwd": f_delta, "rev": r_delta})

    if not valid_deltas:
        return SeqScore(
            variant=variant,
            sequence_disruption=0.0,
            calibrated_seq_percentile=0.0,
            impact_level="no_impact",
            scoring_mode="failed",
            anti_hallucination_flags=flags + ["All Evo2 calls failed"],
        )

    # Best delta = most negative (most deleterious)
    best_delta, best_flank = min(valid_deltas, key=lambda x: x[0])
    disruption = abs(best_delta)
    pct = percentile_like_conditional_ll(best_delta)

    # ── BRCA1/2 TRUNCATING FLOOR (from backend-v2, preserved) ─────────────────
    trunc = any(x in consequence for x in ["stop_gained", "frameshift", "nonsense", "splice"])
    trunc = trunc or any(x in hgvs_p for x in ["*", "FS", "TER"])
    if gene in {"BRCA1", "BRCA2"} and trunc:
        original_pct = pct
        pct = max(pct, 0.80)
        if pct != original_pct:
            flags.append(f"BRCA truncating floor applied: {original_pct:.2f}→{pct:.2f}")

    return SeqScore(
        variant=variant,
        sequence_disruption=disruption,
        delta_ll=best_delta,
        min_delta=best_delta,
        calibrated_seq_percentile=pct,
        impact_level=classify_impact_level(best_delta),
        scoring_mode="evo2_conditional_ensemble",
        best_model="evo2_1b_base",
        best_window_bp=best_flank * 2,
        scoring_strategy={
            "approach": "conditional_ll_multi_window",
            "windows_tested": windows_tested,
            "symmetry": symmetry,
            "modal_app": MODAL_APP_NAME,
        },
        forward_reverse_meta=windows_tested,
        anti_hallucination_flags=flags,
    )


def _curated_prior(gene: str, hgvs_p: str, consequence: str) -> Tuple[float, str]:
    """
    Curated fallback priors for variants where Evo2 cannot be called.
    Identical to backend-v2 logic, but labeled honestly.
    """
    # PGx complete deficiency
    if gene == "DPYD" and any(k in hgvs_p for k in ("*2A", "1905+1", "IVS14+1")):
        return 2.0, "critical"
    if gene == "DPYD" and any(k in hgvs_p for k in ("*13", "2846", "1679")):
        return 1.5, "high"
    if gene == "TPMT" and any(k in hgvs_p for k in ("*3A", "*3B", "*3C", "*2")):
        return 1.5, "high"
    # Strong LOF
    if any(k in consequence for k in ("frameshift", "stop_gained", "nonsense", "splice")):
        return 1.5, "high"
    # Known pathogenic missense hotspot
    if gene == "BRCA1" and "C61" in hgvs_p:
        return 1.2, "high"
    # Missense default prior
    if "missense" in consequence:
        return 0.3, "moderate"
    return 0.1, "low"


async def score_variants_batch(
    variants: List[Dict[str, Any]],
    window_flanks: List[int] = None,
    symmetry: bool = True,
    concurrency: int = 3,
) -> List[SeqScore]:
    """Score multiple variants with bounded concurrency."""
    semaphore = asyncio.Semaphore(concurrency)
    async def _score_one(v):
        async with semaphore:
            return await score_variant(v, window_flanks, symmetry)
    return await asyncio.gather(*[_score_one(v) for v in variants])

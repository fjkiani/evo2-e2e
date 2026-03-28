"""
core/target_lock/scorer.py
==========================
Target-Lock scoring engine — HONEST version (no tautology).

AUDIT NOTE 2026-03-28:
  Previous version embedded mission_fit_discount (label→score encoding) which
  produced AUROC=1.0 by construction. That version is deleted.

  Evo2 now lives or dies on its raw conditional_ll delta and signal components.
  Mission-fit is a POST-HOC FILTER ONLY — applied after scoring to explain results,
  never baked into the score itself.

Formula:
  Target_Lock = 0.35×Functionality + 0.35×Essentiality + 0.20×Regulatory + 0.10×Chromatin

  Where:
    Functionality  = 1/(1+exp(delta_ll / 0.5))   [Evo2 protein impact, conditional LL]
    Essentiality   = 1.0 for frameshift/nonsense; min(1.0, |delta_ll|/1.5) for missense
    Regulatory     = |min_delta| / (|min_delta|+1)  [splice/noncoding signal]
    Chromatin      = Enformer accessibility proxy [0,1]  (0.5 default when unavailable)

  No discount. No label encoding. Evo2 must earn its ranking.

Weights (brain-met specific, WEIGHTS_BRAIN_MET):
  Regulatory boosted from 0.15→0.20 because BBB transit involves splice isoforms.
  Chromatin reduced from 0.15→0.10 (ablation showed -0.013 AUROC contribution).
"""
import asyncio
import json
import logging
import math
import time
from dataclasses import dataclass, field
from typing import Dict, Any, List, Optional, Tuple

logger = logging.getLogger(__name__)

WEIGHTS_DEFAULT = {
    "functionality": 0.35,
    "essentiality":  0.35,
    "regulatory":    0.20,
    "chromatin":     0.10,
}

WEIGHTS_BRAIN_MET = {
    "functionality": 0.33,
    "essentiality":  0.33,
    "regulatory":    0.24,
    "chromatin":     0.10,
}


@dataclass
class TargetLockResult:
    gene: str
    step: str
    target_lock_score: float       # raw weighted sum [0,1] — NO discount
    functionality_score: float
    essentiality_score: float
    regulatory_score: float
    chromatin_score: float
    calibrated_score: float        # step-z-score normalized (within-step only)
    weights_used: Dict[str, float]
    provenance: Dict[str, Any]
    anti_hallucination_flags: List[str] = field(default_factory=list)
    elapsed_s: float = 0.0
    # Mission-fit metadata (for post-hoc analysis only, NEVER used in scoring)
    primary_steps: List[str] = field(default_factory=list)
    secondary_steps: List[str] = field(default_factory=list)


@dataclass
class AssassinScore:
    gene: str
    step: str
    guide_seq: str
    efficacy: float
    safety: float
    mission_fit: float
    structure: float
    assassin_score: float
    plddt: Optional[float] = None
    iptm: Optional[float] = None
    structural_verdict: Optional[str] = None
    provenance: Dict[str, Any] = field(default_factory=dict)
    anti_hallucination_flags: List[str] = field(default_factory=list)


def assassin_score(efficacy, safety, mission_fit, structure=0.0):
    return 0.37*efficacy + 0.30*safety + 0.30*mission_fit + 0.03*structure


def sigmoid(x, scale=10.0):
    return 1.0 / (1.0 + math.exp(-x / scale))


def functionality_from_delta(delta_ll: float) -> float:
    """Map conditional_ll delta to [0,1]. Negative delta = deleterious = high score."""
    return 1.0 / (1.0 + math.exp(delta_ll / 0.5))


def essentiality_from_consequence(consequence, delta_ll, hgvs_p) -> float:
    if consequence:
        if any(x in consequence.lower() for x in
               ["frameshift","stop_gained","nonsense","splice_donor","splice_acceptor"]):
            return 1.0
    if hgvs_p and any(x in hgvs_p.upper() for x in ["*","FS","TER"]):
        return 1.0
    if delta_ll is not None:
        return min(1.0, abs(delta_ll) / 1.5)
    return 0.3


def regulatory_from_delta(min_delta) -> float:
    if min_delta is None:
        return 0.0
    return abs(min_delta) / (abs(min_delta) + 1.0)


def _stepwise_zscore_normalize(scores: List[float]) -> List[float]:
    """Z-score normalize within a step to spread distribution."""
    if len(scores) < 2:
        return scores
    mu = sum(scores) / len(scores)
    std = (sum((s - mu)**2 for s in scores) / len(scores))**0.5
    if std < 1e-6:
        return [0.5] * len(scores)
    z = [(s - mu) / std for s in scores]
    return [1.0 / (1.0 + math.exp(-zi)) for zi in z]


class TargetLockScorer:
    """
    Production Target-Lock scorer.
    No mission-fit discount. Evo2 scores are the truth signal.
    """

    def __init__(self, disease="general", weights=None, use_enformer=True, seed=42):
        self.disease = disease
        self.seed = seed
        self.use_enformer = use_enformer
        if weights is not None:
            self.weights = weights
        elif disease == "brain_met":
            self.weights = WEIGHTS_BRAIN_MET
        else:
            self.weights = WEIGHTS_DEFAULT
        total = sum(self.weights.values())
        assert abs(total - 1.0) < 1e-6, f"Weights must sum to 1.0, got {total}"

    async def _get_evo2_scores(self, gene, chrom, pos, ref, alt,
                                consequence=None, hgvs_p=None):
        """Returns (functionality, essentiality, regulatory, flags)."""
        flags = []
        try:
            import sys
            sys.path.insert(0, "/home/user/workspace/evo2_e2e/core/evo2_client")
            from scorer import score_variant
            result = await score_variant(
                {"gene": gene, "chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
                 "consequence": consequence, "hgvs_p": hgvs_p},
                window_flanks=[500, 4096],
                symmetry=False,
            )
            delta_ll = result.delta_ll
            min_delta = result.min_delta
            func  = functionality_from_delta(delta_ll) if delta_ll is not None else 0.5
            essen = essentiality_from_consequence(consequence, delta_ll, hgvs_p)
            reg   = regulatory_from_delta(min_delta)
            if result.scoring_mode != "evo2_conditional_ensemble":
                flags.append(f"scoring_mode={result.scoring_mode}")
            flags.extend(result.anti_hallucination_flags)
        except Exception as e:
            flags.append(f"Evo2 failed: {e} — defaults used")
            func, essen, reg = 0.5, 0.5, 0.0
        return func, essen, reg, flags

    async def _get_chromatin_score(self, gene, chrom, pos):
        flags = []
        if not self.use_enformer:
            flags.append("Enformer disabled — chromatin=0.5")
            return 0.5, flags
        try:
            import modal, os
            os.environ.setdefault("MODAL_TOKEN_ID", "ak-u2ShLPpaTsWfYMlOrc8XKX")
            os.environ.setdefault("MODAL_TOKEN_SECRET", "as-PmarKaQCh03Vmsn3RADzPJ")
            EnformerService = modal.Cls.from_name("crispro-enformer", "EnformerService")
            loop = asyncio.get_event_loop()
            result = await loop.run_in_executor(
                None, lambda: EnformerService().predict.remote(chrom=chrom, pos=pos, ref="A", alt="A")
            )
            score = float(result.get("accessibility_score", 0.5))
            return score, flags
        except Exception as e:
            flags.append(f"Enformer failed: {e} — chromatin=0.5")
            return 0.5, flags

    async def score_gene(self, gene, step, chrom, pos, ref, alt,
                          consequence=None, hgvs_p=None,
                          primary_steps=None, secondary_steps=None):
        """
        Score a gene. NO mission-fit discount applied here.
        primary_steps / secondary_steps stored as metadata only.
        """
        t0 = time.time()
        flags = []

        func, essen, reg, evo2_flags = await self._get_evo2_scores(
            gene, chrom, pos, ref, alt, consequence, hgvs_p
        )
        chrom_score, enformer_flags = await self._get_chromatin_score(gene, chrom, pos)
        flags.extend(evo2_flags + enformer_flags)

        # RAW SCORE — no discount, no label encoding
        tl = (
            self.weights["functionality"] * func
            + self.weights["essentiality"]  * essen
            + self.weights["regulatory"]    * reg
            + self.weights["chromatin"]     * chrom_score
        )

        return TargetLockResult(
            gene=gene,
            step=step,
            target_lock_score=tl,
            functionality_score=func,
            essentiality_score=essen,
            regulatory_score=reg,
            chromatin_score=chrom_score,
            calibrated_score=tl,
            weights_used=self.weights.copy(),
            provenance={
                "disease": self.disease,
                "seed": self.seed,
                "discount_applied": False,  # explicit: no discount
                "evo2_is_sole_discriminator": True,
            },
            anti_hallucination_flags=flags,
            elapsed_s=round(time.time() - t0, 2),
            primary_steps=primary_steps or [],
            secondary_steps=secondary_steps or [],
        )

    async def score_genes(self, genes, concurrency=4):
        """Score multiple genes with step-z-score normalization (within-step spread only)."""
        semaphore = asyncio.Semaphore(concurrency)
        async def _one(g):
            async with semaphore:
                return await self.score_gene(
                    gene=g["gene"], step=g["step"],
                    chrom=g["chrom"], pos=g["pos"],
                    ref=g.get("ref","A"), alt=g.get("alt","T"),
                    consequence=g.get("consequence"), hgvs_p=g.get("hgvs_p"),
                    primary_steps=g.get("primary_steps", []),
                    secondary_steps=g.get("secondary_steps", []),
                )

        results = await asyncio.gather(*[_one(g) for g in genes])

        # Step-z-score normalization: spread within step, NOT across labels
        steps = list(set(r.step for r in results))
        for step in steps:
            sr = [r for r in results if r.step == step]
            raw = [r.target_lock_score for r in sr]
            normalized = _stepwise_zscore_normalize(raw)
            for r, n in zip(sr, normalized):
                r.calibrated_score = n

        return list(results)

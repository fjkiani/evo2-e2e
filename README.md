# evo2-e2e: Production Metastasis Interception Pipeline

**Brain Metastasis → Breast Metastasis → Any Cancer**

Production implementation of the Target-Lock / Assassin Score framework from the
*Intercepting metastasis: 8-step CRISPR design via multi-modal foundation models* manuscript.

## Architecture

```
evo2-e2e/
├── core/                      # Modality-agnostic scoring engines
│   ├── target_lock/           # Target-Lock scorer (Evo2 + Enformer)
│   ├── assassin_score/        # Assassin composite ranking
│   ├── evo2_client/           # Modal Evo2 service (crispro-evo2-v9, A100)
│   ├── enformer_client/       # Modal Enformer service
│   └── acmg/                  # ACMG/AMP VUS classifier (anti-hallucination)
│
├── brain_met/                 # Brain metastasis module
│   ├── targets/universe.py    # 19 positives + 9 hard negatives, 7 BrM steps
│   ├── guides/                # CRISPR guide design for BrM targets
│   ├── structural/            # AlphaFold3 structural validation
│   └── pipeline/run_brm_pipeline.py  # End-to-end pipeline
│
├── breast_met/                # Breast metastasis (extend universe.py pattern)
│   ├── targets/               # AURORA + TCGA-BRCA gene universe
│   └── guides/
│
├── shared/
│   ├── gene_coords/           # hg38 TSS coordinates (Ensembl validated)
│   ├── ground_truth/          # metastasis_rules_v1.0.1.json
│   └── calibration/           # Gene-specific Evo2 calibration snapshots
│
├── tests/
│   ├── stress/test_antihallucination.py  # 17-test suite, ClinVar ground truth
│   ├── integration/           # End-to-end pipeline tests
│   └── unit/                  # Component tests
│
├── scripts/
│   └── reproduce_all.sh       # One-command reproduction
│
└── data/
    ├── brain_met/             # Real data: GSE237446, GSE205033, MSK-MET, AURORA
    └── validation/            # live_variant_scores.json (Evo2 A100 scores)
```

## Quick Start

```bash
# Install dependencies
pip install -r requirements.txt

# Run full BrM pipeline (fast mode: no Enformer, 2 Evo2 windows)
python brain_met/pipeline/run_brm_pipeline.py --fast

# Run full pipeline with Enformer
python brain_met/pipeline/run_brm_pipeline.py

# Run stress tests (anti-hallucination)
python tests/stress/test_antihallucination.py --fast

# Score specific step
python brain_met/pipeline/run_brm_pipeline.py --step bbb_transit
```

## Deployed Services (Modal, crispro-test workspace)

| Service | App | GPU | Status |
|---------|-----|-----|--------|
| Evo2 1B base | `crispro-evo2-v9` | A100 | ✅ LIVE |
| Enformer | `crispro-enformer` | T4 | ✅ LIVE |
| BrM ICL adapter | `brm-icl-v1` | CPU | ✅ LIVE |

## Live Results (2026-03-28)

```
BrM Clinical Variants — Evo2 Conditional LL (A100, hg38, 8192bp):
  PIK3CA p.H1047R   delta_ll = -0.615  (top penalized — PI3K kinase domain)
  TP53   p.R175H    delta_ll = -0.418  (BrM hotspot, 2× enriched)
  ESR1   p.D538G    delta_ll = -0.402  (ligand-independent ER activation)
  BACE1  p.D289N    delta_ll = +0.002  (near neutral — ACMG PVS1 context needed)
```

## What Was Fixed vs the POC

| POC Issue | Fix |
|-----------|-----|
| Target-Lock saturation (0.352–0.355 band) | Step-z-score normalization in `scorer.py` |
| Chromatin 15% → −0.013 AUROC | Downweighted to 10%; regulatory increased to 20% |
| Off-target: no bulge/CFD weighting | `anti_hallucination_flags` documents gap; CFD TODO |
| 15 guides only in structural validation | Expanded framework: top-5/step reproducible script |
| IDH1/IDH2 used as glioma negatives | Hard negatives added: IDH1, IDH2, FLT3, ABL1, etc. |
| PM2 always fires without gnomAD | Fixed in ACMG classifier (requires real gnomAD data) |
| All Evo2 paper tests failed | New `score_variant_conditional` endpoint (conditional LL) |

## Scoring Methods

### Target-Lock Formula
```
Target_Lock = 0.35×Functionality + 0.35×Essentiality + 0.20×Regulatory + 0.10×Chromatin

Brain-Met variant (WEIGHTS_BRAIN_MET):
  0.33×Functionality + 0.33×Essentiality + 0.24×Regulatory + 0.10×Chromatin
  (Higher regulatory: BBB transit involves splice isoforms)
```

### Assassin Score Formula
```
Assassin = 0.37×Efficacy + 0.30×Safety + 0.30×Mission + 0.03×Structure
```

### Scoring Method
- **Functionality**: `1 / (1 + exp(delta_ll / 0.5))` — conditional LL from Evo2
- **Essentiality**: `1.0` for frameshift/nonsense; `min(1.0, abs(delta_ll)/1.5)` for missense
- **Regulatory**: `|min_delta| / (|min_delta| + 1)` — splice/noncoding impact
- **Chromatin**: Enformer accessibility proxy [0,1]

## Brain Met Target Universe

**7 BrM-specific steps** (vs 8 general metastasis steps in POC):
1. `primary_tumor_escape` — escape from breast/lung/melanoma primary
2. `intravasation` — blood/lymph entry
3. `circulation_survival` — CTC anoikis resistance
4. `bbb_transit` ← **NEW** — blood-brain barrier crossing
5. `cns_colonization` — brain micrometastasis establishment
6. `brain_niche_adaptation` — astrocyte/neuron microenvironment adaptation
7. `brm_angiogenesis` — brain-specific angiogenic switch

**19 positive targets** (evidence from AURORA, MSK-MET, GSE237446, GSE205033, literature)
**9 hard negatives** (IDH1, IDH2 = glioma/NOT BrM; ABL1, FLT3 = hematologic; etc.)

## Reproducibility

- Seed: 42 (matches manuscript)
- Model: `crispro-evo2-v9` (A100, evo2_1b_base, conditional LL scoring)
- All variants scored against real hg38 sequences from UCSC API
- Results timestamped and saved to `data/brain_met/pipeline_results_{ts}.json`
- Anti-hallucination flags in every result object

## Research Use Only

This framework is for research purposes only and has not been validated for clinical use.
All predictions require experimental validation before therapeutic application.

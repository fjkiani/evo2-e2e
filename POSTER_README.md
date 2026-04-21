# CrisPRO BrM Poster — Hunt. Lock. Strike.

> **Computational Prioritization of Brain Metastasis Therapeutic Targets**
> Fahad Kiani | CrisPRO.ai | April 2026

---

## Narrative Arc

The poster tells a 3-act story of computational target discovery for brain metastasis (BrM):

### ACT I: HUNT
*How do we find the right targets in a sea of candidates?*

The CrisPRO pipeline integrates three orthogonal data streams — CRISPRa functional screens (GSE237446), chromatin accessibility (ATAC-seq, GSE205033), and Evo2 evolutionary language model scores — into a single **Assassin Score** that ranks 29 candidate genes by their BrM-specific functional importance. This act establishes the mission, the engine, and the scoring formula.

### ACT II: LOCK
*Which targets survive multi-modal scrutiny?*

The top 10 candidates are stress-tested against 6 independent evidence streams: CRISPRa LFC, ATAC accessibility, Evo2 log-likelihood, RNA differential expression (external), CPTAC proteomics (external), and DepMap essentiality. Two targets — **NOM1** and **BACE1** — emerge as CONFIRMED with convergent evidence. The Evo2 deep-dive validates the model's discriminative power (AUROC=0.72 on 29-gene panel).

### ACT III: STRIKE
*What do we do with the locked targets?*

BrM biology is validated against the largest available cohort (MSK-MET 2021, N=25,775), confirming that known BrM drivers (EGFR, STK11, KEAP1) show expected enrichment patterns. The two confirmed targets are profiled for druggability: NOM1 as a first-in-class nucleolar GTPase opportunity, BACE1 as a repurposing candidate (verubecestat/MK-8931). FAM72A is flagged as a combination candidate via the Aurora kinase pathway.

---

## Panel Guide

| Panel | Title | Act | Data Source | Verified? |
|-------|-------|-----|-------------|-----------|
| A | Mission Fit | ACT I: HUNT | GSE237446 (CRISPRa), GSE205033 (ATAC) | Yes — primary data |
| B | Pipeline Engine | ACT I: HUNT | Pipeline architecture diagram | Yes — conceptual |
| C | Assassin Score Breakdown | ACT I: HUNT | class_aware_composite_final.json | Yes — verified formula |
| D | Evo2 Deep-Dive | ACT II: LOCK | honest_auroc_panel.json (A100 inference) | Yes — AUROC=0.72 verified |
| E | Multi-Modal Evidence Matrix | ACT II: LOCK | 6 streams (see Data Provenance) | Partial — see notes |
| F | BrM Biology Validated | ACT III: STRIKE | step2_validation_results.json (MSK-MET 2021) | Yes — cBioPortal data |
| G | Target Lock | ACT III: STRIKE | Composite from all verified sources | Yes — derived values |
| H | Strike Package | ACT III: STRIKE | Literature + DepMap + CRISPRa | Yes — cited sources |

---

## Asset Inventory

### Poster Files
| File | Location | Size | Description |
|------|----------|------|-------------|
| crispro_brm_poster_v2.png | /mnt/results/ | ~6800x5200px | Full poster, 150 DPI |
| crispro_brm_poster_v2_preview.png | /mnt/results/ | 2400x1800px | Preview version |
| crispro_brm_poster_v2.pdf | /tmp/results-staging/ | — | PDF version |

### Panel PNGs
| File | Location |
|------|----------|
| panel_A_mission_fit.png | /mnt/results/poster_panels_v2/ |
| panel_B_engine.png | /mnt/results/poster_panels_v2/ |
| panel_C_assassin_breakdown.png | /mnt/results/poster_panels_v2/ |
| panel_D_evo2_deepdive.png | /mnt/results/poster_panels_v2/ |
| panel_E_evidence_matrix.png | /mnt/results/poster_panels_v2/ |
| panel_F_brm_validated.png | /mnt/results/poster_panels_v2/ |
| panel_G_target_lock.png | /mnt/results/poster_panels_v2/ |
| panel_H_strike_package.png | /mnt/results/poster_panels_v2/ |

### Data Files
| File | Description |
|------|-------------|
| class_aware_composite_final.json | Assassin Score formula, per-gene components (BACE1 verified) |
| honest_auroc_panel.json | Evo2 per-gene scores, 29-gene panel, AUROC=0.72 |
| step2_validation_results.json | cBioPortal MSK-MET 2021 validation results |
| adversarial_negatives_scored.json | Adversarial attack analysis results |

### Supporting Figures
| File | Description |
|------|-------------|
| honest_auroc_panel.png | Evo2 AUROC panel figure |
| attack_analysis.png | Adversarial attack robustness figure |
| class_aware_composite_auroc.png | Composite AUROC figure |
| step2_prospective_validation.png | BrM biology validation (original, light background) |

---

## Data Provenance

### VERIFIED (primary data on disk, reproducible)
| Data Stream | Source | File/Location | Notes |
|-------------|--------|---------------|-------|
| CRISPRa LFC | GSE237446 (Chafe 2025) | /mnt/datalake/crispr/ | Verified this session |
| ATAC-seq | GSE205033 | /mnt/datalake/atac/ | Verified this session |
| Evo2 LL scores | A100 GPU inference | honest_auroc_panel.json | AUROC=0.72 verified |
| DepMap Chronos (NOM1) | DepMap 23Q4 | /mnt/datalake/depmap/ | NOM1 = -0.736 verified |
| MSK-MET 2021 | cBioPortal | step2_validation_results.json | N=25,775 verified |
| Assassin Score formula | Composite | class_aware_composite_final.json | Formula verified |

### EXTERNAL ANALYSIS (not on disk, cited from literature/agent reports)
| Data Stream | Source | Status |
|-------------|--------|--------|
| RNA padj values | GSE223499 (external analysis) | Files archived, not on disk |
| CPTAC Protein levels | CPTAC BrM cohort (external) | Files archived, not on disk |
| DepMap (non-NOM1 genes) | DepMap 23Q4 | Not individually verified this session |

> **Important:** All "external analysis" data in Panel E is clearly labeled with asterisks (*). Only NOM1 DepMap Chronos = -0.736 was directly verified from /mnt/datalake/depmap/ during this session.

---

## Key Numbers

| Metric | Value | Source | Verified? |
|--------|-------|--------|-----------|
| Evo2 AUROC (29-gene panel) | 0.72 | honest_auroc_panel.json | Yes |
| Enformer AUROC (excluded) | 0.41 | Session analysis | Yes — excluded as non-discriminative |
| NOM1 DepMap Chronos | -0.736 | /mnt/datalake/depmap/ | Yes |
| BACE1 CRISPRa LFC | High (top quartile) | GSE237446 | Yes |
| NOM1 RNA padj | 0.003 | GSE223499 (external) | External analysis |
| NOM1 CPTAC protein | +1.81 | CPTAC BrM (external) | External analysis |
| EGFR BrM OR | 3.91 (p<1e-50) | MSK-MET 2021 | Yes |
| KRAS BrM OR | 0.78 (p<1e-8) | MSK-MET 2021 | Yes |
| MSK-MET cohort size | N=25,775 | cBioPortal | Yes |
| Top 10 targets somatic mutation rate | 0% | MSK-MET 2021 | Yes — by design |

---

## Reproduction Instructions

### Environment
```bash
# Python 3.9+
pip install matplotlib numpy pillow scipy
```

### Data Requirements
```
/mnt/datalake/crispr/     # GSE237446 CRISPRa data
/mnt/datalake/atac/       # GSE205033 ATAC-seq data
/mnt/datalake/depmap/     # DepMap 23Q4 Chronos scores
/mnt/results/             # Output directory
```

### Panel Generation Order
```bash
# Panels A-D were built in a prior session (do not rebuild)
# Panels E-H:

# Panel E: Multi-Modal Evidence Matrix
python poster/panels/build_panel_E.py

# Panel F: BrM Biology Validated
python poster/panels/build_panel_F.py

# Panel G: Target Lock
python poster/panels/build_panel_G.py

# Panel H: Strike Package
python poster/panels/build_panel_H.py

# Assemble full poster
python poster/assemble_poster.py
```

### Output Files
All panels save to `/mnt/results/poster_panels_v2/`
Full poster saves to `/mnt/results/crispro_brm_poster_v2.png`
PDF saves to `/tmp/results-staging/crispro_brm_poster_v2.pdf`

---

## Integrity Audit

### What is VERIFIED (reproducible from data on disk)
- [x] CRISPRa LFC rankings (GSE237446) — top genes confirmed
- [x] ATAC accessibility scores (GSE205033) — verified
- [x] Evo2 log-likelihood scores — A100 inference, AUROC=0.72
- [x] NOM1 DepMap Chronos = -0.736 — directly read from /mnt/datalake/depmap/
- [x] MSK-MET 2021 OR values — from step2_validation_results.json
- [x] Assassin Score formula — from class_aware_composite_final.json
- [x] BACE1 component breakdown — from class_aware_composite_final.json

### What is CLAIMED (external analysis, not independently verified this session)
- [ ] RNA padj values (GSE223499) — reported by RNA analysis agent, files archived
- [ ] CPTAC protein levels — reported by proteomics agent, files archived
- [ ] DepMap Chronos for genes other than NOM1 — not individually verified

### Enformer Exclusion Note
Enformer was evaluated as a 4th evidence stream but excluded from the final pipeline because its AUROC = 0.41 on the 29-gene BrM panel, which is **below random chance** (0.50). Including a non-discriminative model would degrade pipeline performance. This exclusion is documented and the AUROC value is verified.

### Data Integrity Principles Applied
1. No fabricated per-gene component breakdowns beyond what is in class_aware_composite_final.json
2. NOM1 DepMap value = -0.736 (verified), not -0.695 (prior agent claim)
3. All "external analysis" data clearly labeled with asterisks in Panel E
4. Enformer excluded with documented reason
5. Top 10 targets have 0% somatic mutation rate — this is a feature, not a bug (functional dependencies, not somatic drivers)

---

*CrisPRO.ai 2026 | For Research Use Only — not for clinical or diagnostic use*

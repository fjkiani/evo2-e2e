# MANIFEST_NOTES.md — Discrepancy Log

**Date:** 2026-03-30  
**Status:** CANONICAL — supersedes all prior narrative references to dataset provenance

---

## ATAC Discrepancy — RESOLVED

### What the conflict was

Two different GEO accessions were referenced as "the ATAC dataset" in evo2-e2e work:

| Reference | Accession | Where it appeared |
|---|---|---|
| Primary (canonical) | **GSE205033** | `brm_data_loader.py`, `brm_target_scorer.py`, `HONEST_AUDIT.md`, `universe.py`, `README.md` |
| Supplementary (later added) | **GSE129646** | Chat context summary, `gse129646_real_atac.json` (created 2026-03-30) |

### Why the conflict appeared

On 2026-03-30, an external agent delivered `GSE129646_MDAMB231_Par_BrM_LM_atac.counts.txt.gz` as part of a raw data delivery package. This file was parsed and described as "the ATAC data" in session context. However the **production pipeline** (`brm_data_loader.py`) already used `GSE205033` — which was downloaded 2026-03-28, is a 57 MB file with 261,222 peaks, and is the source for `atac_gene_scores.json` (261k peaks) and `atac_gene_level_scores.json` (443 genes).

### What is canonical

**GSE205033** (Biermann et al. *Cell* 2022, PMID 36113464) is the canonical ATAC source for Target-Lock v2.

Reasons:
1. It was used in the actual scored pipeline run that produced `BrM_Top_50_Targets.csv` and the v1 Assassin Scores.
2. It is melanoma BrM patient samples (in vivo, multi-patient) — correct for a BrM target universe.
3. `HONEST_AUDIT.md` (the internal truth document) explicitly lists it as REAL: `"35,241 MBM ATAC peaks from GSE205033 ✅ REAL | Biermann/Izar Cell 2022"`.

**GSE129646** (Cai & Greer 2020) is supplementary. It provides BrM2-vs-Parental differential accessibility for 6 gene-locus windows (MDA-MB-231 in vitro). It was processed into `gse129646_real_atac.json` on 2026-03-30 and is useful for per-gene locus validation but does NOT replace GSE205033 as the scoring signal.

---

## Rank Change — v1 vs v2 (imputed vs real CRISPRa)

The original `BrM_Top_50_Targets.csv` (run 2026-03-28) used an imputed CRISPRa signal derived from a partial MAGeCK run on the raw sgRNA counts. On 2026-03-30, the raw GSE237446 count matrices were parsed directly (brain vs lung guide-level counts).

| Gene | v1 Rank | v1 Assassin | v2 Rank | v2 Assassin | Delta |
|---|---|---|---|---|---|
| FAM72A | 1 | 0.7491 | 2 | 0.6549 | −0.0942 |
| SENP8 | 2 | 0.6857 | 3 | 0.6231 | −0.0626 |
| ATP10D | 3 | 0.6705 | 1 | 0.6743 | +0.0038 |
| BACE1 | 4 | 0.6611 | 8 | 0.5336 | −0.1275 |
| SLC25A32 | 5 | 0.6444 | 4 | 0.6168 | −0.0276 |
| SLC45A4 | 10 | 0.5784 | 5 | 0.5928 | +0.0144 |
| BIN1 | 15 | 0.5668 | 15 | 0.4634 | −0.1034 |

**Interpretation:** The v2 real-data run does not change the strategic picture. ATP10D, FAM72A, SENP8, SLC25A32 remain the top 4 with scores in the 0.62–0.67 range. BACE1 dropped because the real CRISPRa LFC (+7.07) is lower than its imputed value, and its 37.1 TPM GTEx MODERATE flag weights against it. The Evo2 FAM72A synthetic repressor generated on 2026-03-28 remains valid — FAM72A is still #2.

---

## Modal Evo2 Service — Version Lock

The live scoring service is `crispro-evo2-v9` in the `crispro-test` Modal workspace.

| Parameter | Value |
|---|---|
| App name | `crispro-evo2-v9` |
| Model | `evo2_1b_base` |
| GPU | A100 (SM8.0) |
| Scoring method | conditional log-likelihood (position-specific) |
| TE version | 1.13 (pinned — cuDNN 8 compat) |
| torch version | 2.3.0+cu121 |
| Patches | `patches.py` (torch.load + FP8 disable) |
| Live verification | 12/12 clinical variants scored 2026-03-28 |

**Do not bump TE version** without migrating to a cuDNN 9+ base image.  
**Do not pass `weights_only=False` natively** into vortex model_load classes — retain global patch.

---

## What changed, why, and what is now canonical

| Item | Was | Is now | Why |
|---|---|---|---|
| ATAC source | Ambiguous (GSE205033 vs GSE129646 in chat) | **GSE205033 (canonical)** | Used in production pipeline; in vivo patient data |
| ATAC role of GSE129646 | Described as "the real ATAC data" | Supplementary locus validation | Delivered 2026-03-30 after production run |
| #1 Assassin target | FAM72A (0.7491, v1) | ATP10D (0.6743, v2 real) | Real CRISPRa +12.16 confirms ATP10D screen hit |
| #1 safest target | FAM72A (0.28 TPM MINIMAL) | FAM72A still safest | GTEx v8 real data confirmed |
| CRISPRa signal | Imputed from partial MAGeCK | Real: 56,160 guides × 23 samples | GSE237446 raw counts parsed 2026-03-30 |
| data_loader_manifest.json | Missing | **Created: evo2_brm_manifest.json** | This document |
| MSK-MET brain filter | Described as applied | NOT applied (clinical file absent) | data_clinical_sample.txt not available |

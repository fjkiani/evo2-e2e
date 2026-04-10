# Evo2 Scoring — Use-Case Manifest

**Date:** 2026-04-10  
**Status:** CANONICAL — supersedes all prior informal notes on Evo2 signal interpretation  
**Manifest version:** aligns with `evo2_brm_manifest.json v1.1.0`

---

## Regime B Contract (validated)

**Service:** Modal `crispro-evo2-v9` → `Evo2Service.score_variant_conditional(wt_seq, mut_seq, var_idx)`  
**Returns:** `method: conditional_ll`, `delta_ll`, `ll_wt_conditional`, `ll_mut_conditional`

**Sequence contract (validated against TP53 R175H anchor):**
```
Assembly:  GRCh38
Source:    Ensembl REST API, forward genomic strand (no RC)
Window:    8192 bp
  start  = pos - 4095
  end    = pos + 4096
  length = 8192
var_idx  = 4095   ← variant at index 4095 (0-based)
```

**DO NOT use** `pos ± 4096` (produces 8193 bp, var_idx=4096, wrong conditional position — TP53 collapses from −0.4177 to −0.0086).  
**DO NOT** reverse-complement for minus-strand genes — breaks March 28 TP53 parity.

---

## Anchor Receipt (TP53 R175H)

| Field | Value |
|-------|-------|
| Gene | TP53 |
| Variant | p.R175H |
| GRCh38 | chr17:7674220 C>T |
| delta_ll | **−0.4177124500274658** |
| ll_wt_conditional | −1.2569301128387451 |
| ll_mut_conditional | −1.6746425628662110 |
| method | conditional_ll |
| model | evo2_1b_base |
| service | crispro-evo2-v9 (A100, Modal) |
| scored | 2026-03-28 |

Validate any new run against this anchor before trusting signs elsewhere.

---

## Regime A Contract (marginal LL — app-path)

**Service:** `score_variant_multi` → `alt_ll − ref_ll` across windows [1024, 2048, 4096, 8192]  
**Returns:** `min_delta` (most negative across windows), `window_used`  
**Scale:** ~25–50× smaller than Regime B for same variant  
**Symmetry fix:** deployed 2026-04-10 — `disable_symmetry=True` now uses `forward_min` directly  

**DO NOT mix scales or thresholds between Regime A and Regime B.**

---

## Use Case Classification

### ✅ VALIDATED — Regime B
**Rare pathogenic hotspots on conserved positions**

- delta_ll negative = evolutionarily surprising at scored locus = sequence-level disruption
- Signal is real and scaled: TP53 R175H −0.4177, comparable to FOXO3 −0.4673
- Correct for: TP53, BRCA1, PIK3CA active sites, known ClinVar pathogenic coding variants
- Weight in Target-Lock v2: 5% (Pearson r=0.016 vs CRISPRa brain LFC — orthogonal signal)

---

### ⚠️ NOT VALID WITHOUT CALIBRATION — Regime B

**Ancestral-allele variants and common coding changes**

Example: APOE alleles under identical Regime B contract (forward, var_idx=4095):

| Variant | rsID | GRCh38 | Allele | delta_ll | Clinical meaning |
|---------|------|--------|--------|----------|-----------------|
| APOE ε4 | rs429358 | 19:44908684 | T→C | **+0.7635** | AD risk ↑, longevity ↓ |
| APOE ε2 | rs7412 | 19:44908822 | C→T | **−0.8826** | AD risk ↓, longevity ↑ |

**ε2 − ε4 separation: −1.646 delta_ll** (magnitude ~1.65) — largest separation in this panel.

**Interpretation:** ε4 is the *ancestral* human allele. Across millions of species, the Cys at APOE position 112 is the evolutionary consensus — Evo2 finds it *expected* (positive delta_ll). ε3 (current reference) and ε2 are derived alleles that arose later in human evolution — Evo2 finds them *surprising* (negative delta_ll). The clinical AD risk of ε4 is a derived human population effect invisible to cross-species sequence constraint.

**Rule:** delta_ll sign does not map to clinical harm for common variants or ancestral alleles without a gene-specific calibration study.

---

### ❌ NOT VALID — Regime B (without separate evidence)

**Common longevity-associated SNPs**

| Variant | rsID | GRCh38 | delta_ll | Longevity association |
|---------|------|--------|----------|----------------------|
| FOXO3 | rs2802292 | 6:108587315 | **−0.4673** | Protective (epidemiological) |

FOXO3 scores more negative than TP53 R175H under Regime B. The longevity T allele is less probable than ancestral G given evolutionary context. The longevity association is a human cohort epidemiological finding, not a sequence-surprise effect. Common SNPs (MAF >5%) are by definition not evolutionarily surprising — moderate delta_ll values reflect position-level conservation, not variant-level pathogenicity.

**Rule:** Longevity associations require clinical cohort calibration. Evo2 conditional LL is not a longevity signal.

---

## What Evo2 Conditional LL Actually Measures

```
log P(mut_i | context_{1:i-1}) - log P(wt_i | context_{1:i-1})

Negative delta_ll = mut base is LESS probable than wt given left context
                  = evolutionarily SURPRISING at this position
                  = correlates with pathogenicity ONLY for rare variants
                    on conserved positions (TP53, BRCA1 etc.)

Positive delta_ll = mut base is MORE probable than wt
                  = may indicate ANCESTRAL allele (APOE ε4)
                  = does NOT mean clinically protective
```

The model has no concept of: human lifespan, disease phenotype, population founder effects, or which allele is "modern." It has deep knowledge of cross-species sequence conservation and nucleotide context probability.

---

## Operational Rules

1. **Anchor every new gene** against TP53 R175H before interpreting signs.
2. **Do not use delta_ll sign to infer clinical harm** for common variants (MAF >1%).
3. **Ancestral alleles** produce positive delta_ll — check allele frequency and evolutionary history before interpreting.
4. **Regime A and Regime B scales are incompatible.** Do not threshold one with the other's cutoffs.
5. **Longevity SNPs:** Evo2 is not a longevity model. Use epidemiological cohort data.
6. **RUO:** All scores are Research Use Only. No clinical interpretation without calibration study.

---

## Scoring Regime Quick Reference

| Question | Regime | Method | Scale |
|----------|--------|--------|-------|
| Is this cancer hotspot likely disruptive? | B | conditional_ll | −0.2 to −0.5 expected |
| Does this rare variant disrupt a conserved domain? | B | conditional_ll | Negative = yes |
| Is APOE ε4 more pathogenic than ε2? | Neither alone | requires calibration | Signs invert |
| Does this SNP affect longevity? | Neither | cohort epidemiology | Evo2 uninformative |
| Fast triage across a variant panel | A | marginal_ll_multi | −0.01 to −0.02 scale |

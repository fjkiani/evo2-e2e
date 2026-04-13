# PhD Handoff — CrisPRO BrM Pipeline
## Full Build-On Instructions — April 13, 2026

**Current state:** Class-aware composite AUROC = 0.7278. Target reached. All receipts on disk.  
**Repo:** `fjkiani/evo2-e2e` + `fjkiani/crispro-backend-v3`  
**Modal workspace:** `crispro-test` (credentials in `~/.modal.toml`)

---

## Active Endpoints (all live, no action needed)

| Service | URL | GPU | What it does |
|---|---|---|---|
| Evo2 conditional LL | `https://crispro-test--score-variant-conditional.modal.run` | A100 | POST `{wt_seq, mut_seq, var_idx}` → `delta_ll` |
| Evo2 mean LL | `https://crispro-test--crispro-evo2-v9-evo2service-score-variant.modal.run` | A100 | POST same → mean LL (Regime A) |
| Enformer brain DNase | `https://crispro-test--enformer-accessibility.modal.run` | T4 | POST `{sequence, track_type}` → brain accessibility |
| v3 backend | Render (crispro-backend-v3) | CPU | Full oncology pipeline + `/brm/*` endpoints |

**Scoring contract (never deviate):**
```
window   = 8192 bp
start    = pos - 4095
end      = pos + 4096
var_idx  = 4095
strand   = forward GRCh38 Ensembl
anchor   = TP53 R175H chr17:7674220 C>T → -0.4177124500274658
```
Validate the anchor before trusting any new variant run.

---

## Where Everything Lives

```
evo2-e2e/
  evo2_brm_manifest.json          ← v1.4.0 — canonical data provenance
  MANIFEST_NOTES.md               ← ATAC discrepancy log + key decisions
  HONEST_AUDIT.md                 ← tautology confession, AUROC history
  docs/evo2_scoring_use_case_manifest.md  ← what Evo2 can/cannot do
  data/brain_met/
    BrM_Assassin_Scored_v3_realdata.csv   ← 50 targets, all signals real
    gse237446_real_lfc.json               ← CRISPRa 56,083 genes
    gse205033_real_atac_mbm_ecm.json      ← ATAC 14,269 genes (50kb window)
    atac_fill_missing_onc.json            ← 5 missing ONC genes (500kb window)
    mskmet_brain_burden_real.json         ← MSK-MET 2,921 brain samples
    gtex_brain_safety.json                ← GTEx v8 50 targets
    fam72a_synthetic_repressor.fasta      ← Evo2-generated repressor (4 candidates)
    af3_cas9_rnp_*.json                   ← AF3 Cas9 RNP specs Top 5
  data/validation/
    live_variant_scores.json              ← 12 BrM + 3 longevity variants scored
    honest_auroc_panel.json               ← 29-gene panel (9 scored positives)
    class_aware_composite_final.json      ← AUROC=0.7278, all signals
    adversarial_negatives_scored.json     ← 4 adversarial negatives + Evo2 scores
    regime_b_longevity_run.json           ← APOE e2/e4, FOXO3 Regime B scores

crispro-backend-v3/
  layer0_contracts/evo2_scoring.py        ← Layer 0 contract (RUO, typed)
  layer0_contracts/patient_state.py       ← SomaticMutation.sequence_window field
  layer1_engines/evo2_scorer/__init__.py  ← thin httpx client to Modal
  layer1_engines/brm_engine/__init__.py   ← Assassin Score engine
  layer3_orchestrator/pipeline.py         ← Phase 6.5 EVO2_SCORING_ENABLED
  layer4_output/tumor_board_formatter.py  ← evo2_evidence in output
  tests/test_ayesha_e2e.py               ← MERGE GATE (must pass before push)
  tests/test_layer0_evo2_contracts.py    ← 30 tests, RUO invariants
  tests/test_evo2_scorer_engine.py       ← 22 tests, sequence_window
  tests/test_pipeline_evo2_flagged.py    ← 15 tests, flag behavior
```

---

## Open Items — Priority Order

### P0 — These block paper submission

**[UNRESOLVED-002] T0 pellet normalization**
- What: GSE237446 T0 pellet (GSM7616089) has no supplementary count file
- Why it matters: current CRISPRa LFC = brain vs lung with no T0 baseline
- Fix: Download SRA, align, count, run MAGeCK 3-condition
```bash
# SRA accession: SRX21037341
fastq-dump --split-files SRX21037341
# Calabrese A library: download from Addgene #92379
bowtie2 -x calabrese_a -U SRR*.fastq | samtools view -bS > t0.bam
mageck count -l calabrese_a.csv -n brm_3cond \
  --sample-label "T0,Brain,Lung" \
  --fastq t0.bam brain_*.bam lung_*.bam
mageck rra -k brm_3cond.count.txt -n brm_rra
```
- Output needed: replace `gse237446_real_lfc.json` with MAGeCK RRA scores

**[UNRESOLVED-003] AlphaFold3 4-chain submissions**
- Specs ready in `data/brain_met/af3_cas9_rnp_*.json`
- FAM72A, ATP10D, BACE1, BIN1 not yet submitted (SLC25A32 partial iptm=0.42)
- Go to: alphafoldserver.com → New Prediction → Upload JSON
- 10 jobs/day limit. Submit one per day for 4 days.
- When iptm returns: update `structure` component in Assassin Score formula
- BACE1 and FAM72A results will change Assassin Score ranking

**Prospective AUROC validation**
- Current AUROC (0.7278) is on the training panel. Not prospective.
- Fix: Pull TCGA-BRCA brain metastasis patients (n≈40) from cBioPortal
- Score them blind through Target-Lock v2
- Compare enrichment of top-ranked genes in BrM vs primary cohort
- This converts "useful signal" into a citable validation number

---

### P1 — These improve the paper materially

**Fix CRISPRa signal in scorer.py**
- Current: `brm_target_scorer.py` uses raw LFC
- Fix: use `abs(lfc)` for the composite, routing by gene class
- One-line change, but requires re-running the full 56,083-gene scoring pass
- Re-run `scripts/run_brm_pipeline.py` and save new `BrM_Top_50_Targets.csv`

**sequence_window population for v3**
- `SomaticMutation.sequence_window` is defined but callers must pre-populate
- Fix: build an Ensembl REST client that fetches sequences from HGVSp/HGVSc
```python
# Template for the fetcher:
def populate_sequence_window(mutation, window=8192):
    """Fetch hg38 sequence for a somatic mutation via Ensembl REST."""
    gene_id = get_ensembl_gene_id(mutation.gene)
    variant_pos = get_variant_position(gene_id, mutation.hgvs_c)
    chrom, pos = variant_pos
    wt_seq = fetch_ncbi_window(chrom, pos, window)
    alt_base = mutation.hgvs_p.split(">")[1] if ">" in mutation.hgvs_p else None
    if alt_base:
        var_idx = window // 2
        mut_seq = wt_seq[:var_idx] + alt_base + wt_seq[var_idx+1:]
        mutation.sequence_window = {"wt_seq": wt_seq, "mut_seq": mut_seq, "var_idx": var_idx}
    return mutation
```
- Wire this upstream of the v3 PatientState submission
- Then set `EVO2_SCORING_ENABLED=true` in Render env and re-test Ayesha

**[UNRESOLVED-004] SENP8 synthetic repressor**
- SENP8 is #3 in final Assassin Score, 1.4 TPM LOW safety
- Run the same Evo2 steered generation as FAM72A:
```bash
MODAL_PROFILE=crispro-test python3 << 'EOF'
# Same as fam72a_evo2_steered.py but:
# 1. Fetch SENP8 TSS: NC_000003.12, chr3q22.3
# 2. Scan promoter for SP1/E2F1/CREB1 motifs (JASPAR scan)
# 3. Run svc.generate.remote() with SENP8 exon1 context
# 4. Save to data/brain_met/senp8_synthetic_repressor.fasta
EOF
```

**Harder adversarial negatives with known brain expression**
- Current negatives (MYC, KRAS etc.) are too easy — AUROC improved with adversarials
- The true hard negatives are brain-expressed non-BrM genes: GRIN2A, PTPRT, NF1
- Score all three through conditional LL (endpoint is live)
- Recompute AUROC on 32-gene panel (29 + 3 hard adversarials)
- If AUROC holds at >0.70: the architecture is robust
- If AUROC drops below 0.65: the ATAC signal is doing the heavy lifting

---

### P2 — These are the next capability expansions

**Breast met pipeline (Track B — started, not complete)**
- `breast_met/` skeleton exists. GSE184869 is already downloaded.
- PhD found: no validated breast-to-brain CRISPRa screen in GEO
- Substitute signal: GSE184869 upregulation as primary (BrM vs primary breast matched RNA-seq)
- Target universe: ESR1, PIK3CA, ERBB2, PTEN, CDH1, MYC, CCND1 as positives
- Run Target-Lock v2 adapted for breast: swap CRISPRa→RNA as primary signal
- Report: breast BrM AUROC vs brain BrM AUROC (cross-cancer comparison)

**Enformer as systematic 6th signal**
- Currently used for 29-gene panel only
- Extend to all 50 top Assassin Score targets
- Enformer predicts 5000 tracks from 196,608bp windows
- Brain-specific tracks: H3K27ac (brain), DNase (brain), ATAC (brain)
- Track indices from Enformer paper Table S1:
  - DNase brain: tracks 121-244
  - H3K27ac brain: tracks 684-723
- Score all 50 targets' TSS ±100kb
- Report: Enformer vs GSE205033 ATAC correlation (validates both signals)

**Clinical AUROC with VUS panel**
- Current 12 scored clinical variants: TP53, PIK3CA, BRCA1, PTEN etc.
- Expand to 50 variants from MSK-MET brain_muts (use `HGVSp_Short` field)
- Populate sequence windows via Ensembl
- Score all 50 through conditional LL endpoint
- Compare to ClinVar pathogenicity labels → calibrated AUROC for clinical use

---

## Scoring Contract Reference

### Evo2 Conditional LL (Regime B)
```python
import urllib.request, json

def score_variant_regime_b(chrom, pos, ref_pos_base, alt_base, window=8192):
    """Score a variant through crispro-evo2-v9 conditional LL."""
    ACC = {"chr17":"NC_000017.11","chr1":"NC_000001.11",...}  # full map in fetch code
    start = pos - 4095; stop = pos + 4096
    url = (f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
           f"?db=nuccore&id={ACC[chrom]}&rettype=fasta&retmode=text"
           f"&seq_start={start}&seq_stop={stop}")
    wt_seq = ''.join(l for l in urllib.request.urlopen(url).read().decode().split('\n')
                     if not l.startswith('>')).upper()
    var_idx = 4095
    mut_seq = wt_seq[:var_idx] + alt_base + wt_seq[var_idx+1:]
    
    payload = json.dumps({"wt_seq": wt_seq, "mut_seq": mut_seq, "var_idx": var_idx}).encode()
    req = urllib.request.Request(
        "https://crispro-test--score-variant-conditional.modal.run",
        data=payload, headers={"Content-Type": "application/json"}
    )
    return json.loads(urllib.request.urlopen(req, timeout=120).read())
    # Returns: {"delta_ll": float, "ll_wt_conditional": float, "ll_mut_conditional": float,
    #           "method": "conditional_ll", "model": "evo2_1b_base"}

# Validate anchor before any run:
result = score_variant_regime_b("chr17", 7674220, "C", "T")
assert abs(result["delta_ll"] - (-0.4177124500274658)) < 0.01, "ANCHOR MISMATCH — stop"
```

### Enformer Brain Accessibility
```python
def score_enformer_dnase(chrom, tss_pos, window=5000):
    """Score TSS brain DNase accessibility via Enformer."""
    ACC = {"chr17":"NC_000017.11",...}
    start=max(1,tss_pos-window//2); stop=tss_pos+window//2
    url = (f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
           f"?db=nuccore&id={ACC[chrom]}&rettype=fasta&retmode=text"
           f"&seq_start={start}&seq_stop={stop}")
    seq = ''.join(l for l in urllib.request.urlopen(url).read().decode().split('\n')
                  if not l.startswith('>')).upper()
    payload = json.dumps({"sequence": seq, "track_type": "dnase"}).encode()
    req = urllib.request.Request(
        "https://crispro-test--enformer-accessibility.modal.run",
        data=payload, headers={"Content-Type": "application/json"}
    )
    return json.loads(urllib.request.urlopen(req, timeout=120).read())
    # Returns: {"accessibility_score": float, "track_type": "dnase", "n_tracks": 124}
```

### Class-Aware Composite Formula (AUROC=0.7278)
```python
def class_aware_score(gene, abs_crispr_norm, atac_norm, evo2_norm, enformer_norm):
    """Route gene through mechanistically appropriate signal weights."""
    TSG_GENES = {"TP53","BRCA1","PTEN","SMARCA4","KMT2C","CDKN2A","SMAD4","APC","CDH1"}
    if gene in TSG_GENES:
        # CRISPRa activation magnitude is the correct TSG signal
        return 0.60*abs_crispr_norm + 0.20*atac_norm + 0.15*evo2_norm + 0.05*enformer_norm
    else:
        # ATAC chromatin accessibility is the correct ONC/functional signal
        return 0.25*abs_crispr_norm + 0.35*atac_norm + 0.20*evo2_norm + 0.20*enformer_norm

# NOTE: all signals must be min-max normalized TOGETHER before applying weights
# NOTE: evo2 direction: MORE NEGATIVE delta_ll → HIGHER score (negate before norming)
# NOTE: CRISPRa uses ABSOLUTE VALUE (|lfc|), not signed LFC
```

---

## Ayesha Gate — Non-Negotiable

Every push to `crispro-backend-v3` must pass:
```bash
cd crispro-backend-v3
pip install -r requirements.txt  # structlog must be installed
python -m pytest tests/test_ayesha_e2e.py tests/test_layer0_evo2_contracts.py \
  tests/test_evo2_scorer_engine.py tests/test_pipeline_evo2_flagged.py \
  tests/test_layer1/test_synthetic_lethality/ -q
# Must show: 310 passed, 0 failed
```

If Ayesha fails → do not push. Fix first.

---

## What Evo2 Can and Cannot Do (Locked)

### CAN (validated):
- Rare pathogenic hotspots at conserved positions (TP53, BRCA1, PIK3CA active sites)
- Discriminate BrM oncogenes from general oncogenes (AUROC=0.70 in oncogene class)
- Generate synthetic repressor sequences via steered generation (FAM72A done, SENP8 pending)

### CANNOT (proven):
- Discriminate BrM-relevant TSG loss from non-BrM TSG loss (AUROC=0.56)
- Detect functional cascade mechanisms (BACE1→EGFR→MEK is invisible to sequence constraint)
- Score longevity/common population SNPs meaningfully (FOXO3 at -0.467 ≠ longevity signal)
- Infer clinical direction from delta_ll sign for common/ancestral alleles (APOE ε4 is positive because it is ancestral — not protective)

### Weight in Target-Lock v2: **5% — LOCKED**
Do not increase without AUROC > 0.70 on a prospective holdout. Rationale locked in manifest v1.4.0.

---

## AUROC History — Never Lose This

| AUROC | Date | What it was | Why it changed |
|---|---|---|---|
| 0.976 | March 2026 | Tautology | mission_fit_discount encoded labels directly |
| 0.98 | March 2026 | Same tautology | discount-alone AUROC = 1.0 |
| 0.58 | March 2026 | Honest attempt | Hard negatives were hematologic (too easy) |
| **0.5778** | April 13 | Honest — Track A | Solid-tumor negatives, consistent conditional LL |
| **0.7278** | April 13 | Class-aware composite | ATAC fill + Enformer + gene class routing |

The 0.7278 is not comparable to the 0.5778 — it uses a different formula and different signal composition. Both are honest. The 0.5778 is the single-signal Evo2 baseline. The 0.7278 is the multi-modal class-aware composite. Both should appear in the paper.

---

## Manifest Versioning

`evo2_brm_manifest.json` is currently v1.4.0.

When to bump:
- **Minor (1.x.0):** new data signal, closed unresolved item, new scored panel
- **Major (x.0.0):** formula weight change, threshold change, negative set replacement

After any AUROC run: update `honest_auroc_v2` block in manifest.  
After paper submission: freeze the manifest. No edits after submission.

---

## Remaining Open Items (from manifest v1.4.0)

| ID | Item | Status | Effort |
|---|---|---|---|
| UNRESOLVED-002 | T0 pellet normalization | Open | 1-2 days (needs SRA download) |
| UNRESOLVED-003 | AF3 4-chain submissions (FAM72A, ATP10D, BACE1, BIN1) | Open | 4 days (10/day limit) |
| UNRESOLVED-004 | SENP8 synthetic repressor | Open | 2 hours |
| NEW | Prospective AUROC validation (TCGA-BRCA holdout) | Open | 3-5 days |
| NEW | sequence_window Ensembl fetcher for v3 | Open | 1 day |
| NEW | Breast met pipeline (Track B) | Open | 3-5 days |

---

## What NOT To Do

1. **Do not use raw CRISPRa LFC in the composite.** TSGs have negative LFC (activation suppresses BrM). Use `abs(lfc)` for the composite signal.

2. **Do not mix Regime A and Regime B scores.** Mean LL (Regime A, `score_variant_multi`) and conditional LL (Regime B, `score_variant_conditional`) are different computations. Spearman ρ=0.35 between them. Never mix in the same AUROC.

3. **Do not use GSE129646 for ATAC.** That is MDA-MB-231 in vitro (Cai & Greer 2020). The canonical ATAC source is GSE205033 (Biermann/Izar Cell 2022, melanoma patient-derived BrM cell lines).

4. **Do not score longevity variants through delta_ll and interpret the sign as clinical direction.** APOE ε4 is positive because it is the ancestral allele. The model is not wrong. The question is wrong.

5. **Do not push to crispro-backend-v3 without running Ayesha first.** 310 tests, 0 failures. Non-negotiable.

6. **Do not inflate maturity labels.** VALIDATED means external prospective validation. The BrM pipeline is currently IMPLEMENTED with promising internal validation. It becomes VALIDATED after the prospective TCGA-BRCA holdout.

7. **Do not change the Evo2 5% weight** without a new calibration study showing AUROC > 0.70 on a prospective holdout. The current 5% is locked in the manifest with empirical justification (r=0.016 Pearson, AUROC=0.5778 alone).

# Agent Audit — April 2, 2026

Files received: GSE237446 (×2 new), H3K27ac bedgraphs (×2), data_mutations (×3 new copies)
Audit method: byte-level hash comparison, direct column inspection, line counts, format verification

---

## PROMPT A — GSE237446 CRISPRa Screen

### What the agent said
- Brain file: 576 KB, 12 brain samples
- Lung file: 1.5 MB, 11 lung samples
- T0 pellet: no standalone file, T0 column may be embedded
- GeneID format: gene symbol embedded before `_CalA_`

### What audit found

**Files are IDENTICAL to the March 30 copies.** Both the `-1-` renamed files match byte-for-byte:
```
BRAINS old vs new: IDENTICAL (md5: d1ca7e5f...)
LUNGS  old vs new: IDENTICAL (md5: 804290a6...)
```
No new data was delivered. The agent re-downloaded the same files.

**Structure confirmed correct:**
- Brain: 56,160 guides × 12 samples — all columns named `Brain 1` through `Brain 12`
- Lung:  56,160 guides × 11 samples — all columns named `Lung 1` through `Lung 11`
- 85.7% of brain rows are all-NA (sparse, expected for CRISPRa screen)
- 24.9% of lung rows are all-NA

**T0 pellet:** No T0 column in either matrix. Agent's assessment is correct. GSM7616089 has no supplementary file.

**GeneID format — CORRECTION to previous audit:**

Prior report (March 30) stated "GeneID is GENESYMBOL_CalA_GUIDEINDEX — confirmed." This was *directionally* correct but needed refinement. The format is:

- **99.9% of guides** (56,083/56,160): GENE_SYMBOL_CalA_GUIDE_INDEX (e.g. `A1BG_CalA_1234`)
- **0.1% of guides** (77/56,160): NUMERIC_CalA_GUIDE_INDEX (e.g. `44621_CalA_12123`)

The numeric IDs (range 44621–44896) do not match standard human Entrez gene IDs (which are scattered 1–100k). These 77 guides are likely internal library controls or lncRNA/non-coding entries. They affect <0.1% of genes.

**Previous pipeline impact:** Negligible. 77 guides out of 56,160 had unparseable gene symbols. All went into a `44621` etc. bin which would not appear in any top-N target list. The ranked gene results are unaffected.

**Status: NO DATA ISSUE. Previous pipeline results valid.**

---

## PROMPT B — H3K27ac bedgraphs (GSM3530739, GSM3530746)

### What the agent said
- GSE124379 H3K27ac ChIP-seq
- MDA231 (parental) + LM2 (lung met) cell lines
- No BrM2 cell line present
- Not peak calls — need MACS2 to threshold
- No ATAC-seq in this dataset

### What audit found

**Agent's description is factually correct.** The bedgraph format and cell lines are accurately described.

**The files themselves:**
- GSM3530739 (MDA231): 33,961,279 intervals, signal = integer read pileup (1 per read position)
- GSM3530746 (LM2): 32,777,488 intervals
- Format: `chr start end signal` — raw coverage, NOT fold-enrichment, NOT RPM-normalized
- Chromosomes are hg19 (chr1, chr2, etc. — no 'chr' prefix would indicate hg38)

**CRITICAL FINDINGS — THREE PROBLEMS:**

**Problem 1: Wrong GEO accession.** The agent cites GSE124379 (Li et al. 2019, H3K27ac ChIP-seq). The canonical ATAC source locked in `evo2_brm_manifest.json` is **GSE205033** (Biermann/Izar Cell 2022, melanoma BrM ATAC-seq). These are entirely different datasets from different papers.

**Problem 2: Wrong data type.** These are H3K27ac ChIP-seq pileup tracks. They cannot substitute for ATAC-seq open chromatin data. H3K27ac marks active enhancers; ATAC marks accessible chromatin. They are correlated but not equivalent. Using H3K27ac in place of ATAC for Target-Lock v2 would require re-validating the ATAC component weight (0.20) against a different signal type.

**Problem 3: Wrong cell lines.** MDA231 (parental) and LM2 (lung met) are not brain metastasis models. The canonical ATAC source (GSE205033) uses melanoma BrM *patient* samples. LM2 is a lung metastasis model. No BrM2 column exists in this dataset — the agent confirmed this but still delivered the files as the ATAC replacement.

**These files serve no purpose in the current pipeline.** They should not be integrated.

**The canonical ATAC source remains GSE205033** (already processed in `data/brain_met/atac_gene_scores.json`).

**Status: FILES REJECTED. Do not integrate. GSE205033 remains canonical.**

---

## PROMPT C — MSK-MET data_mutations.txt

### What the agent said
- 77.8 MB file from Zenodo 5801902
- Columns confirmed: Hugo_Symbol, Tumor_Sample_Barcode, HGVSp_Short, etc.
- Brain filter: use `METASTATIC_SITE` column
- Python filter: `clinical['METASTATIC_SITE'].str.contains('Brain', na=False)`

### What audit found

**THREE NEW FILES ARE IDENTICAL TO THE MARCH 30 COPY:**
```
OLD (Mar 30):   3d13d561785b8cd6...
NEW -1 (Apr 2): 3d13d561785b8cd6...
NEW -2 (Apr 2): 3d13d561785b8cd6...
NEW -3 (Apr 2): 3d13d561785b8cd6...
```
The agent downloaded and delivered 3 identical copies of the same file that already existed. No new data.

**Column structure confirmed correct:**
- 122 columns, all tab-separated
- HGVSp_Short: ✅ present (needed for sequence_window in v3 Evo2 pipeline)
- First mutation: TP53 p.R273C sample=P-0028912-T01-IM6

**CRITICAL ERROR in agent's filter instruction:**

The agent says to filter with:
```python
brain_samples = clinical[
    clinical['METASTATIC_SITE'].str.contains('Brain', na=False)
]['SAMPLE_ID'].tolist()
```

This is **wrong table, wrong column.** `METASTATIC_SITE` is not in `data_mutations.txt`. It is in `data_clinical_sample.txt`. `data_mutations.txt` has 122 columns and `METASTATIC_SITE` is not one of them.

**Correct approach (two-step join):**
```python
# Step 1: Get brain-met sample IDs from the clinical file
import pandas as pd
clinical = pd.read_csv('data_clinical_sample.txt', sep='\t', comment='#')
brain_samples = clinical[
    clinical['METASTATIC_SITE'].str.contains('Brain', na=False)
]['SAMPLE_ID'].tolist()

# Step 2: Filter mutations by those sample IDs
mutations = pd.read_csv('data_mutations.txt', sep='\t')
brain_muts = mutations[mutations['Tumor_Sample_Barcode'].isin(brain_samples)]
```

`data_clinical_sample.txt` is the missing file. Without it, the brain-specific filter cannot be applied. This was `UNRESOLVED-001` in `evo2_brm_manifest.json` and remains unresolved.

**Status: NO NEW DATA. Filter instruction has wrong-table error. Unresolved-001 still open.**

---

## Summary Table

| Prompt | Files Received | New Data? | Agent Accurate? | Action Required |
|--------|---------------|-----------|-----------------|-----------------|
| A: GSE237446 | 2 files (renamed copies) | ❌ Identical to Mar 30 | ✅ Mostly correct | Fix numeric GeneID note (0.1%) |
| B: H3K27ac bedgraphs | 2 files (192MB, 185MB) | ✅ New files | ⚠️ Factually correct but wrong dataset | Reject files. GSE205033 stays canonical. |
| C: MSK-MET mutations | 3 files | ❌ All identical to Mar 30 | ❌ Wrong filter table | Get data_clinical_sample.txt |

## Open Items (unchanged from evo2_brm_manifest.json)

- **UNRESOLVED-001**: `data_clinical_sample.txt` still not available — MSK-MET brain filter still cannot be applied
- **UNRESOLVED-002**: GSE237446 T0 pellet still has no count data (SRA reads only)
- **UNRESOLVED-003**: AF3 jobs for FAM72A, ATP10D, BACE1, BIN1 pending submission
- **UNRESOLVED-004**: SENP8 repressor generation pending
- **NEW UNRESOLVED-005**: H3K27ac bedgraphs (GSE124379) delivered in place of ATAC — do not integrate

## What to tell the agent

1. Re-download is not audit. Three identical copies of `data_mutations.txt` is not progress.
2. GSE124379 H3K27ac ChIP-seq tracks are not a substitute for ATAC-seq. Wrong paper, wrong data type, wrong cell line. Do not integrate.
3. The `METASTATIC_SITE` brain filter belongs in `data_clinical_sample.txt`, not `data_mutations.txt`. The filter code you provided will crash or silently return nothing on the mutations file.
4. Obtain `data_clinical_sample.txt` from Zenodo 5801902 to resolve UNRESOLVED-001.
5. GSE205033 remains the canonical ATAC source. No action needed on ATAC unless `data_clinical_sample.txt` is delivered and brain-specific ATAC is being sought from a different dataset.

"""
breast_met/targets/universe.py
================================
Breast Metastasis Target Universe.

Extends the same BrMGene pattern from brain_met/targets/universe.py.
The architecture is identical — only the biology differs.

Key breast met biology:
  - ER+ BrM: ESR1 acquired mutations (endocrine resistance)
  - HER2+: ERBB2 amplification drives aggressive brain colonization
  - TNBC: BRCA1/2, immune-cold, high MMP expression
  - Bone met: different from brain — RANkL/OPG axis, osteolytic
  - Liver met: different from brain — portal circulation, Kupffer cells
  - Brain met from breast: AURORA cohort data (9 clinical variants)

Datasets:
  - AURORA breast BrM study
  - TCGA-BRCA (breast primary mutations)
  - Yale T-DXd ADC resistance cohort (scripts/yale_tdzd/)
  - Our scored variants: ESR1 Y537S/D538G, BRCA1 T1685A, PIK3CA H1047R/E545K

Steps (breast-met specific, 6 steps):
  1. primary_escape_breast   — escape from breast primary
  2. lymphatic_dissemination — lymph node involvement (different from general)
  3. distant_hematogenous    — blood-borne spread to distant sites
  4. site_selection          — organ tropism selection (brain vs bone vs liver)
  5. organ_colonization      — establish in target organ
  6. therapy_resistance      — acquired resistance to endocrine/HER2/chemo therapy
"""
from dataclasses import dataclass, field
from typing import List, Optional, Dict

@dataclass
class BrcaMGene:
    """Breast metastasis gene — same schema as BrMGene for pipeline compatibility."""
    symbol: str
    chrom: str
    pos: int
    ref: str
    alt: str
    primary_steps: List[str]
    secondary_steps: List[str]
    is_positive: bool
    label_source: str
    notes: str
    aurora_variant: Optional[str] = None
    her2_relevant: bool = False
    er_relevant: bool = False
    tnbc_relevant: bool = False
    brca_relevant: bool = False


BREAST_MET_UNIVERSE: List[BrcaMGene] = [

    # ─── ER+ BREAST MET DRIVERS ─────────────────────────────────────────────
    BrcaMGene(
        symbol="ESR1", chrom="chr6", pos=152098791, ref="C", alt="G",
        primary_steps=["therapy_resistance", "organ_colonization"],
        secondary_steps=["site_selection"],
        is_positive=True, er_relevant=True,
        label_source="AURORA",
        notes="ESR1 Y537S: endocrine resistance acquired in ER+ BrM; ligand-independent ER",
        aurora_variant="p.Y537S",
    ),
    BrcaMGene(
        symbol="ESR1", chrom="chr6", pos=152098791, ref="T", alt="A",
        primary_steps=["therapy_resistance"],
        secondary_steps=["organ_colonization"],
        is_positive=True, er_relevant=True,
        label_source="AURORA",
        notes="ESR1 D538G: ligand-independent ER activation in BrM; acquired variant",
        aurora_variant="p.D538G",
    ),
    BrcaMGene(
        symbol="PIK3CA", chrom="chr3", pos=179234297, ref="A", alt="G",
        primary_steps=["primary_escape_breast", "organ_colonization"],
        secondary_steps=["therapy_resistance"],
        is_positive=True, er_relevant=True,
        label_source="MSK-MET",
        notes="PI3K H1047R: ER+ BrM; top Evo2 score (delta_ll=-0.615); common in ER+ breast",
        aurora_variant="p.H1047R",
    ),
    BrcaMGene(
        symbol="PIK3CA", chrom="chr3", pos=179203765, ref="A", alt="G",
        primary_steps=["primary_escape_breast"],
        secondary_steps=["therapy_resistance"],
        is_positive=True, er_relevant=True,
        label_source="MSK-MET",
        notes="PI3K E545K: helical domain hotspot in ER+ breast",
        aurora_variant="p.E545K",
    ),

    # ─── HER2+ BREAST MET ───────────────────────────────────────────────────
    BrcaMGene(
        symbol="ERBB2", chrom="chr17", pos=39688484, ref="A", alt="G",
        primary_steps=["primary_escape_breast", "distant_hematogenous"],
        secondary_steps=["organ_colonization"],
        is_positive=True, her2_relevant=True,
        label_source="NCT04740918",  # HER2+ BrM trial
        notes="HER2 amplification/mutation driver in aggressive BrM breast cancer",
    ),

    # ─── BRCA1/2 (TNBC + Hereditary) ────────────────────────────────────────
    BrcaMGene(
        symbol="BRCA1", chrom="chr17", pos=43057051, ref="A", alt="G",
        primary_steps=["primary_escape_breast"],
        secondary_steps=["organ_colonization", "therapy_resistance"],
        is_positive=True, brca_relevant=True, tnbc_relevant=True,
        label_source="ClinVar/AURORA",
        notes="BRCA1 VUS T1685A; hereditary BrM; ACMG scoring target; Evo2 delta_ll=-0.095",
        aurora_variant="p.T1685A",
    ),
    BrcaMGene(
        symbol="PTEN", chrom="chr10", pos=89692905, ref="G", alt="A",
        primary_steps=["organ_colonization", "therapy_resistance"],
        secondary_steps=["primary_escape_breast"],
        is_positive=True,
        label_source="MSK-MET",
        notes="PTEN R130Q LOF: immune evasion + PI3K; MSK-MET enriched in BrM",
        aurora_variant="p.R130Q",
    ),
    BrcaMGene(
        symbol="TP53", chrom="chr17", pos=7674220, ref="C", alt="T",
        primary_steps=["primary_escape_breast", "organ_colonization"],
        secondary_steps=["therapy_resistance"],
        is_positive=True,
        label_source="AURORA",
        notes="TP53 R175H: 2× enriched in BrM; AURORA hotspot; Evo2 delta_ll=-0.418",
        aurora_variant="p.R175H",
    ),
    BrcaMGene(
        symbol="KMT2C", chrom="chr7", pos=151835898, ref="C", alt="T",
        primary_steps=["primary_escape_breast"],
        secondary_steps=["organ_colonization"],
        is_positive=True,
        label_source="AURORA",
        notes="KMT2C R4854* truncation → H3K4me1 loss → MMP3 → invasion; AURORA cohort",
        aurora_variant="p.R4854*",
    ),

    # ─── ADC RESISTANCE (Yale T-DXd cohort) ─────────────────────────────────
    BrcaMGene(
        symbol="ABCB1", chrom="chr7", pos=87531302, ref="G", alt="A",
        primary_steps=["therapy_resistance"],
        secondary_steps=[],
        is_positive=True, her2_relevant=True,
        label_source="PMID:32816860",  # ADC efflux resistance
        notes="MDR1/P-gp efflux pump; T-DXd resistance via payload pumping; Kill Chain DRUG_EFFLUX",
    ),
    BrcaMGene(
        symbol="SLFN11", chrom="chr17", pos=31962468, ref="C", alt="T",
        primary_steps=["therapy_resistance"],
        secondary_steps=[],
        is_positive=True, tnbc_relevant=True,
        label_source="PMID:32816860",
        notes="SLFN11 silencing → PARPi+platinum resistance; Yale T-DXd cohort marker",
    ),

    # ─── HARD NEGATIVES ────────────────────────────────────────────────────
    BrcaMGene(
        symbol="IDH1", chrom="chr2", pos=209113112, ref="C", alt="T",
        primary_steps=[], secondary_steps=[],
        is_positive=False,
        label_source="Hard negative — glioma, not breast met",
        notes="IDH1 R132H = glioma driver, not breast cancer biology",
    ),
    BrcaMGene(
        symbol="FLT3", chrom="chr13", pos=28033987, ref="A", alt="T",
        primary_steps=[], secondary_steps=[],
        is_positive=False,
        label_source="Hard negative — AML, not breast",
        notes="AML driver; no breast cancer biology",
    ),
    BrcaMGene(
        symbol="ABL1", chrom="chr9", pos=130713881, ref="G", alt="A",
        primary_steps=[], secondary_steps=[],
        is_positive=False,
        label_source="Hard negative — CML",
        notes="BCR-ABL CML driver; no breast metastasis role",
    ),
    BrcaMGene(
        symbol="JAK2", chrom="chr9", pos=5073770, ref="G", alt="T",
        primary_steps=[], secondary_steps=[],
        is_positive=False,
        label_source="Hard negative — MPN",
        notes="JAK2 V617F myeloproliferative; no breast met biology",
    ),
    BrcaMGene(
        symbol="RUNX1", chrom="chr21", pos=34792956, ref="C", alt="T",
        primary_steps=[], secondary_steps=[],
        is_positive=False,
        label_source="Hard negative — AML/MDS",
        notes="RUNX1 leukemia driver; not a solid tumor metastasis gene",
    ),
]

POSITIVE_GENES = [g for g in BREAST_MET_UNIVERSE if g.is_positive]
NEGATIVE_GENES = [g for g in BREAST_MET_UNIVERSE if not g.is_positive]

UNIVERSE_STATS = {
    "total": len(BREAST_MET_UNIVERSE),
    "positive": len(POSITIVE_GENES),
    "negative": len(NEGATIVE_GENES),
    "er_relevant": sum(1 for g in BREAST_MET_UNIVERSE if g.er_relevant),
    "her2_relevant": sum(1 for g in BREAST_MET_UNIVERSE if g.her2_relevant),
    "brca_relevant": sum(1 for g in BREAST_MET_UNIVERSE if g.brca_relevant),
    "with_aurora": sum(1 for g in BREAST_MET_UNIVERSE if g.aurora_variant),
    "steps": 6,
}


def as_scoring_input(gene: BrcaMGene, step: str) -> dict:
    """Same interface as brain_met for pipeline compatibility."""
    return {
        "gene": gene.symbol,
        "step": step,
        "chrom": gene.chrom,
        "pos": gene.pos,
        "ref": gene.ref,
        "alt": gene.alt,
        "hgvs_p": gene.aurora_variant,
        "consequence": None,
        "primary_steps": gene.primary_steps,
        "secondary_steps": gene.secondary_steps,
    }


if __name__ == "__main__":
    print(f"Breast Met Target Universe: {UNIVERSE_STATS}")
    print(f"\nER+ relevant: {[g.symbol for g in BREAST_MET_UNIVERSE if g.er_relevant and g.is_positive]}")
    print(f"HER2 relevant: {[g.symbol for g in BREAST_MET_UNIVERSE if g.her2_relevant and g.is_positive]}")
    print(f"Hard negatives: {[g.symbol for g in NEGATIVE_GENES]}")

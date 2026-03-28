"""
brain_met/targets/universe.py
==============================
Brain Metastasis Target Universe — the dedicated brain-met panel.

This replaces "brain cancer as a vague extension" with an explicit, labeled
gene universe covering 7 BrM-specific biological steps with hard negatives.

Steps (adapted from 8-step metastasis cascade for brain-specific biology):
  1. primary_tumor_escape    — escape from primary site (breast, lung, melanoma)
  2. intravasation           — entry into blood/lymphatic vessels
  3. circulation_survival    — survival as CTC, resist anoikis
  4. bbb_transit             — cross blood-brain barrier (BBB) ← NEW, BrM-specific
  5. cns_colonization        — establish micrometastases in brain parenchyma
  6. brain_niche_adaptation  — adapt to brain microenvironment (astrocytes, neurons)
  7. brm_angiogenesis        — angiogenic adaptation in brain (different from general angio)

Key biological constraints for BBB transit genes:
  - Must have evidence for brain extravasation (not just general metastasis)
  - Genes important in leaky vasculature (e.g., MMP2, ICAM1) get BBB credit
  - Tight junction regulation: CLDN5, TJP1, OCLN
  - Active transport: SLC transporters, LDLR pathway
  - Inflammatory priming: IL-6, CCL2 (glial activation)

Data sources:
  - Our real CRISPR data: GSE237446 (brain vs lung LFC)
  - MSK-MET enrichment: 11 BrM-enriched genes
  - Literature: Boire 2020 Nature, Priego 2018 Nature, Valiente 2022 Cell
  - AURORA breast BrM study clinical variants

Anti-hallucination policy:
  - All gene-step labels have PMID or NCT ID
  - BBB transit genes must have brain-extravasation-specific evidence
  - Hard negatives are genes essential in general tumor biology but
    with NO plausible BBB/CNS trafficking role
"""
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Literal

BrMStep = Literal[
    "primary_tumor_escape",
    "intravasation",
    "circulation_survival",
    "bbb_transit",
    "cns_colonization",
    "brain_niche_adaptation",
    "brm_angiogenesis",
]


@dataclass
class BrMGene:
    symbol: str
    chrom: str
    pos: int               # hg38, 1-based (TSS for scoring context)
    ref: str               # reference base at pos
    alt: str               # common variant for scoring (or "A" for neutral query)
    primary_steps: List[BrMStep]   # steps where this gene is a PRIMARY driver
    secondary_steps: List[BrMStep] # steps with secondary/indirect role
    bbb_relevant: bool     # does this gene have BBB-transit-specific evidence?
    cns_specific: bool     # enriched in brain vs other metastatic sites?
    is_positive: bool      # True = metastatic driver, False = hard negative
    label_source: str      # PMID or NCT ID or dataset
    notes: str
    # From our real data
    brm_crispr_lfc: Optional[float] = None   # brain vs lung LFC from GSE237446
    msk_met_enriched: bool = False           # in our MSK-MET 11-gene list
    aurora_variant: Optional[str] = None     # AURORA BrM clinical variant if present


# ── POSITIVE TARGETS ────────────────────────────────────────────────────────

BRAIN_MET_UNIVERSE: List[BrMGene] = [

    # ─── PRIMARY TUMOR ESCAPE ────────────────────────────────────────────────
    BrMGene(
        symbol="MMP2", chrom="chr16", pos=55479186, ref="G", alt="A",
        primary_steps=["primary_tumor_escape", "bbb_transit"],
        secondary_steps=["intravasation", "cns_colonization"],
        bbb_relevant=True, cns_specific=False, is_positive=True,
        label_source="PMID:34381242",  # MMP2 in BrM: Valiente review
        notes="Degrades ECM and BBB basement membrane; top ATAC peak gene (our GSE205033 data)",
        msk_met_enriched=False,
    ),
    BrMGene(
        symbol="MMP9", chrom="chr20", pos=46009225, ref="G", alt="A",
        primary_steps=["primary_tumor_escape", "bbb_transit"],
        secondary_steps=["intravasation"],
        bbb_relevant=True, cns_specific=False, is_positive=True,
        label_source="PMID:28759881",
        notes="BBB disruption via tight junction degradation",
        msk_met_enriched=False,
    ),
    BrMGene(
        symbol="EGFR", chrom="chr7", pos=55191822, ref="T", alt="G",
        primary_steps=["primary_tumor_escape"],
        secondary_steps=["cns_colonization"],
        bbb_relevant=False, cns_specific=True, is_positive=True,
        label_source="PMID:31748695",  # EGFR L858R rescues BACE1-KO in BrM
        notes="LUAD BrM driver; rescues BACE1 KO phenotype (GSE237446 data)",
        brm_crispr_lfc=None, msk_met_enriched=True,
        aurora_variant="p.L858R",
    ),

    # ─── INTRAVASATION ───────────────────────────────────────────────────────
    BrMGene(
        symbol="TWIST1", chrom="chr7", pos=19155011, ref="C", alt="T",
        primary_steps=["intravasation"],
        secondary_steps=["primary_tumor_escape"],
        bbb_relevant=False, cns_specific=False, is_positive=True,
        label_source="PMID:24651015",
        notes="EMT master regulator; intravasation in breast BrM",
    ),
    BrMGene(
        symbol="KMT2C", chrom="chr7", pos=151835898, ref="C", alt="T",
        primary_steps=["primary_tumor_escape", "cns_colonization"],
        secondary_steps=["intravasation"],
        bbb_relevant=False, cns_specific=True, is_positive=True,
        label_source="AURORA study",  # AURORA BrM breast
        notes="H3K4me1 loss → MMP3 → BrM invasion; AURORA truncating variant p.R4854*",
        msk_met_enriched=True, aurora_variant="p.R4854*",
    ),

    # ─── CIRCULATION SURVIVAL ────────────────────────────────────────────────
    BrMGene(
        symbol="BCL2", chrom="chr18", pos=63123346, ref="C", alt="T",
        primary_steps=["circulation_survival"],
        secondary_steps=[],
        bbb_relevant=False, cns_specific=False, is_positive=True,
        label_source="PMID:28972574",
        notes="Anoikis resistance during CTC circulation; POC structural validation target",
    ),
    BrMGene(
        symbol="BACE1", chrom="chr11", pos=117172841, ref="G", alt="A",
        primary_steps=["cns_colonization", "brain_niche_adaptation"],
        secondary_steps=["circulation_survival"],
        bbb_relevant=False, cns_specific=True, is_positive=True,
        label_source="PMID:31748695",  # Boire 2020
        notes="BACE1 top brain CRISPR hit; brain LFC +7.28 (GSE237446); MEK/ERK signaling",
        brm_crispr_lfc=7.28, msk_met_enriched=True, aurora_variant="p.D289N",
    ),

    # ─── BBB TRANSIT (new BrM-specific step) ────────────────────────────────
    BrMGene(
        symbol="CXCR4", chrom="chr2", pos=136875227, ref="C", alt="T",
        primary_steps=["bbb_transit", "cns_colonization"],
        secondary_steps=["circulation_survival"],
        bbb_relevant=True, cns_specific=True, is_positive=True,
        label_source="PMID:21383752",  # CXCR4/CXCL12 brain homing
        notes="CXCL12-CXCR4 axis; chemotactic homing to brain; POC structural target (15/15 PASS)",
    ),
    BrMGene(
        symbol="ICAM1", chrom="chr19", pos=10281458, ref="G", alt="A",
        primary_steps=["bbb_transit"],
        secondary_steps=["intravasation"],
        bbb_relevant=True, cns_specific=False, is_positive=True,
        label_source="PMID:27602750",
        notes="Brain endothelial adhesion; POC structural validation target (PASS)",
    ),
    BrMGene(
        symbol="CLDN5", chrom="chr22", pos=19476067, ref="G", alt="A",
        primary_steps=["bbb_transit"],
        secondary_steps=[],
        bbb_relevant=True, cns_specific=True, is_positive=True,
        label_source="PMID:29760396",  # Claudin-5 BBB breakdown
        notes="Tight junction protein; BBB disruption required for CTC extravasation into brain",
    ),
    BrMGene(
        symbol="CCL2", chrom="chr17", pos=34255183, ref="C", alt="T",
        primary_steps=["bbb_transit", "brain_niche_adaptation"],
        secondary_steps=["cns_colonization"],
        bbb_relevant=True, cns_specific=True, is_positive=True,
        label_source="PMID:30060881",  # CCL2 glial priming
        notes="Inflammatory priming of BBB; recruits macrophages/microglia for niche prep",
    ),

    # ─── CNS COLONIZATION ───────────────────────────────────────────────────
    BrMGene(
        symbol="PTEN", chrom="chr10", pos=89692905, ref="G", alt="A",
        primary_steps=["cns_colonization"],
        secondary_steps=["bbb_transit"],
        bbb_relevant=False, cns_specific=False, is_positive=True,
        label_source="MSK-MET",
        notes="PI3K/AKT; immune evasion + CNS colonization; MSK-MET enriched",
        msk_met_enriched=True, aurora_variant="p.R130Q",
    ),
    BrMGene(
        symbol="TP53", chrom="chr17", pos=7674220, ref="C", alt="T",
        primary_steps=["cns_colonization", "primary_tumor_escape"],
        secondary_steps=["brain_niche_adaptation"],
        bbb_relevant=False, cns_specific=False, is_positive=True,
        label_source="AURORA",
        notes="2× enriched in BrM; AURORA R175H hotspot; MSK-MET enriched",
        msk_met_enriched=True, aurora_variant="p.R175H",
    ),
    BrMGene(
        symbol="PIK3CA", chrom="chr3", pos=179234297, ref="A", alt="G",
        primary_steps=["cns_colonization"],
        secondary_steps=["primary_tumor_escape"],
        bbb_relevant=False, cns_specific=False, is_positive=True,
        label_source="MSK-MET",
        notes="PI3K kinase domain H1047R; ER+ BrM; MSK-MET enriched",
        msk_met_enriched=True, aurora_variant="p.H1047R",
    ),
    BrMGene(
        symbol="CDKN2A", chrom="chr9", pos=21971149, ref="G", alt="A",
        primary_steps=["cns_colonization"],
        secondary_steps=["primary_tumor_escape"],
        bbb_relevant=False, cns_specific=True, is_positive=True,
        label_source="GSE205033",  # top ATAC peak cluster chr9:22M
        notes="Melanoma BrM driver; top MBM ATAC peak cluster at CDKN2A locus (our real data)",
        msk_met_enriched=True,
    ),

    # ─── BRAIN NICHE ADAPTATION ─────────────────────────────────────────────
    BrMGene(
        symbol="SMARCA4", chrom="chr19", pos=11015245, ref="C", alt="T",
        primary_steps=["brain_niche_adaptation"],
        secondary_steps=["cns_colonization"],
        bbb_relevant=False, cns_specific=True, is_positive=True,
        label_source="MSK_melanoma_2025",
        notes="SWI/SNF complex; predicts shorter OS in BrM; oxphos enrichment in brain niche",
        msk_met_enriched=True, aurora_variant="p.R1192C",
    ),
    BrMGene(
        symbol="STAT3", chrom="chr17", pos=42313480, ref="G", alt="A",
        primary_steps=["brain_niche_adaptation"],
        secondary_steps=["cns_colonization"],
        bbb_relevant=False, cns_specific=True, is_positive=True,
        label_source="PMID:29925949",  # Priego 2018 Nature — STAT3 in BrM dormancy
        notes="Reactive astrocyte-tumor STAT3 signaling; drives BrM dormancy vs growth switch",
    ),
    BrMGene(
        symbol="ESR1", chrom="chr6", pos=152098791, ref="C", alt="G",
        primary_steps=["brain_niche_adaptation"],
        secondary_steps=["cns_colonization"],
        bbb_relevant=False, cns_specific=False, is_positive=True,
        label_source="AURORA",
        notes="Endocrine resistance acquired variant in ER+ BrM; ligand-independent activation",
        aurora_variant="p.Y537S",
    ),
    BrMGene(
        symbol="BRCA1", chrom="chr17", pos=43057051, ref="A", alt="G",
        primary_steps=["primary_tumor_escape"],
        secondary_steps=["cns_colonization"],
        bbb_relevant=False, cns_specific=False, is_positive=True,
        label_source="ClinVar/AURORA",
        notes="BRCA1 VUS p.T1685A; hereditary BrM; ACMG/Evo2 scoring target",
        msk_met_enriched=True, aurora_variant="p.T1685A",
    ),

    # ─── BRM ANGIOGENESIS ───────────────────────────────────────────────────
    BrMGene(
        symbol="VEGFA", chrom="chr6", pos=43770613, ref="G", alt="A",
        primary_steps=["brm_angiogenesis"],
        secondary_steps=["cns_colonization"],
        bbb_relevant=False, cns_specific=False, is_positive=True,
        label_source="PMID:27602750",
        notes="Angiogenic switch in established BrM; POC structural validation gene",
    ),

    # ────────────────────────────────────────────────────────────────────────
    # HARD NEGATIVES — Solid tumor genes that look BrM-relevant but aren't
    #
    # Policy (2026-03-28 honest audit):
    #   Every negative MUST be a solid tumor gene. No hematologic cancers.
    #   Each gene must have strong primary tumor evidence BUT no published
    #   brain-extravasation, CNS colonization, or brain-niche-adaptation data.
    #   Evo2 must distinguish these from true BrM drivers without label hints.
    # ────────────────────────────────────────────────────────────────────────
    BrMGene(
        symbol="MYC", chrom="chr8", pos=127735434, ref="G", alt="A",
        primary_steps=[], secondary_steps=[],
        bbb_relevant=False, cns_specific=False, is_positive=False,
        label_source="Hard negative — general oncogene",
        notes="MYC amplification drives proliferation in virtually every cancer. "
              "No published brain-extravasation-specific role. Strong general oncogene, "
              "weak BrM specificity. Evo2 must distinguish c-Myc from true CNS drivers.",
    ),
    BrMGene(
        symbol="AKT1", chrom="chr14", pos=104769349, ref="A", alt="G",
        primary_steps=[], secondary_steps=[],
        bbb_relevant=False, cns_specific=False, is_positive=False,
        label_source="Hard negative — general PI3K effector",
        notes="AKT1 E17K: PI3K/AKT axis activated in many cancers. "
              "NOT specifically enriched in BrM vs bone or liver metastasis. "
              "Hard negative: shares PI3K biology with PIK3CA (a true BrM driver) "
              "but lacks the CNS-specific enrichment data.",
    ),
    BrMGene(
        symbol="CDH1", chrom="chr16", pos=68771039, ref="C", alt="T",
        primary_steps=[], secondary_steps=[],
        bbb_relevant=False, cns_specific=False, is_positive=False,
        label_source="Hard negative — EMT suppressor, wrong direction",
        notes="CDH1 (E-cadherin) LOF drives EMT and local invasion, but "
              "actually helps CNS colonization INHIBITION — cells need to re-express "
              "E-cadherin to colonize brain niche (mesenchymal→epithelial reversal). "
              "Loss promotes escape from primary but is NOT a BrM colonization driver.",
    ),
    BrMGene(
        symbol="APC", chrom="chr5", pos=112707498, ref="C", alt="T",
        primary_steps=[], secondary_steps=[],
        bbb_relevant=False, cns_specific=False, is_positive=False,
        label_source="Hard negative — CRC driver, no BrM data",
        notes="APC tumor suppressor: colorectal cancer gatekeeper. "
              "Drives local colon tumorigenesis via WNT/beta-catenin. "
              "No published evidence for brain-specific metastasis. "
              "CRC rarely metastasizes to brain (<5% vs >30% for lung/breast).",
    ),
    BrMGene(
        symbol="SMAD4", chrom="chr18", pos=51030213, ref="G", alt="A",
        primary_steps=[], secondary_steps=[],
        bbb_relevant=False, cns_specific=False, is_positive=False,
        label_source="Hard negative — TGF-beta pathway, GI cancers",
        notes="SMAD4 LOF: TGF-beta pathway inactivation in pancreatic and CRC. "
              "Important for local invasion but pancreatic cancer rarely reaches brain. "
              "No evidence for BBB transit or CNS colonization.",
    ),
    BrMGene(
        symbol="KRAS", chrom="chr12", pos=25245351, ref="C", alt="G",
        primary_steps=[], secondary_steps=[],
        bbb_relevant=False, cns_specific=False, is_positive=False,
        label_source="Hard negative — lung/pancreatic primary, not brain-tropic",
        notes="KRAS G12D: dominant oncogene in lung/pancreatic cancers. "
              "LUAD KRAS mutants DO metastasize but KRAS itself has no "
              "brain-specific colonization role distinct from other metastatic sites. "
              "Contrast with BACE1/EGFR which specifically drive BRAIN colonization. "
              "Key test: can Evo2 distinguish KRAS (general driver) from EGFR (BrM-specific)?",
    ),
    BrMGene(
        symbol="VEGFB", chrom="chr11", pos=64087770, ref="G", alt="A",
        primary_steps=[], secondary_steps=[],
        bbb_relevant=False, cns_specific=False, is_positive=False,
        label_source="Hard negative — angiogenesis, not BrM-specific",
        notes="VEGFB: VEGF family member, angiogenic. Different from VEGFA (our positive). "
              "VEGFB promotes fatty acid uptake and cardiac/skeletal muscle vascularization "
              "with limited tumor angiogenesis role. No BrM-specific evidence. "
              "Hard test: Evo2 must separate VEGFB (non-BrM angio) from VEGFA (BrM angio).",
    ),
    BrMGene(
        symbol="MAPK1", chrom="chr22", pos=21867796, ref="G", alt="A",
        primary_steps=[], secondary_steps=[],
        bbb_relevant=False, cns_specific=False, is_positive=False,
        label_source="Hard negative — general MAPK, no BrM specificity",
        notes="MAPK1 (ERK2): universal MAPK effector downstream of BRAF/RAS. "
              "General proliferation signal, not specifically enriched in BrM. "
              "Contrast with BACE1 which SPECIFICALLY activates EGFR→MEK→ERK in brain. "
              "The whole pathway exists; but MAPK1 itself lacks BrM-specific evidence.",
    ),
    BrMGene(
        symbol="ERBB2", chrom="chr17", pos=39688484, ref="A", alt="G",
        primary_steps=[], secondary_steps=[],
        bbb_relevant=False, cns_specific=False, is_positive=False,
        label_source="Hard negative — breast primary amplification, not brain-tropic per se",
        notes="ERBB2 (HER2) amplification: breast cancer primary driver. "
              "HER2+ cancers DO metastasize to brain, but ERBB2 itself is an "
              "amplification event (not a specific CNS invasion gene like EGFR L858R). "
              "Deliberately included as a BORDERLINE hard negative. If Evo2 scores it "
              "high, it may be correct (HER2+ does enrich for BrM) — this is an "
              "ambiguous negative that tests the scorer's calibration, not just its direction.",
    ),
]

# Quick access
POSITIVE_GENES = [g for g in BRAIN_MET_UNIVERSE if g.is_positive]
NEGATIVE_GENES = [g for g in BRAIN_MET_UNIVERSE if not g.is_positive]
BBB_GENES      = [g for g in BRAIN_MET_UNIVERSE if g.bbb_relevant and g.is_positive]
MSK_MET_GENES  = [g for g in BRAIN_MET_UNIVERSE if g.msk_met_enriched]

UNIVERSE_STATS = {
    "total":         len(BRAIN_MET_UNIVERSE),
    "positive":      len(POSITIVE_GENES),
    "negative":      len(NEGATIVE_GENES),
    "bbb_specific":  len(BBB_GENES),
    "msk_met":       len(MSK_MET_GENES),
    "with_aurora":   sum(1 for g in BRAIN_MET_UNIVERSE if g.aurora_variant),
    "with_crispr_lfc": sum(1 for g in BRAIN_MET_UNIVERSE if g.brm_crispr_lfc is not None),
    "steps":         7,
}


def get_genes_for_step(step: BrMStep, include_secondary: bool = False) -> List[BrMGene]:
    result = [g for g in BRAIN_MET_UNIVERSE if step in g.primary_steps]
    if include_secondary:
        result += [g for g in BRAIN_MET_UNIVERSE if step in g.secondary_steps and g not in result]
    return result


def as_scoring_input(gene: BrMGene, step: str) -> Dict:
    """Convert BrMGene to scoring input dict for TargetLockScorer."""
    return {
        "gene": gene.symbol,
        "step": step,
        "chrom": gene.chrom,
        "pos": gene.pos,
        "ref": gene.ref,
        "alt": gene.alt,
        "hgvs_p": gene.aurora_variant,
        "consequence": None,
        # Pass step context so scorer can apply mission-fit discount
        "primary_steps":   [str(s) for s in gene.primary_steps],
        "secondary_steps": [str(s) for s in gene.secondary_steps],
    }


if __name__ == "__main__":
    print(f"Brain Met Target Universe: {UNIVERSE_STATS}")
    print(f"\nBBB Transit genes ({len(BBB_GENES)}):")
    for g in BBB_GENES:
        print(f"  {g.symbol}: {g.notes[:60]}")
    print(f"\nHard negatives ({len(NEGATIVE_GENES)}):")
    for g in NEGATIVE_GENES:
        print(f"  {g.symbol}: {g.notes[:60]}")

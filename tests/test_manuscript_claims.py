"""Assert every numeric claim in the manuscript against the results it describes.

Three rounds of manual review found stale numbers in the manuscript that the
analysis had moved underneath: "32,052 ortholog pairs" against 31,993, "24
MAGMA-significant genes" against 25, "16 cross-species concordant transcripts"
against 2, "88 CIPN features" against 87. Every one of them would have been
caught here in a single run, and would have stayed caught as the analysis
changed.

Each claim pairs a regex locating the sentence in the .tex with a function
computing what the number should be from the tracked outputs. A claim whose
pattern no longer matches fails loudly: if a sentence is rewritten, its entry
belongs in this file too, and silence would defeat the purpose.

results/ is gitignored and regenerable, so claims whose inputs are absent skip
rather than fail. A fresh clone reports skips, not a broken suite; run the
pipeline first to get real coverage.
"""

from __future__ import annotations

import math
import re
from collections.abc import Callable
from dataclasses import dataclass, field
from functools import cache
from pathlib import Path

import pandas as pd
import pytest
from scipy import stats

REPO_ROOT = Path(__file__).resolve().parent.parent
MANUSCRIPT = REPO_ROOT / "manuscript" / "manuscript.tex"
ABBREVIATED_SUMMARY = REPO_ROOT / "manuscript" / "abbreviated_summary.txt"
COVER_LETTER = REPO_ROOT / "manuscript" / "cover_letter.md"
PRISMA_CHECKLIST = REPO_ROOT / "manuscript" / "prisma_2020_checklist.md"
SUPPLEMENTARY = REPO_ROOT / "manuscript" / "supplementary_body.tex"

#: Documents whose numbers this registry checks, keyed by the name a claim
#: names in `source`. The mutation harness enumerates this map too, so a
#: document added here is mutation-tested without being named a second time.
#:
#: The Elsevier Highlights left with the Journal of Pain on 2026-09-21. Brain
#: Communications asks instead for a 50-word abbreviated summary, which states
#: no number (test_abbreviated_summary_meets_the_journal_format keeps it that
#: way), so it needs no entry here.
DOCUMENTS: dict[str, Path] = {
    "manuscript": MANUSCRIPT,
    # \input by supplementary.tex, the standalone document uploaded beside the
    # paper; since the move to Brain Communications it carries the detailed
    # methods and results the 6,000-word limit sent out of the paper. It is its
    # own entry because the mutation harness perturbs files: a number it could
    # not find in manuscript.tex would be reported unperturbable and checked by
    # nothing.
    "supplementary": SUPPLEMENTARY,
    # The Brain Communications cover letter restates headline numbers, and an
    # editor reads it first. Markdown, not LaTeX: "%" is literal here.
    "cover_letter": COVER_LETTER,
}

TX_META = REPO_ROOT / "results" / "meta" / "transcriptomics"
BY_SPECIES = TX_META / "by_species"
STRAT = TX_META / "stratified"
CROSS = REPO_ROOT / "results" / "cross_species" / "transcriptomics"
GENOMICS = REPO_ROOT / "results" / "meta" / "genomics"
PROTEOMICS = REPO_ROOT / "results" / "per_study" / "proteomics" / "PXD013362"
PROT_PER_STUDY = REPO_ROOT / "results" / "per_study" / "proteomics"
PROT_META = REPO_ROOT / "results" / "meta" / "proteomics"
PRISMA = REPO_ROOT / "literature" / "prisma"
PEPTIDE_MAP = REPO_ROOT / "conf" / "analysis" / "peptide_gene_map.csv"
MET_MAP = REPO_ROOT / "conf" / "analysis" / "metabolite_id_map.csv"
MET_CROSSCHECK = (REPO_ROOT / "results" / "meta" / "metabolomics"
                  / "corpus_chebi_crosscheck.csv")
MET_META = REPO_ROOT / "results" / "meta" / "metabolomics"
MET_PER_STUDY = REPO_ROOT / "results" / "per_study" / "metabolomics"


class MissingInput(Exception):
    """Raised when a claim's underlying results file has not been generated."""


def _read(path: Path, **kwargs) -> pd.DataFrame:
    if not path.exists():
        raise MissingInput(str(path.relative_to(REPO_ROOT)))
    return pd.read_csv(path, **kwargs)


@cache
def _pooled() -> pd.DataFrame:
    return _read(TX_META / "pooled_effects.csv",
                 usecols=["feature_id", "k", "padj_pooled"])


@cache
def _species(slug: str) -> pd.DataFrame:
    return _read(BY_SPECIES / f"{slug}_pooled.csv", usecols=["feature_id", "k", "padj"])


@cache
def _stratum(name: str) -> pd.DataFrame:
    return _read(STRAT / f"{name}_pooled.csv", usecols=["feature_id", "k", "padj"])


@cache
def _stratum_effects(name: str) -> pd.DataFrame:
    """A stratum with its pooled effect sizes, for claims that quote a log2FC."""
    return _read(STRAT / f"{name}_pooled.csv", usecols=["feature_id", "yi", "padj"])


@cache
def _prot_pooled() -> pd.DataFrame:
    return _read(PROT_META / "pooled_effects.csv",
                 usecols=["feature_id", "k", "yi_pooled", "I2", "padj_pooled"])


def _prot_family() -> pd.DataFrame:
    """The features BH correction was applied over: k >= meta.min_studies."""
    return _prot_pooled()[_prot_pooled()["padj_pooled"].notna()]


def _prot_gene(gene: str, col: str) -> float:
    row = _prot_pooled()[_prot_pooled()["feature_id"] == gene]
    if row.empty:
        raise MissingInput(f"{gene} absent from pooled proteomics")
    return float(row[col].iloc[0])


@cache
def _prot_study(acc: str) -> pd.DataFrame:
    return _read(PROT_PER_STUDY / acc / f"{acc}_effects.csv",
                 usecols=["feature_id", "padj"])


@cache
def _snl_stats() -> dict:
    """Power diagnostics for the SNL stratum, read from its pooled table."""
    d = _read(STRAT / "SNL_pooled.csv",
              usecols=["feature_id", "k", "yi", "I2", "padj"])
    fam = d[d["padj"].notna()]
    sig = fam[fam["padj"] < 0.05]
    return {
        "pct_sig": 100.0 * len(sig) / len(fam),
        "median_abs_g_sig": float(sig["yi"].abs().median()),
        "median_abs_g_all": float(fam["yi"].abs().median()),
        "median_i2_sig": float(sig["I2"].median()),
    }


@cache
def _triage_metab() -> pd.DataFrame:
    return _read(PRISMA / "metabolomics_file_triage.csv",
                 usecols=["accession", "verdict", "reason_code"])


@cache
def _metab_pooled() -> pd.DataFrame:
    return _read(MET_META / "pooled_effects.csv",
                 usecols=["feature_id", "k", "padj_pooled"])


@cache
def _metab_study(acc: str) -> pd.DataFrame:
    return _read(MET_PER_STUDY / acc / f"{acc}_effects.csv", usecols=["se"])


@cache
def _metab_units() -> tuple[str, ...]:
    """Study units actually entering the metabolomic pool.

    Derived from disk and filtered through the superseded registry rather than
    listed here. A hardcoded list was wrong the moment the arm grew: it still
    named the three MetaboLights studies after six Workbench studies had
    joined, so the weight-share claim computed a three-study denominator and
    reported 90% where the truth was 38% — the claim would have passed while
    testing a cohort that no longer exists.
    """
    # `iterdir` on a missing directory raises FileNotFoundError, which is not
    # MissingInput and so does not reach the skip path: a fresh clone with no
    # results/ errored the whole test instead of skipping, which CI caught.
    if not MET_PER_STUDY.is_dir():
        raise MissingInput(str(MET_PER_STUDY.relative_to(REPO_ROOT)))
    superseded = set(_read(REPO_ROOT / "conf" / "analysis" / "superseded_studies.csv",
                           usecols=["study_id"])["study_id"])
    units = sorted(d.name for d in MET_PER_STUDY.iterdir()
                   if d.is_dir() and (d / f"{d.name}_effects.csv").exists())
    if not units:
        raise MissingInput(f"no per-study effect tables under {MET_PER_STUDY.name}")
    return tuple(u for u in units if u not in superseded)


@cache
def _metab_padj(feature_id: str) -> float:
    """Adjusted p-value for one pooled metabolite."""
    pooled = _metab_pooled()
    row = pooled[pooled["feature_id"] == feature_id]
    assert len(row) == 1, f"{feature_id}: {len(row)} rows in the pooled table"
    return float(row["padj_pooled"].iloc[0])


def _metab_weight_share(acc: str) -> float:
    """Percent of fixed-effect weight, approximated at each study's median SE.

    Exact per-feature shares vary; the manuscript quotes the round figure this
    approximates, so the claim carries a tolerance rather than demanding the
    precision the sentence does not use.
    """
    w = {a: 1.0 / float(_metab_study(a)["se"].median()) ** 2 for a in _metab_units()}
    return 100.0 * w[acc] / sum(w.values())


@cache
def _triage() -> pd.DataFrame:
    return _read(PRISMA / "proteomics_triage.csv", usecols=["accession", "verdict"])


@cache
def _concordance(slug: str) -> pd.DataFrame:
    return _read(CROSS / f"concordance_{slug}.csv")


@cache
def _cs_summary() -> pd.DataFrame:
    return _read(CROSS / "summary.csv")


def _summary_row(species: str) -> pd.Series:
    df = _cs_summary()
    row = df[df["species"] == species]
    if row.empty:
        raise MissingInput(f"summary.csv has no row for {species}")
    return row.iloc[0]


def _sens_binomial_p(arm: str, species: str) -> float:
    """The binomial p Table 4 prints for one sensitivity arm."""
    row = _sens_row(arm, species)
    return float(stats.binomtest(int(row["n_concordant"]),
                                 int(row["n_ortholog_pairs"]), 0.5).pvalue)


def _binomial_p(species: str) -> float:
    """Two-sided binomial test of concordance against 50%.

    Recomputed rather than read: summary.csv carries the correlation p-values,
    not this one. scripts/build_concordance_table.py computes Table 4's column
    the same way, so the prose and the table cannot disagree about it.
    """
    row = _summary_row(species)
    return float(stats.binomtest(int(row["n_concordant"]),
                                 int(row["n_ortholog_pairs"]), 0.5).pvalue)


@cache
def _median_k(which: str) -> int:
    """Median contributing-study count, over the BH family or its significant set.

    The whole-cohort significant set is now pooled from far more studies than
    the family as a whole; the Results sentence says so, so both numbers are
    read from the file rather than asserted.
    """
    d = _read(TX_META / "pooled_effects.csv", usecols=["k", "padj_pooled"])
    fam = d[d["padj_pooled"].notna()]
    sub = fam[fam["padj_pooled"] < 0.05] if which == "significant" else fam
    return int(sub["k"].median())


@cache
def _bonferroni_survivors(stratum: str) -> int:
    """Features in a stratum whose raw p clears Bonferroni over its own family.

    The manuscript contrasts the BH count against this. It has to be derived
    from the raw p-values: `padj` is a BH quantity and comparing it to a
    Bonferroni threshold would be a category error.
    """
    d = _read(STRAT / f"{stratum}_pooled.csv", usecols=["pval", "padj"])
    fam = d[d["padj"].notna()]
    return int((fam["pval"] < 0.05 / len(fam)).sum())


@cache
def _species_studies(slug: str) -> int:
    """Distinct studies contributing to a species pool."""
    d = _read(BY_SPECIES / f"{slug}_pooled.csv", usecols=["study_ids"])
    return d["study_ids"].astype(str).str.split(";").explode().str.strip().nunique()


#: Identifier forms that are not gene names: repository placeholders and
#: array/accession addresses. The abstract quotes the complement of this as
#: "named genes", and Methods states the rule, so the number is reproducible
#: rather than a judgement. Kept separate from _UNMAPPED_ID, which asks a
#: different question (was this study mapped at all).
_UNNAMED_FEATURE = re.compile(
    r"^(LOC\d+"                  # unnamed NCBI loci
    r"|Gm\d+"                    # predicted mouse genes
    r"|ENS[A-Z]*\d{6,}"          # Ensembl accessions
    r"|AABR\d|RGD\d+"           # rat clone/RGD identifiers
    r"|ASHGV\d+|A_\d+_P\d+"     # array addresses
    r"|\(\+\)|\(-\)"           # array spike-in controls
    r"|[A-Z]{2}\d{6}\.\d"        # clone-based accessions
    r"|\d+_[sx]?_?at|\d+_st"      # Affymetrix probeset addresses
    r"|.*Rik\d*l?$"              # Riken clones
    r"|.*-ps\d*$"                # processed pseudogenes
    r"|Mir\d|MIR\d|Snor|LINC\d)",  # unnamed RNA classes
    re.IGNORECASE)


DEPTH = TX_META / "depth_sensitivity.csv"
DEPTH_SHARED = TX_META / "depth_sensitivity_shared.csv"


@cache
def _depth(species: str, shared: bool = False) -> dict:
    """Subsampling results for one species, from 06h_depth_sensitivity.R.

    None of these numbers was claimed before 2026-08-30, and 06h had not been
    re-run since 2026-08-13, so the whole depth paragraph was verified against
    a corpus that had stopped existing -- the rat correction family it
    described was 47,950 against 24,380 today, and the median-k asymmetry the
    argument rested on (13 v 3) had become 15 v 16.
    """
    d = _read(DEPTH_SHARED if shared else DEPTH)
    s = d[d["species"] == species]
    sub = s[s["seed"].notna()]
    out = {"mean": float(sub["pct_sig"].mean()),
           "lo": float(sub["pct_sig"].min()),
           "hi": float(sub["pct_sig"].max()),
           "pct_subsampled": 100.0 * float(sub["n_subsampled"].iloc[0])
                             / float(sub["n_features"].iloc[0])}
    obs = s[s["seed"].isna()]
    if len(obs):
        out["observed"] = float(obs["pct_sig"].iloc[0])
    return out


@cache
def _median_per_study_se(species: str) -> float:
    """Median SE across a species' per-study tables.

    Pooling depth cannot affect this, which is why the paragraph uses it to
    separate precision from depth.
    """
    assign = _read(BY_SPECIES / "species_assignment.csv")
    studies = assign.loc[assign["species"] == species, "study_id"]
    frames = []
    per_study = REPO_ROOT / "results" / "per_study" / "transcriptomics"
    for study in studies:
        for suffix in ("_effects_normalized.csv", "_effects.csv"):
            path = per_study / f"{study}{suffix}"
            if path.exists():
                frames.append(_read(path, usecols=["se"])["se"])
                break
    if not frames:
        raise MissingInput(f"no per-study tables for {species}")
    se = pd.concat(frames).dropna()
    return float(se[se > 0].median())


@cache
def _species_median_k(slug: str) -> int:
    d = _read(BY_SPECIES / f"{slug}_pooled.csv", usecols=["k", "padj"])
    return int(d[d["padj"].notna()]["k"].median())


@cache
def _species_sig_pct(slug: str) -> float:
    """A species pool's significant share of its own correction family.

    The paragraph quotes this for all three species and it was claimed for
    none of them; on 2026-08-31 the Results gave rat's as 5.27% eight lines
    after giving it as 6.26%, which is the drift the registry exists to stop.
    """
    d = _species(slug)
    fam = d[d["padj"].notna()]
    if fam.empty:
        raise MissingInput(f"{slug} pool has an empty correction family")
    return 100.0 * float((fam["padj"] < 0.05).sum()) / len(fam)


#: Median tau-squared below this counts as "zero" for the depth argument. The
#: pools sit four orders of magnitude apart around it: 1e-6 at k=3 in both
#: species against 3e-3 for rat at k=4.
_TAU2_RISEN = 1e-5


@cache
def _tau2_rise_k(slug: str) -> int:
    """First pooling depth at which a species' median tau-squared has risen.

    The sentence reads this off a curve, which is exactly the kind of number
    that goes stale unread: mouse's rise moved from beyond k=20 to k=16 when
    the corpus grew to 82 units, and nothing would have said so.
    """
    d = _read(BY_SPECIES / f"{slug}_pooled.csv", usecols=["k", "tau2"])
    med = d[d["k"] >= 3].groupby("k")["tau2"].median()
    risen = med[med >= _TAU2_RISEN]
    if risen.empty:
        raise MissingInput(f"{slug} median tau2 never reaches {_TAU2_RISEN}")
    return int(risen.index.min())


#: The five strata of the original design that reach no corrected result. The
#: sentence quoting their feature range names them, so the list is explicit
#: rather than derived -- a stratum that started pooling should break the
#: sentence, not silently leave the range.
_ORIGINAL_EMPTY_STRATA = ("CFA", "neuropathic_human", "lbp_human",
                          "nociplastic_human", "invitro_human")


@cache
def _all_strata() -> list[str]:
    """Every stratum 06b wrote.

    Raises rather than returning [] when there are none. A glob over a missing
    directory is empty, not an error, so without this the claims counting
    strata reported "manuscript says 17, results give 0" on a checkout with no
    results/ -- a red suite where the contract is a skip. Every other helper
    gets this free from _read.
    """
    found = sorted(p.name.replace("_pooled.csv", "")
                   for p in STRAT.glob("*_pooled.csv"))
    if not found:
        raise MissingInput(str(STRAT.relative_to(REPO_ROOT)) + "/*_pooled.csv")
    return found


@cache
def _stratum_merged_arms(name: str) -> int:
    """Studies in a stratum whose case or control arm merges several levels.

    Reads the contrast audit sheet, so it tracks a re-run of
    scripts/audit_contrasts.py rather than a hand count.
    """
    pooled = _studies_in(STRAT / f"{name}_pooled.csv")
    audit = _contrast_audit()
    rows = audit[audit["accession"].isin(pooled)]
    return int(rows["flags"].fillna("").str.contains("arms_merge").sum())


@cache
def _sni_tissue() -> pd.DataFrame:
    return _read(STRAT / "sni_tissue_meta.csv", usecols=["feature_id", "padj"])


@cache
def _sni_tissue_top_n() -> int:
    """How many candidates 06d ranks down to before fitting."""
    src = (REPO_ROOT / "pipeline" / "06d_sni_tissue_stratified.R").read_text()
    m = re.search(r"head\(sni_pooled_idx\$feature_id,\s*(\d+)\)", src)
    if not m:
        raise MissingInput("06d_sni_tissue_stratified.R: no head(...) cap found")
    return int(m.group(1))


@cache
def _q_permutation() -> pd.DataFrame:
    """Permutation null for the between-model Q test (100 shuffles).

    perm 0 is the observed statistic; perm > 0 the permuted ones. Written by
    scripts/permute_between_model_q.R.
    """
    return _read(STRAT / "between_model_Q_permutation.csv")


@cache
def _hvh(design: str) -> pd.DataFrame:
    """One design of the human-vs-human positive control.

    scripts/human_vs_human_concordance.py. "replicate" is the tightest arm:
    two halves of the endometriosis stratum, so species, compartment and
    condition are all held constant.
    """
    d = _read(CROSS / "human_vs_human.csv")
    return d[d["design"] == design]


def _q_excess_over_null(nearest: int) -> float:
    """Features differing by pain model beyond what label-shuffling produces.

    The permutation runs on a subsample, so it yields a rate rather than a
    count: the excess is (observed - mean permuted) share applied to the whole
    testable set. Returned rounded to `nearest`, because the number is quoted
    rounded ("some 8,600") and a claim must expect what the sentence says.
    Binding the raw 8,573 with a tolerance wide enough to admit 8,600 would
    admit 8,601 as well, and the mutation harness perturbs by one -- the claim
    would survive its own mutation and check nothing.
    """
    perm = _q_permutation()
    observed = float(perm.query("perm == 0")["pct"].iloc[0])
    null = float(perm.query("perm > 0")["pct"].mean())
    testable = len(_read(STRAT / "between_model_Q.csv", usecols=["padj_QM"]))
    return round((observed - null) / 100.0 * testable / nearest) * nearest


def _q_null_attributable(nearest: int) -> float:
    """Model-dependent features the permuted null already accounts for.

    The complement of `_q_excess_over_null` within the significant set, and
    rounded the same way and for the same reason. The Results said "around
    5,000 of the 14,267" -- a count from two regenerations earlier, restating
    a figure the same paragraph already gives, and read by no claim.
    """
    q = _read(STRAT / "between_model_Q.csv", usecols=["padj_QM"])
    significant = int((q["padj_QM"] < 0.05).sum())
    return round((significant - _q_excess_over_null(1)) / nearest) * nearest


@cache
def _gwas_replication() -> pd.DataFrame:
    """Lead-variant genes ranked by how many independent GWAS report them."""
    return _read(GENOMICS / "gwas_gene_meta_replication.csv",
                 usecols=["gene", "n_studies", "min_p"])


@cache
def _gwas_fig_effects() -> pd.DataFrame:
    """Pooled transcriptomic estimates for the genes the figure shows.

    Selection mirrors pipeline/15_gwas_tx_overlap.py: the top 15 by
    replication count, ties broken on the minimum p-value.
    """
    rep = _gwas_replication().sort_values(["n_studies", "min_p"],
                                          ascending=[False, True]).head(15)
    pooled = _read(TX_META / "pooled_effects.csv",
                   usecols=["feature_id", "k", "yi_pooled", "padj_pooled"])
    out = pooled[pooled["feature_id"].isin(rep["gene"])]
    if len(out) != len(rep):
        raise MissingInput(
            f"{len(rep) - len(out)} of the figure's genes are absent from the "
            "pooled table"
        )
    return out


@cache
def _gwas_fetched() -> int:
    """GWAS whose harmonized summary statistics were retrieved and analysed.

    The manifest records all 138 retrievable studies with a disposition; the
    analysed set is the fetched one. Reading the manifest rather than a
    literal keeps the count honest if a study is later excluded.
    """
    d = _read(GENOMICS / "sumstats_manifest.csv", usecols=["disposition"])
    return int((d["disposition"] == "fetched").sum())


@cache
def _between_model_q() -> pd.DataFrame:
    return _read(STRAT / "between_model_Q.csv",
                 usecols=["feature_id", "QM_pval", "padj_QM"])


@cache
def _contrast_audit() -> pd.DataFrame:
    return _read(PRISMA / "transcriptomics_contrast_audit.csv",
                 usecols=["accession", "separating_fields", "flags"])


@cache
def _egger() -> pd.DataFrame:
    return _read(TX_META / "egger_test.csv", usecols=["k", "bias_flagged"])


@cache
def _egger_z() -> pd.DataFrame:
    """Egger z for the flagged features only -- the set the sentence describes."""
    d = _read(TX_META / "egger_test.csv", usecols=["egger_z", "bias_flagged"])
    return d[d["bias_flagged"].astype(bool)]


@cache
def _effect_estimate_count() -> int:
    """Rows in the normalized per-study table: one effect estimate each."""
    d = _read(REPO_ROOT / "results" / "per_study" / "transcriptomics"
              / "effects_normalized_all.csv", usecols=["study_id"])
    return len(d)


@cache
def _pooled_study_count() -> int:
    """Study units contributing to the whole-cohort pool."""
    d = _read(TX_META / "pooled_effects.csv", usecols=["study_ids"])
    return d["study_ids"].astype(str).str.split(";").explode().str.strip().nunique()


@cache
def _named_gene_fraction() -> float:
    """Percentage of pooled features whose identifier is a gene name."""
    ids = _read(TX_META / "pooled_effects.csv",
                usecols=["feature_id"], dtype=str)["feature_id"].astype(str)
    return 100.0 * float((~ids.str.match(_UNNAMED_FEATURE)).mean())


@cache
def _search_log_last() -> pd.Series:
    """The most recent literature-search run, as PRISMA reports it."""
    return _read(PRISMA / "search_log.csv").iloc[-1]


@cache
def _literature_included() -> int:
    """Records the literature screen retained -- the figure's final box."""
    d = _read(PRISMA / "screening.csv", usecols=["screening_decision"])
    return int((d["screening_decision"] == "include").sum())


@cache
def _prisma_verdicts() -> pd.Series:
    return (_read(PRISMA / "transcriptomics_candidates.csv", dtype=str)
            ["verdict"].value_counts())


SC_META = REPO_ROOT / "results" / "meta" / "single_cell"
SC_MANIFEST = (REPO_ROOT / "data" / "interim" / "single_cell"
               / "pseudobulk_manifest.csv")


@cache
def _sc_pooled() -> pd.DataFrame:
    """The pseudobulk single-cell pool.

    Its own stratum, never merged into the bulk transcriptomic pool: the
    comparison worth making is bulk against pseudobulk on the same genes.
    """
    return _read(SC_META / "pooled_effects.csv",
                 usecols=["feature_id", "k", "study_ids", "padj_pooled"])


@cache
def _sc_manifest() -> pd.DataFrame:
    """One row per pseudobulked sample, with the cells each contributed."""
    return _read(SC_MANIFEST, usecols=["subject", "cells_kept"])


@cache
def _sc_studies() -> int:
    """Study units contributing to the single-cell pool."""
    return (_sc_pooled()["study_ids"].astype(str).str.split(";")
            .explode().str.strip().nunique())


@cache
def _top_q_i2() -> pd.DataFrame:
    """The five features leading the between-model Q test, with their I2.

    The sentence naming them also characterised their heterogeneity, and that
    half was read by nothing: it said all five exceeded 99% while the feature
    it names first, ARHGAP36, sits at 15.
    """
    d = _read(STRAT / "between_model_Q.csv",
              usecols=["feature_id", "padj_QM", "QM_pval", "I2_total"])
    return d.sort_values(["padj_QM", "QM_pval"]).head(5)


@cache
def _endometrial_human_units() -> int:
    """Human pool units whose tissue compartment is endometrial.

    Read from conf/analysis/human_tissue_compartments.csv rather than from the
    endometriosis stratum, because the two answered differently while seven
    rescued endometriosis studies sat in no stratum at all.
    """
    comp = _read(REPO_ROOT / "conf" / "analysis" / "human_tissue_compartments.csv",
                 dtype=str).set_index("study_id")["compartment"]
    d = _read(BY_SPECIES / "Homo_sapiens_pooled.csv", usecols=["study_ids"])
    units = (d["study_ids"].astype(str).str.split(";").explode()
             .str.strip().unique())
    return sum(1 for u in units
               if comp.get(u.split("_")[0]) == "endometrial")


@cache
def _amendment_admitted() -> int:
    """Studies the 2026-08-28 repository amendment brought in and kept.

    The Methods disclose the amendment's composition, and the sentence was
    written when it stood at 24. Nothing read either figure, so it survived the
    rescue that took the cohort to 34.
    """
    d = _read(PRISMA / "transcriptomics_candidates.csv", dtype=str)
    kept = d["verdict"] == "include"
    amended = d["retrieval"].astype(str).str.contains("amendment", na=False)
    return int((kept & amended).sum())


@cache
def _table_s1() -> pd.DataFrame:
    return _read(REPO_ROOT / "manuscript" / "supplementary"
                 / "Table_S1_transcriptomic_studies.csv")


def _n_sig(df: pd.DataFrame, col: str = "padj") -> int:
    fam = df[df[col].notna()]
    return int((fam[col] < 0.05).sum())


def _concordant_transcripts() -> set[str]:
    genes: set[str] = set()
    for slug in ("musculus", "norvegicus"):
        d = _concordance(slug)
        genes |= set(d.loc[d["both_significant"] & d["concordant"], "human_feature_id"])
    return genes


# ---------------------------------------------------------------------------
# Claim registry
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Claim:
    """One numeric assertion made by the manuscript.

    `pattern` must contain exactly one capture group holding the number as the
    LaTeX writes it, and must match exactly once in the document, so that a
    claim cannot be silently duplicated into two disagreeing sentences.

    `source` names the document the pattern is matched against, one of the
    keys of DOCUMENTS. It is per-claim rather than per-registry because the
    "matches exactly once" rule is only meaningful within one document: a
    document submitted beside the paper may restate one of its figures by
    design (the Highlights did, for the Journal of Pain), and counting that
    restatement as a duplicate would forbid the document from existing.
    """

    name: str
    pattern: str
    expected: Callable[[], float]
    tol: float = 0.0          # absolute tolerance; 0 demands exactness
    note: str = field(default="")
    source: str = "manuscript"




SENS = CROSS / "sensitivity"


@cache
def _sens_summary(arm: str) -> pd.DataFrame:
    """One arm of the cross-species sensitivity analysis (Table 4).

    The arms differ only in which human studies are pooled; the rodent pools
    are identical across all three, which is what makes the comparison a
    statement about compartment rather than about depth.
    """
    # "all_human" is the main analysis, not a variant of it: the sensitivity
    # arms exist to say what the published pool and a tissue-matched pool give
    # against the same rodent pools. Copying it under sensitivity/ would leave
    # two files that must agree and eventually will not.
    path = CROSS / "summary.csv" if arm == "all_human" else SENS / arm / "summary.csv"
    return _read(path)


def _sens_row(arm: str, species: str) -> pd.Series:
    d = _sens_summary(arm)
    return d[d["species"] == species].iloc[0]




@cache
def _stratum_studies(name: str) -> int:
    """Distinct study units contributing to a stratum.

    Table 2 prints Studies, Features and Pooled beside Sig, and until 2026-08-29
    only Sig was claimed for every stratum (and Pooled for SNL alone). That is
    how SNL's Features column kept its pre-2026-08-16 value of 287,084 through
    two green claim runs: nothing read it.
    """
    d = _read(STRAT / f"{name}_pooled.csv", usecols=["study_ids"])
    return d["study_ids"].astype(str).str.split(";").explode().str.strip().nunique()


def _sens_claims() -> list[Claim]:
    """Table 4's cells: three arms x two species x three columns.

    Built in a loop rather than written out. The eighteen differ only in the
    arm label and the species, and in the hand-written version the arm's study
    count was baked into each pattern -- so when the human pool went from 24
    studies to 23 on 2026-08-30, six claims stopped matching at once and
    reported "no longer appears in the manuscript" rather than a wrong number.
    The count is now a claim of its own, and the row patterns do not care.
    """
    arms = [("published", r"Published \(5\)", "published_pool"),
            ("tissue", r"Tissue-matched \(10\)", "tissue_matched"),
            ("all", r"All human \(\d+\)", "all_human")]
    species = [("mouse", "Mus musculus"), ("rat", "Rattus norvegicus")]

    out: list[Claim] = []
    for arm_tag, label, slug in arms:
        for sp_tag, sp in species:
            row = label + r" & \\emph\{" + sp + r"\}"
            out += [
                Claim(f"sens_{arm_tag}_{sp_tag}_pairs",
                      row + r" & ([\d,]+) &",
                      lambda slug=slug, sp=sp:
                          int(_sens_row(slug, sp)["n_ortholog_pairs"])),
                Claim(f"sens_{arm_tag}_{sp_tag}_pct",
                      row + r" & [\d,]+ & ([\d.]+)\\%",
                      lambda slug=slug, sp=sp:
                          float(_sens_row(slug, sp)["pct_concordant"]),
                      tol=0.05),
                Claim(f"sens_{arm_tag}_{sp_tag}_rho",
                      row + r" & [\d,]+ & [\d.]+\\% & [^&]+ & (?:\$-\$)?([\d.]+) ",
                      lambda slug=slug, sp=sp:
                          abs(float(_sens_row(slug, sp)["spearman_r"])),
                      tol=0.0005),
            ]
    # The whole-corpus arm's study count, which the row patterns deliberately
    # no longer pin. It moves with the human pool: 24 before GSE262037 was
    # held, 23 after.
    out.append(Claim("sens_all_human_studies",
                     r"All human \((\d+)\) & \\emph\{Mus musculus\}",
                     lambda: _species_studies("Homo_sapiens")))
    return out


CLAIMS: list[Claim] = [
    # -- stratified table (Table 2), all columns ----------------------------
    Claim("strat_cci_studies",
          r"CCI & (\d+) & [\d,]+ & [\d,]+ &",
          lambda slug="CCI": _stratum_studies(slug)),
    Claim("strat_cci_features",
          r"CCI & \d+ & ([\d,]+) & [\d,]+ &",
          lambda slug="CCI": len(_stratum(slug))),
    Claim("strat_cci_pooled",
          r"CCI & \d+ & [\d,]+ & ([\d,]+) &",
          lambda slug="CCI": int(_stratum(slug)["padj"].notna().sum())),
    Claim("strat_sni_studies",
          r"SNI & (\d+) & [\d,]+ & [\d,]+ &",
          lambda slug="SNI": _stratum_studies(slug)),
    Claim("strat_sni_features",
          r"SNI & \d+ & ([\d,]+) & [\d,]+ &",
          lambda slug="SNI": len(_stratum(slug))),
    Claim("strat_sni_pooled",
          r"SNI & \d+ & [\d,]+ & ([\d,]+) &",
          lambda slug="SNI": int(_stratum(slug)["padj"].notna().sum())),
    Claim("strat_snl_studies",
          r"SNL & (\d+) & [\d,]+ & [\d,]+ &",
          lambda slug="SNL": _stratum_studies(slug)),
    Claim("strat_snl_features",
          r"SNL & \d+ & ([\d,]+) & [\d,]+ &",
          lambda slug="SNL": len(_stratum(slug))),
    Claim("strat_snl_pooled",
          r"SNL & \d+ & [\d,]+ & ([\d,]+) &",
          lambda slug="SNL": int(_stratum(slug)["padj"].notna().sum())),
    Claim("strat_cfa_studies",
          r"CFA & (\d+) & [\d,]+ & [\d,]+ &",
          lambda slug="CFA": _stratum_studies(slug)),
    Claim("strat_cfa_features",
          r"CFA & \d+ & ([\d,]+) & [\d,]+ &",
          lambda slug="CFA": len(_stratum(slug))),
    Claim("strat_cfa_pooled",
          r"CFA & \d+ & [\d,]+ & ([\d,]+) &",
          lambda slug="CFA": int(_stratum(slug)["padj"].notna().sum())),
    Claim("strat_cipn_studies",
          r"CIPN & (\d+) & [\d,]+ & [\d,]+ &",
          lambda slug="CIPN": _stratum_studies(slug)),
    Claim("strat_cipn_features",
          r"CIPN & \d+ & ([\d,]+) & [\d,]+ &",
          lambda slug="CIPN": len(_stratum(slug))),
    Claim("strat_cipn_pooled",
          r"CIPN & \d+ & [\d,]+ & ([\d,]+) &",
          lambda slug="CIPN": int(_stratum(slug)["padj"].notna().sum())),
    Claim("strat_neuro_human_studies",
          r"Neuropathic \(human\) & (\d+) & [\d,]+ & [\d,]+ &",
          lambda slug="neuropathic_human": _stratum_studies(slug)),
    Claim("strat_neuro_human_features",
          r"Neuropathic \(human\) & \d+ & ([\d,]+) & [\d,]+ &",
          lambda slug="neuropathic_human": len(_stratum(slug))),
    Claim("strat_neuro_human_pooled",
          r"Neuropathic \(human\) & \d+ & [\d,]+ & ([\d,]+) &",
          lambda slug="neuropathic_human": int(_stratum(slug)["padj"].notna().sum())),
    Claim("strat_lbp_studies",
          r"LBP \(human\) & (\d+) & [\d,]+ & [\d,]+ &",
          lambda slug="lbp_human": _stratum_studies(slug)),
    Claim("strat_lbp_features",
          r"LBP \(human\) & \d+ & ([\d,]+) & [\d,]+ &",
          lambda slug="lbp_human": len(_stratum(slug))),
    Claim("strat_lbp_pooled",
          r"LBP \(human\) & \d+ & [\d,]+ & ([\d,]+) &",
          lambda slug="lbp_human": int(_stratum(slug)["padj"].notna().sum())),
    Claim("strat_nocip_studies",
          r"Nociplastic \(human\) & (\d+) & [\d,]+ & [\d,]+ &",
          lambda slug="nociplastic_human": _stratum_studies(slug)),
    Claim("strat_nocip_features",
          r"Nociplastic \(human\) & \d+ & ([\d,]+) & [\d,]+ &",
          lambda slug="nociplastic_human": len(_stratum(slug))),
    Claim("strat_nocip_pooled",
          r"Nociplastic \(human\) & \d+ & [\d,]+ & ([\d,]+) &",
          lambda slug="nociplastic_human": int(_stratum(slug)["padj"].notna().sum())),
    Claim("strat_invitro_studies",
          r"\\emph\{In vitro\} \(human\) & (\d+) & [\d,]+ & [\d,]+ &",
          lambda slug="invitro_human": _stratum_studies(slug)),
    Claim("strat_invitro_features",
          r"\\emph\{In vitro\} \(human\) & \d+ & ([\d,]+) & [\d,]+ &",
          lambda slug="invitro_human": len(_stratum(slug))),
    Claim("strat_invitro_pooled",
          r"\\emph\{In vitro\} \(human\) & \d+ & [\d,]+ & ([\d,]+) &",
          lambda slug="invitro_human": int(_stratum(slug)["padj"].notna().sum())),
    # Sig column: the one column of Table 2 still unguarded after 2026-08-29.
    # SNI's Sig and CIPN's top gene were both stale, and nothing read them.
    Claim("strat_cci_sig",
          r"CCI & \d+ & [\d,]+ & [\d,]+ & ([\d,]+) &",
          lambda slug="CCI": _n_sig(_stratum(slug))),
    Claim("strat_sni_sig",
          r"SNI & \d+ & [\d,]+ & [\d,]+ & ([\d,]+) &",
          lambda slug="SNI": _n_sig(_stratum(slug))),
    Claim("strat_snl_sig",
          r"SNL & \d+ & [\d,]+ & [\d,]+ & ([\d,]+) &",
          lambda slug="SNL": _n_sig(_stratum(slug))),
    Claim("strat_cfa_sig",
          r"CFA & \d+ & [\d,]+ & [\d,]+ & ([\d,]+) &",
          lambda slug="CFA": _n_sig(_stratum(slug))),
    Claim("strat_cipn_sig",
          r"CIPN & \d+ & [\d,]+ & [\d,]+ & ([\d,]+) &",
          lambda slug="CIPN": _n_sig(_stratum(slug))),
    Claim("strat_neuro_human_sig",
          r"Neuropathic \(human\) & \d+ & [\d,]+ & [\d,]+ & ([\d,]+) &",
          lambda slug="neuropathic_human": _n_sig(_stratum(slug))),
    Claim("strat_lbp_sig",
          r"LBP \(human\) & \d+ & [\d,]+ & [\d,]+ & ([\d,]+) &",
          lambda slug="lbp_human": _n_sig(_stratum(slug))),
    Claim("strat_nocip_sig",
          r"Nociplastic \(human\) & \d+ & [\d,]+ & [\d,]+ & ([\d,]+) &",
          lambda slug="nociplastic_human": _n_sig(_stratum(slug))),
    Claim("strat_invitro_sig",
          r"\\emph\{In vitro\} \(human\) & \d+ & [\d,]+ & [\d,]+ & ([\d,]+) &",
          lambda slug="invitro_human": _n_sig(_stratum(slug))),
    # The eight strata the 2026-08-28 search amendment added. 06b has written
    # them since 2026-08-29; Table 2 showed only the original nine until
    # 2026-08-30, which is how a poolable human stratum (endometriosis: 14
    # studies, 126 significant) stayed out of a paper asserting that no human
    # subtype supports meta-analytic inference.
    Claim("strat_endo_human_studies",
          r"Endometriosis \(human\) & (\d+) & [\d,]+ & [\d,]+ &",
          lambda slug="endometriosis_human": _stratum_studies(slug)),
    Claim("strat_endo_human_features",
          r"Endometriosis \(human\) & \d+ & ([\d,]+) & [\d,]+ &",
          lambda slug="endometriosis_human": len(_stratum(slug))),
    Claim("strat_endo_human_pooled",
          r"Endometriosis \(human\) & \d+ & [\d,]+ & ([\d,]+) &",
          lambda slug="endometriosis_human": int(_stratum(slug)["padj"].notna().sum())),
    Claim("strat_endo_human_sig",
          r"Endometriosis \(human\) & \d+ & [\d,]+ & [\d,]+ & ([\d,]+) &",
          lambda slug="endometriosis_human": _n_sig(_stratum(slug))),
    Claim("strat_endo_rodent_studies",
          r"Endometriosis \(rodent\) & (\d+) & [\d,]+ & [\d,]+ &",
          lambda slug="endometriosis_rodent": _stratum_studies(slug)),
    Claim("strat_endo_rodent_features",
          r"Endometriosis \(rodent\) & \d+ & ([\d,]+) & [\d,]+ &",
          lambda slug="endometriosis_rodent": len(_stratum(slug))),
    Claim("strat_endo_rodent_pooled",
          r"Endometriosis \(rodent\) & \d+ & [\d,]+ & ([\d,]+) &",
          lambda slug="endometriosis_rodent": int(_stratum(slug)["padj"].notna().sum())),
    Claim("strat_endo_rodent_sig",
          r"Endometriosis \(rodent\) & \d+ & [\d,]+ & [\d,]+ & ([\d,]+) &",
          lambda slug="endometriosis_rodent": _n_sig(_stratum(slug))),
    Claim("strat_fibro_studies",
          r"Fibromyalgia \(human\) & (\d+) & [\d,]+ & [\d,]+ &",
          lambda slug="fibromyalgia_human": _stratum_studies(slug)),
    Claim("strat_fibro_features",
          r"Fibromyalgia \(human\) & \d+ & ([\d,]+) & [\d,]+ &",
          lambda slug="fibromyalgia_human": len(_stratum(slug))),
    Claim("strat_fibro_pooled",
          r"Fibromyalgia \(human\) & \d+ & [\d,]+ & ([\d,]+) &",
          lambda slug="fibromyalgia_human": int(_stratum(slug)["padj"].notna().sum())),
    Claim("strat_fibro_sig",
          r"Fibromyalgia \(human\) & \d+ & [\d,]+ & [\d,]+ & ([\d,]+) &",
          lambda slug="fibromyalgia_human": _n_sig(_stratum(slug))),
    Claim("strat_ibs_studies",
          r"IBS \(human\) & (\d+) & [\d,]+ & [\d,]+ &",
          lambda slug="ibs_human": _stratum_studies(slug)),
    Claim("strat_ibs_features",
          r"IBS \(human\) & \d+ & ([\d,]+) & [\d,]+ &",
          lambda slug="ibs_human": len(_stratum(slug))),
    Claim("strat_ibs_pooled",
          r"IBS \(human\) & \d+ & [\d,]+ & ([\d,]+) &",
          lambda slug="ibs_human": int(_stratum(slug)["padj"].notna().sum())),
    Claim("strat_ibs_sig",
          r"IBS \(human\) & \d+ & [\d,]+ & [\d,]+ & ([\d,]+) &",
          lambda slug="ibs_human": _n_sig(_stratum(slug))),
    Claim("strat_neuropathy_human_studies",
          r"Neuropathy \(human\) & (\d+) & [\d,]+ & [\d,]+ &",
          lambda slug="neuropathy_human": _stratum_studies(slug)),
    Claim("strat_neuropathy_human_features",
          r"Neuropathy \(human\) & \d+ & ([\d,]+) & [\d,]+ &",
          lambda slug="neuropathy_human": len(_stratum(slug))),
    Claim("strat_neuropathy_human_pooled",
          r"Neuropathy \(human\) & \d+ & [\d,]+ & ([\d,]+) &",
          lambda slug="neuropathy_human": int(_stratum(slug)["padj"].notna().sum())),
    Claim("strat_neuropathy_human_sig",
          r"Neuropathy \(human\) & \d+ & [\d,]+ & [\d,]+ & ([\d,]+) &",
          lambda slug="neuropathy_human": _n_sig(_stratum(slug))),
    Claim("strat_neuropathy_rodent_studies",
          r"Neuropathy \(rodent\) & (\d+) & [\d,]+ & [\d,]+ &",
          lambda slug="neuropathy_rodent": _stratum_studies(slug)),
    Claim("strat_neuropathy_rodent_features",
          r"Neuropathy \(rodent\) & \d+ & ([\d,]+) & [\d,]+ &",
          lambda slug="neuropathy_rodent": len(_stratum(slug))),
    Claim("strat_neuropathy_rodent_pooled",
          r"Neuropathy \(rodent\) & \d+ & [\d,]+ & ([\d,]+) &",
          lambda slug="neuropathy_rodent": int(_stratum(slug)["padj"].notna().sum())),
    Claim("strat_neuropathy_rodent_sig",
          r"Neuropathy \(rodent\) & \d+ & [\d,]+ & [\d,]+ & ([\d,]+) &",
          lambda slug="neuropathy_rodent": _n_sig(_stratum(slug))),
    Claim("strat_other_human_studies",
          r"Other \(human\) & (\d+) & [\d,]+ & [\d,]+ &",
          lambda slug="other_human_amendment": _stratum_studies(slug)),
    Claim("strat_other_human_features",
          r"Other \(human\) & \d+ & ([\d,]+) & [\d,]+ &",
          lambda slug="other_human_amendment": len(_stratum(slug))),
    Claim("strat_other_human_pooled",
          r"Other \(human\) & \d+ & [\d,]+ & ([\d,]+) &",
          lambda slug="other_human_amendment": int(_stratum(slug)["padj"].notna().sum())),
    Claim("strat_other_human_sig",
          r"Other \(human\) & \d+ & [\d,]+ & [\d,]+ & ([\d,]+) &",
          lambda slug="other_human_amendment": _n_sig(_stratum(slug))),
    Claim("strat_other_rodent_studies",
          r"Other \(rodent\) & (\d+) & [\d,]+ & [\d,]+ &",
          lambda slug="other_rodent_amendment": _stratum_studies(slug)),
    Claim("strat_other_rodent_features",
          r"Other \(rodent\) & \d+ & ([\d,]+) & [\d,]+ &",
          lambda slug="other_rodent_amendment": len(_stratum(slug))),
    Claim("strat_other_rodent_pooled",
          r"Other \(rodent\) & \d+ & [\d,]+ & ([\d,]+) &",
          lambda slug="other_rodent_amendment": int(_stratum(slug)["padj"].notna().sum())),
    Claim("strat_other_rodent_sig",
          r"Other \(rodent\) & \d+ & [\d,]+ & [\d,]+ & ([\d,]+) &",
          lambda slug="other_rodent_amendment": _n_sig(_stratum(slug))),
    # -- cross-species sensitivity (Table 4) --------------------------------
    *_sens_claims(),
    # -- corpus -------------------------------------------------------------
    Claim("features_total",
          r"across ([\d,]+) unique features passed quality control",
          lambda: len(_pooled())),
    Claim("studies_total_abstract",
          r"We meta-analysed (\d+) transcriptomic study units",
          lambda: _pooled_study_count()),
    Claim("features_total_abstract",
          r"transcriptomic study units \(([\d,]+) features",
          lambda: len(_pooled())),
    Claim("named_gene_pct",
          r"([\d.]+)\\% of them named genes",
          lambda: _named_gene_fraction(), tol=0.05),
    Claim("bh_family",
          r"This inferential family comprises the ([\d,]+)\s*\n?features measured",
          lambda: int(_pooled()["padj_pooled"].notna().sum())),
    Claim("significant_all",
          r"identified ([\d,]+) features significant at",
          lambda: _n_sig(_pooled(), "padj_pooled")),
    Claim("significant_median_k",
          r"they are pooled from a median of (\d+)\s*\ncontributing studies",
          lambda: _median_k("significant")),
    Claim("family_median_k",
          r"contributing studies, against (\d+) for the correction family",
          lambda: _median_k("family")),
    Claim("below_threshold",
          r"studies; the remaining ([\d,]+) of",
          lambda: len(_pooled()) - int(_pooled()["padj_pooled"].notna().sum())),
    # The denominator of that same sentence. Nothing read it, and it drifted to
    # 159,640 -- a figure inconsistent with its own sentence, since the two
    # parts of the split must add up to the corpus total.
    Claim("below_threshold_total",
          r"the remaining [\d,]+ of ([\d,]+) features fell below",
          lambda: len(_pooled())),
    # The study count in the corpus sentence. It said 82 while the Abstract,
    # the section heading, the Introduction, the Discussion and the PRISMA
    # caption all said 95, because this was the one of the six no claim read.
    Claim("corpus_study_units",
          r"Of these, (\d+) study units encompassing",
          lambda: _pooled_study_count()),

    # -- the pseudobulk single-cell arm -------------------------------------
    # Its own stratum, and until now it stated no number the registry read --
    # a figure no claim reads is outside every check this project has.
    Claim("sc_units",
          r"direction. (Ten) single-cell and single-nucleus series",
          lambda: 10.0 if _sc_studies() == 10 else -1.0,
          note="the word 'Ten' is the claim; a different count must rewrite it"),
    Claim("sc_samples",
          r"giving ([\d,]+) samples over",
          lambda: len(_sc_manifest())),
    Claim("sc_subjects",
          r"samples over ([\d,]+) subjects",
          lambda: int(_sc_manifest()["subject"].nunique())),
    Claim("sc_cells",
          r"subjects and\s+([\d,]+) cells,",
          lambda: int(_sc_manifest()["cells_kept"].sum())),
    Claim("sc_features",
          r"it estimates ([\d,]+) features",
          lambda: len(_sc_pooled())),
    Claim("sc_family",
          r"of which ([\d,]+) are measured in \$\\geq\$3 units",
          lambda: int(_sc_pooled()["padj_pooled"].notna().sum())),
    Claim("sc_median_k",
          r"at a median of (\w+) contributing units",
          lambda: int(_sc_pooled()[_sc_pooled()["padj_pooled"].notna()]["k"].median())),
    Claim("sc_ppm1l_padj",
          r"\(\\padj\\,=\\,([\d.]+), pooled from",
          lambda: float(_sc_pooled()
                        .set_index("feature_id").loc["Ppm1l", "padj_pooled"]),
          tol=0.0005),
    Claim("sc_ppm1l_k",
          r"pooled from (\w+) units\)",
          lambda: int(_sc_pooled().set_index("feature_id").loc["Ppm1l", "k"])),

    # -- between-model Q test ----------------------------------------------
    Claim("q_top_high_i2",
          r"(Four) of the five show near-complete heterogeneity",
          lambda: float(int((_top_q_i2()["I2_total"] > 99).sum()))),
    Claim("q_top_i2_exception",
          r"\\gene\{ARHGAP36\} does not, at I\\textsuperscript\{2\}\\,=\\,(\d+)",
          lambda: round(float(_top_q_i2()
                              .set_index("feature_id").loc["ARHGAP36", "I2_total"]))),
    Claim("q_testable",
          r"The between-model Q test identified [\d,]+ of ([\d,]+)",
          lambda: len(_read(STRAT / "between_model_Q.csv", usecols=["padj_QM"]))),
    Claim("q_testable_abstract",
          r"Of the ([\d,]+) features testable in \$\\geq\$2",
          lambda: len(_read(STRAT / "between_model_Q.csv", usecols=["padj_QM"]))),
    Claim("q_significant_abstract",
          r"models, ([\d,]+) showed significant between-model",
          lambda: int((_read(STRAT / "between_model_Q.csv",
                             usecols=["padj_QM"])["padj_QM"] < 0.05).sum())),
    Claim("q_significant",
          r"The between-model Q test identified ([\d,]+) of",
          lambda: int((_read(STRAT / "between_model_Q.csv",
                             usecols=["padj_QM"])["padj_QM"] < 0.05).sum())),
    # Table 2's Q row and the Discussion state the same two figures. Both were
    # stale and stale differently: the row said 33,853 tested and 3,685
    # significant, the sentence 62,792 and 13,192, and the Discussion restated
    # the row's 3,685. Nothing read any of the three.
    Claim("q_table_testable",
          r"Between-model Q & --- & ([\d,]+) &",
          lambda: len(_between_model_q())),
    Claim("q_table_significant",
          r"Between-model Q & --- & [\d,]+ & --- & \\textbf\{([\d,]+)\}",
          lambda: int((_between_model_q()["padj_QM"] < 0.05).sum())),
    Claim("q_significant_discussion",
          r"([\d,]+) features differ significantly between pain models",
          lambda: int((_between_model_q()["padj_QM"] < 0.05).sum())),
    # The permutation null. Without it "14,267 features differ by model" is a
    # statement about study-to-study variation as much as about models: 7.4% of
    # features look model-dependent under randomly assigned labels.
    Claim("q_null_pct",
          r"differences, against ([\d.]+)\\%",
          lambda: float(_q_permutation().query("perm > 0")["pct"].mean()),
          tol=0.05),
    Claim("q_null_p",
          r"\(\\pval\\,=\\,([\d.]+), 100 permutations\)",
          lambda: (int((_q_permutation().query("perm > 0")["pct"]
                        >= _q_permutation().query("perm == 0")["pct"].iloc[0]).sum()) + 1)
                  / (len(_q_permutation().query("perm > 0")) + 1),
          tol=0.005),
    # The Results state the empirical p as well as the Abstract, and only the
    # Abstract's copy was read: q_null_p needs the ", 100 permutations)" tail.
    Claim("q_null_p_results",
          r"exceeded by \w+ of 100 permutations \(\\pval\\,=\\,([\d.]+)\)",
          lambda: (int((_q_permutation().query("perm > 0")["pct"]
                        >= _q_permutation().query("perm == 0")["pct"].iloc[0]).sum()) + 1)
                  / (len(_q_permutation().query("perm > 0")) + 1),
          tol=0.005),
    Claim("q_null_pct_results",
          r"still makes ([\d.]+)\\% of features \(95th percentile",
          lambda: float(_q_permutation().query("perm > 0")["pct"].mean()),
          tol=0.05),
    Claim("q_null_p95",
          r"\(95th percentile ([\d.]+)\\%\)",
          lambda: float(_q_permutation().query("perm > 0")["pct"].quantile(0.95)),
          tol=0.05),
    Claim("q_null_attributable",
          r"some ([\d,]+) of the significant features are the\s*\nbetween-study variation",
          lambda: _q_null_attributable(100)),
    Claim("q_model_excess",
          r"model-attributable excess is closer to ([\d,]+)",
          lambda: _q_excess_over_null(100)),
    # The Discussion's copy of the excess. Unclaimed, and it said 8,600 while
    # the Results and the Highlights both said 8,800 from the same quantity.
    Claim("q_model_excess_discussion",
          r"some ([\d,]+) of them beyond a permuted-label null",
          lambda: _q_excess_over_null(100)),
    Claim("q_observed_pct",
          r"the observed rate is ([\d.]+)\\%\s*\nagainst that null",
          lambda: float(_q_permutation().query("perm == 0")["pct"].iloc[0]),
          tol=0.05),
    Claim("q_perms_exceeding",
          r"exceeded by (\w+) of 100 permutations",
          lambda: int((_q_permutation().query("perm > 0")["pct"]
                       >= _q_permutation().query("perm == 0")["pct"].iloc[0]).sum())),
    # The within-species positive control, which is what withdrew the species
    # reading of the cross-species null.
    Claim("hvh_random_pct",
          r"two halves agree on effect direction ([\d.]+)\\% of the time",
          lambda: float(_hvh("random")["pct_concordant"].mean()), tol=0.05),
    Claim("hvh_replicate_pct",
          r"endometriosis stratum agree ([\d.]+)\\% of the time",
          lambda: float(_hvh("replicate")["pct_concordant"].mean()), tol=0.05),
    Claim("hvh_replicate_pct_abstract",
          r"halves of a single human\s+condition agree on effect direction ([\d.]+)\\% of the time",
          lambda: float(_hvh("replicate")["pct_concordant"].mean()), tol=0.05),
    Claim("hvh_precise_pct",
          r"agree at ([\d.]+)\\%, no better than the set as a whole",
          lambda: float(_hvh("replicate")["pct_concordant_precise10"].mean()),
          tol=0.05),
    # Two statements the rewritten passage now makes that the old one did not,
    # both of which reversed: the control used to find no feature significant
    # in both halves of any split, and concordance used to fall with depth.
    Claim("hvh_both_significant",
          r"and (\d+) features across the nine\s*\nsplits reach",
          lambda: int(sum(int(_hvh(d)["n_both_significant"].sum())
                          for d in ("random", "compartment", "replicate")))),
    Claim("hvh_depth_pct",
          r"reaching ([\d.]+)\\% on features measured in",
          lambda: float(_hvh("replicate")["pct_concordant_k3"].mean()), tol=0.05),

    # -- pain-model strata --------------------------------------------------
    Claim("snl_significant",
          r"with ([\d,]+) features significant among the [\d,]+ pooled",
          lambda: _n_sig(_stratum("SNL"))),
    Claim("snl_family",
          r"features significant among the ([\d,]+) pooled in \$\\geq\$3 studies",
          lambda: int(_stratum("SNL")["padj"].notna().sum())),
    Claim("cci_significant",
          r"The CCI model[\s\S]{0,60}?yielded ([\d,]+) cross-study\s+significant features",
          lambda: _n_sig(_stratum("CCI"))),
    Claim("cipn_significant",
          r"The CIPN stratum \(\w+ studies\) yielded ([\d,]+) cross-study\s+significant features",
          lambda: _n_sig(_stratum("CIPN"))),
    Claim("sni_significant",
          r"yielded only ([\d,]+) cross-study significant hits",
          lambda: _n_sig(_stratum("SNI"))),
    Claim("sni_bonferroni",
          r"of which just (\w+)\s+survive Bonferroni correction",
          lambda: _bonferroni_survivors("SNI")),
    Claim("sni_studies_prose",
          r"most-studied model\s*\nat (\d+) study units",
          lambda: _stratum_studies("SNI")),
    # 06d's tissue meta-regression. Never claimed, and it had gone from 38
    # genes to 16 -- with the microglial markers the sentence was built on
    # falling to padj 0.95 -- because 06d reads 06b's SNI pool and the tissue
    # map, and both moved when GSE197233 was split per region.
    Claim("sni_tissue_genes",
          r"identified (\d+) genes with\s*\nsignificant tissue-dependent expression",
          lambda: int((_sni_tissue()["padj"] < 0.05).sum())),
    # The denominator, which the sentence did not state. 16 of 116 candidates
    # is a different claim from 16 of the transcriptome, and the missing
    # denominator is how the microglial reading drifted unchallenged.
    Claim("sni_tissue_tested",
          r"expression among the ([\d,]+) it could test",
          lambda: len(_sni_tissue())),
    Claim("sni_tissue_tested_recap",
          r"\\gene\{C1qa\} and \\gene\{C1qb\} are among\s+those ([\d,]+) and are null",
          lambda: len(_sni_tissue())),
    # 06d's candidate cap, read out of 06d rather than transcribed: Methods
    # describes the pre-filter, so a change to the cap silently makes Methods
    # wrong. Same reason build_table_s1.py parses STUDY_MODEL out of 06b.
    Claim("sni_tissue_top_n",
          r"fall outside the largest ([\d,]+) pooled effects\s+that the candidate set admits",
          lambda: _sni_tissue_top_n()),
    # How many strata yield nothing, and the feature range over the five the
    # sentence is about. The range read 21,914--47,060, an upper bound no
    # stratum has carried for two regenerations.
    Claim("strata_total",
          r"of the (\w+) strata in Table~\\ref\{tab:stratified\} yielded no",
          lambda: len(_all_strata())),
    Claim("strata_no_result",
          r"(\w+) of the \w+ strata in Table~\\ref\{tab:stratified\} yielded no",
          lambda: sum(1 for n in _all_strata()
                      if int(_stratum(n)["padj"].notna().sum()) == 0)),
    Claim("zero_strata_features_min",
          r"so across ([\d,]+)--[\d,]+ features\s+per stratum",
          lambda: min(len(_stratum(n)) for n in _ORIGINAL_EMPTY_STRATA)),
    Claim("zero_strata_features_max",
          r"so across [\d,]+--([\d,]+) features\s+per stratum",
          lambda: max(len(_stratum(n)) for n in _ORIGINAL_EMPTY_STRATA)),
    # The endometriosis stratum is the one human subtype that reaches poolable
    # depth. Until 2026-08-29 the paper said no human subtype does.
    Claim("endo_studies",
          r"an endometriosis stratum of\s+(\d+) studies",
          lambda: _stratum_studies("endometriosis_human")),
    Claim("endo_pooled",
          r"studies that pools ([\d,]+) features",
          lambda: int(_stratum("endometriosis_human")["padj"].notna().sum())),
    Claim("endo_significant",
          r"features and yields (\d+) significant",
          lambda: _n_sig(_stratum("endometriosis_human"))),
    # How many of the stratum's studies merge levels of a separating field.
    # Stated because two of them merge in a condition that is not
    # endometriosis (see conf/analysis/phenotype_scope.csv), so the count is a
    # composition caveat rather than a diagnostic.
    Claim("endo_merged_arms",
          r"substituting another: (\w+) of the\s+\d+ merge levels of a separating field",
          lambda: _stratum_merged_arms("endometriosis_human")),
    Claim("endo_studies_recap",
          r"substituting another: \w+ of the\s+(\d+) merge levels of a separating field",
          lambda: _stratum_studies("endometriosis_human")),

    # -- cross-species ------------------------------------------------------
    Claim("pairs_mouse",
          r"orthologous gene pairs from mouse\s*\n?\(([\d,]+)\)",
          lambda: int(_summary_row("Mus musculus")["n_ortholog_pairs"])),
    Claim("pairs_rat",
          r"and rat \(([\d,]+)\) to human\.",
          lambda: int(_summary_row("Rattus norvegicus")["n_ortholog_pairs"])),
    Claim("pairs_total",
          r"species separately and mapped ([\d,]+) orthologous gene pairs",
          lambda: int(_cs_summary()["n_ortholog_pairs"].sum())),
    Claim("pairs_total_discussion",
          r"the scarcity itself: across ([\d,]+) ortholog pairs",
          lambda: int(_cs_summary()["n_ortholog_pairs"].sum())),
    Claim("concordance_mouse",
          r"same sign of meta-analytic effect size,\s+was ([\d.]+)\\% for mouse",
          lambda: _summary_row("Mus musculus")["pct_concordant"], tol=0.05),
    Claim("concordance_rat",
          r"was [\d.]+\\% for mouse and ([\d.]+)\\% for\s*\n?rat",
          lambda: _summary_row("Rattus norvegicus")["pct_concordant"], tol=0.05),
    # The table and the Discussion carried pre-regeneration values -- 48.0/49.6%
    # with rho -0.068/-0.034 -- while the abstract and Results, which were
    # claim-tested, carried the correct 48.7/48.5%. The stale pair were also
    # swapped between species. Nothing covered the table, the rho values or the
    # Discussion range, so they drifted silently through a regeneration.
    # The rest of Table 3. Only its Concordant column was claimed, so its
    # pairs, percentage and rho columns restated four prose figures with
    # nothing checking that they still agreed with them.
    Claim("rho_mouse",
          r"\(mouse:\s*\n\$\\rho\$\\,=\\,([\d.]+),",
          lambda: abs(float(_summary_row("Mus musculus")["spearman_r"])),
          tol=0.0005),
    Claim("rho_rat",
          r"; rat:\s*\n\$\\rho\$\\,=\\,([\d.]+),",
          lambda: abs(float(_summary_row("Rattus norvegicus")["spearman_r"])),
          tol=0.0005),
    Claim("coverage_mouse",
          r"ortholog: ([\d,]+) of [\d,]+ for mouse",
          lambda: int(_summary_row("Mus musculus")["n_mapped_to_human"])),
    Claim("coverage_mouse_denom",
          r"ortholog: [\d,]+ of ([\d,]+) for mouse",
          lambda: int(_summary_row("Mus musculus")["n_animal_features"])),
    Claim("coverage_rat",
          r"and ([\d,]+) of [\d,]+ for\s+rat",
          lambda: int(_summary_row("Rattus norvegicus")["n_mapped_to_human"])),
    # The rest of the concordance sentence. Its percentages and Spearman rho
    # carried claims; its binomial p-values, Spearman p-values and Pearson
    # coefficients did not, and on 2026-08-31 all six were stale -- rat's
    # Pearson r had gone from -0.002 to +0.006, so the sentence still called it
    # "negative and not significant" when it was neither.
    Claim("binom_p_mouse",
          r"mouse falls above it, binomial\s*\n\\pval\\,=([\d.]+)"
          r"\$\\times\$10\\textsuperscript\{-3\};",
          lambda: _binomial_p("Mus musculus") * 1e3, tol=0.05),
    Claim("binom_p_rat",
          r"rat sits further above,\s*\n\\pval\\,=([\d.]+)"
          r"\$\\times\$10\\textsuperscript\{-8\}\)",
          lambda: _binomial_p("Rattus norvegicus") * 1e8, tol=0.05),
    Claim("rho_p_mouse",
          r"\$\\rho\$\\,=\\,[\d.]+, \\pval\\,=([\d.]+)"
          r"\$\\times\$10\\textsuperscript\{-7\}; rat:",
          lambda: float(_summary_row("Mus musculus")["spearman_p"]) * 1e7, tol=0.05),
    Claim("rho_p_rat",
          r"\$\\rho\$\\,=\\,[\d.]+, \\pval\\,=([\d.]+)\$\\times\$10\\textsuperscript\{-20\}\)",
          lambda: float(_summary_row("Rattus norvegicus")["spearman_p"]) * 1e20,
          tol=0.05),
    Claim("pearson_r_mouse",
          r"\(mouse:\s*\n\\textit\{r\}\\,=\\,([\d.]+),",
          lambda: float(_summary_row("Mus musculus")["pearson_r"]), tol=0.0005),
    Claim("pearson_p_mouse",
          r"\\textit\{r\}\\,=\\,[\d.]+, \\pval\\,=([\d.]+)"
          r"\$\\times\$10\\textsuperscript\{-10\}; rat:",
          lambda: float(_summary_row("Mus musculus")["pearson_p"]) * 1e10, tol=0.05),
    Claim("pearson_r_rat",
          r"; rat:\s*\n\\textit\{r\}\\,=\\,((?:\$-\$)?[\d.]+), \\pval\\,=[\d.]+\$",
          lambda: float(_summary_row("Rattus norvegicus")["pearson_r"]),
          tol=0.0005),
    Claim("pearson_p_rat",
          r"; rat:\s*\n\\textit\{r\}\\,=\\,(?:\$-\$)?[\d.]+, \\pval\\,=([\d.]+)"
          r"\$\\times\$10\\textsuperscript\{-11\}\)",
          lambda: float(_summary_row("Rattus norvegicus")["pearson_p"]) * 1e11,
          tol=0.05),
    # The abstract restated the SNL count and both concordance percentages with
    # no claim behind them, so they survived a regeneration that moved every
    # one of them. Registered now.
    Claim("abstract_snl_significant",
          r"model yielded the most\s*\n\s*significant features \((\d[\d,]*)\)",
          lambda: _n_sig(_stratum("SNL"))),
    # The tissue-matched arm as the prose states it. Table 4's row carries the
    # same figures and is claimed, but the sentence restating them was not --
    # and it said "neither significant" of a mouse coefficient that is now
    # significant, and quoted the published-pool rho in the Discussion while
    # describing the tissue-matched pool.
    Claim("tm_pct_mouse",
          r"takes concordance to ([\d.]+)\\% for mouse",
          lambda: float(_sens_row("tissue_matched", "Mus musculus")["pct_concordant"]),
          tol=0.05),
    Claim("tm_pct_rat",
          r"and ([\d.]+)\\%\s*\nfor rat \(\\pval",
          lambda: float(_sens_row("tissue_matched", "Rattus norvegicus")["pct_concordant"]),
          tol=0.05),
    Claim("tm_binom_mouse",
          r"\(binomial\s*\n\\pval\\,=([\d.]+)\$\\times\$10\\textsuperscript\{-3\}\)",
          lambda: _sens_binomial_p("tissue_matched", "Mus musculus") * 1e3,
          tol=0.05),
    Claim("tm_binom_rat",
          r"for rat \(\\pval\\,=\\,([\d.]+)\)",
          lambda: _sens_binomial_p("tissue_matched", "Rattus norvegicus"),
          tol=0.005),
    Claim("tm_rho_mouse",
          r"\(\$-\$([\d.]+) for mouse and \$-\$[\d.]+ for rat;",
          lambda: abs(float(_sens_row("tissue_matched",
                                      "Mus musculus")["spearman_r"])),
          tol=0.0005),
    Claim("tm_rho_rat",
          r"\(\$-\$[\d.]+ for mouse and \$-\$([\d.]+) for rat;",
          lambda: abs(float(_sens_row("tissue_matched",
                                      "Rattus norvegicus")["spearman_r"])),
          tol=0.0005),
    Claim("tm_pairs_mouse_prose",
          r"which\s+is a statement about the ([\d,]+) pairs behind the test",
          lambda: int(_sens_row("tissue_matched",
                                "Mus musculus")["n_ortholog_pairs"])),
    # The abstract and the Discussion restate the arm's two percentages.
    Claim("abstract_tm_rat",
          r"for mouse \(([\d.]+)\\% and",
          lambda: float(_sens_row("tissue_matched",
                                  "Rattus norvegicus")["pct_concordant"]),
          tol=0.05),
    Claim("abstract_tm_mouse",
          r"for mouse \([\d.]+\\% and ([\d.]+)\\%\)",
          lambda: float(_sens_row("tissue_matched", "Mus musculus")["pct_concordant"]),
          tol=0.05),
    Claim("discussion_tm_low",
          r"concordance of ([\d.]+)--[\d.]+\\% in tissue-comparable",
          lambda: min(float(_sens_row("tissue_matched", sp)["pct_concordant"])
                      for sp in ("Mus musculus", "Rattus norvegicus")),
          tol=0.05),
    Claim("discussion_tm_high",
          r"concordance of [\d.]+--([\d.]+)\\% in tissue-comparable",
          lambda: max(float(_sens_row("tissue_matched", sp)["pct_concordant"])
                      for sp in ("Mus musculus", "Rattus norvegicus")),
          tol=0.05),
    Claim("discussion_tm_rho_mouse",
          r"\(\$\\rho\$ of \$-\$([\d.]+) for mouse",
          lambda: abs(float(_sens_row("tissue_matched",
                                      "Mus musculus")["spearman_r"])),
          tol=0.0005),
    Claim("discussion_tm_rho_rat",
          r"for mouse and \$-\$([\d.]+) for rat\) suggest",
          lambda: abs(float(_sens_row("tissue_matched",
                                      "Rattus norvegicus")["spearman_r"])),
          tol=0.0005),
    Claim("abstract_concordance_rat",
          r"concordance was\s*\n\s*([\d.]+)\\% for rat",
          lambda: float(_summary_row("Rattus norvegicus")["pct_concordant"]),
          tol=0.05),
    Claim("abstract_concordance_mouse",
          r"for rat and ([\d.]+)\\% for mouse, falling",
          lambda: float(_summary_row("Mus musculus")["pct_concordant"]),
          tol=0.05),
    # The SNL power caveat. Four numbers that would otherwise drift silently,
    # which is how the abstract's figures survived a regeneration untouched.
    Claim("snl_pct_significant",
          r"That (\d+)\\% of pooled features reach\s+significance",
          lambda: _snl_stats()["pct_sig"], tol=0.5),
    Claim("snl_median_g_sig",
          r"median \$\|g\|\$\\,=\\,([\d.]+) against",
          lambda: _snl_stats()["median_abs_g_sig"], tol=0.005),
    Claim("snl_median_g_all",
          r"against ([\d.]+) across the",
          lambda: _snl_stats()["median_abs_g_all"], tol=0.005),
    Claim("snl_median_i2",
          r"heterogeneous \(median I\\textsuperscript\{2\}\\,=\\,(\d+)\)",
          lambda: _snl_stats()["median_i2_sig"], tol=0.5),
    Claim("concordant_transcripts",
          r"There are (\d+) cross-species\s+concordant transcripts",
          lambda: len(_concordant_transcripts())),
    # The paragraph restates the pair a third time to compute that 16%. Both
    # halves are stratified outputs and move whenever 06b is re-run.
    # Both species, separately: a single claim on the union would pass while
    # one side silently gained a pair.
    Claim("mouse_both_significant",
          r"across ([\d,]+) mouse and [\d,]+ rat ortholog",
          lambda: len(_concordance("musculus"))),
    Claim("rat_both_significant_pairs",
          r"across [\d,]+ mouse and ([\d,]+) rat ortholog",
          lambda: len(_concordance("norvegicus"))),

    # -- human pool ---------------------------------------------------------
    # -- depth and precision (06h) ------------------------------------------
    Claim("median_k_mouse",
          r"a median of (\d+) contributing studies for mouse",
          lambda: _species_median_k("Mus_musculus"), source="supplementary"),
    Claim("median_k_rat",
          r"contributing studies for mouse against (\d+) for\s*\nrat",
          lambda: _species_median_k("Rattus_norvegicus"), source="supplementary"),
    Claim("subsampled_pct_mouse_features",
          r"touching ([\d.]+)\\% of mouse features",
          lambda: _depth("Mus musculus")["pct_subsampled"], tol=0.05, source="supplementary"),
    Claim("subsampled_pct_rat_features",
          r"of mouse features against ([\d.]+)\\% of rat's",
          lambda: _depth("Rattus norvegicus")["pct_subsampled"], tol=0.05, source="supplementary"),
    Claim("subsampled_mouse_mean",
          r"rises from [\d.]+\\% to ([\d.]+)\\% \(five",
          lambda: _depth("Mus musculus")["mean"], tol=0.005, source="supplementary"),
    Claim("subsampled_rat_mean",
          r"and rat's falls from [\d.]+\\% to ([\d.]+)\\%",
          lambda: _depth("Rattus norvegicus")["mean"], tol=0.005, source="supplementary"),
    Claim("shared_mouse_mean",
          r"leaves ([\d.]+)\\% \([\d.]+--[\d.]+\\%\) for mouse against",
          lambda: _depth("Mus musculus", shared=True)["mean"], tol=0.005, source="supplementary"),
    Claim("shared_rat_mean",
          r"for mouse against ([\d.]+)\\%\s*\n\([\d.]+--[\d.]+\\%\) for rat",
          lambda: _depth("Rattus norvegicus", shared=True)["mean"], tol=0.005,
          source="supplementary"),
    Claim("median_se_mouse",
          r"is\s*\n([\d.]+) for mouse against [\d.]+ for rat",
          lambda: _median_per_study_se("Mus musculus"), tol=0.0005, source="supplementary"),
    Claim("median_se_rat",
          r"[\d.]+ for mouse against ([\d.]+) for rat",
          lambda: _median_per_study_se("Rattus norvegicus"), tol=0.0005, source="supplementary"),
    Claim("mouse_significant",
          r"only (\w+) of [\d,]+ mouse correction-family",
          lambda: _n_sig(_species("Mus_musculus"))),
    # The replicate ranges, and the three per-species significant fractions the
    # paragraph compares. None was claimed: the argument that a roughly
    # fivefold gap "survives every control available to us" rests on the ranges
    # not overlapping, and on 2026-08-31 every one of the eight was stale.
    Claim("subsampled_mouse_lo",
          r"seeded replicates, range ([\d.]+)--[\d.]+\\%\)",
          lambda: _depth("Mus musculus")["lo"], tol=0.005, source="supplementary"),
    Claim("subsampled_mouse_hi",
          r"seeded replicates, range [\d.]+--([\d.]+)\\%\)",
          lambda: _depth("Mus musculus")["hi"], tol=0.005, source="supplementary"),
    Claim("subsampled_rat_lo",
          r"falls from [\d.]+\\% to [\d.]+\\%\s*\n\(([\d.]+)--[\d.]+\\%\)",
          lambda: _depth("Rattus norvegicus")["lo"], tol=0.005, source="supplementary"),
    Claim("subsampled_rat_hi",
          r"falls from [\d.]+\\% to [\d.]+\\%\s*\n\([\d.]+--([\d.]+)\\%\)",
          lambda: _depth("Rattus norvegicus")["hi"], tol=0.005, source="supplementary"),
    Claim("shared_mouse_lo",
          r"leaves [\d.]+\\% \(([\d.]+)--[\d.]+\\%\) for mouse",
          lambda: _depth("Mus musculus", shared=True)["lo"], tol=0.005, source="supplementary"),
    Claim("shared_mouse_hi",
          r"leaves [\d.]+\\% \([\d.]+--([\d.]+)\\%\) for mouse",
          lambda: _depth("Mus musculus", shared=True)["hi"], tol=0.005, source="supplementary"),
    Claim("shared_rat_lo",
          r"for mouse against [\d.]+\\%\s*\n\(([\d.]+)--[\d.]+\\%\) for rat",
          lambda: _depth("Rattus norvegicus", shared=True)["lo"], tol=0.005,
          source="supplementary"),
    Claim("shared_rat_hi",
          r"for mouse against [\d.]+\\%\s*\n\([\d.]+--([\d.]+)\\%\) for rat",
          lambda: _depth("Rattus norvegicus", shared=True)["hi"], tol=0.005,
          source="supplementary"),
    Claim("mouse_baseline_pct",
          r"the mouse significant fraction rises from ([\d.]+)\\% to",
          lambda: _species_sig_pct("Mus_musculus"), tol=0.005, source="supplementary"),
    Claim("rat_baseline_pct",
          r"and rat's falls from ([\d.]+)\\% to",
          lambda: _species_sig_pct("Rattus_norvegicus"), tol=0.005, source="supplementary"),
    Claim("human_species_pct",
          r"correction family \(([\d.]+)\\%\), a rate between",
          lambda: _species_sig_pct("Homo_sapiens"), tol=0.005),
    Claim("tau2_rise_k_rat",
          r"thereafter, from\s*\n\$k\$\\,=\\,(\d+) in rat but",
          lambda: _tau2_rise_k("Rattus_norvegicus"), source="supplementary"),
    Claim("tau2_rise_k_mouse",
          r"in rat but only from \$k\$\\,=\\,(\d+) in mouse",
          lambda: _tau2_rise_k("Mus_musculus"), source="supplementary"),

    # Two more copies of figures the section already states. Both were
    # unclaimed and both were stale: the Results said the human arm had 23
    # studies where the sentence above it said 24, and the tissue paragraph
    # said 13 endometrial studies against a stratum that had grown to 21.
    Claim("human_pool_studies_results",
          r"productive of the\s+three: its (\d+) human studies yield",
          lambda: _species_studies("Homo_sapiens")),
    Claim("endometrial_in_human_pool",
          r"The human pool contains\s*\n(\d+) endometrial studies",
          lambda: _endometrial_human_units()),
    Claim("human_significant",
          r"(\d+) significant features of [\d,]+ in their\s*\ncorrection family",
          lambda: _n_sig(_species("Homo_sapiens"))),
    Claim("human_family",
          r"significant features of ([\d,]+) in their",
          lambda: int(_species("Homo_sapiens")["padj"].notna().sum())),
    Claim("human_pool_studies",
          r"escapes these limits only by combining all (\d+) human\s+studies",
          lambda: _species_studies("Homo_sapiens")),
    Claim("human_pool_studies_discussion",
          r"exists only because all (\d+) human studies\s+are",
          lambda: _species_studies("Homo_sapiens")),
    # Methods restates the depth of all three species pools. It said 26 mouse,
    # 26 rat and 5 human studies against pools of 31, 28 and 34, the human
    # figure predating the search amendment, and nothing read the sentence.
    Claim("species_studies_methods_mouse",
          r"within each species \((\d+) mouse,",
          lambda: _species_studies("Mus_musculus")),
    Claim("species_studies_methods_rat",
          r"\(\d+ mouse, (\d+) rat and \d+ human",
          lambda: _species_studies("Rattus_norvegicus")),
    Claim("species_studies_methods_human",
          r"\(\d+ mouse, \d+ rat and (\d+) human\s+studies\)",
          lambda: _species_studies("Homo_sapiens")),
    Claim("human_family_limitations",
          r"which yields ([\d,]+) features measured in\s+\$\\geq\$3 studies",
          lambda: int(_species("Homo_sapiens")["padj"].notna().sum())),
    Claim("mouse_family_size",
          r"only \w+ of ([\d,]+) mouse correction-family",
          lambda: int(_species("Mus_musculus")["padj"].notna().sum())),

    # -- PRISMA disposition (Methods) ---------------------------------------
    # Composition of the search amendment, disclosed in the Methods.
    Claim("amendment_endometrial",
          r"and (\d+) of the \d+ admitted are endometrial",
          lambda: _stratum_studies("endometriosis_human"), source="supplementary"),
    Claim("amendment_admitted",
          r"and \d+ of the (\d+) admitted are endometrial",
          lambda: _amendment_admitted(), source="supplementary"),
    Claim("prisma_include",
          r"Of the [\d,]+ candidates, (\d+) are\s+included",
          lambda: int(_prisma_verdicts().get("include", 0))),
    Claim("prisma_exclude",
          r"are\s+included, (\d+) excluded",
          lambda: int(_prisma_verdicts().get("exclude", 0))),
    Claim("prisma_pending",
          r"excluded and (\d+) remain pending",
          lambda: int(_prisma_verdicts().get("pending", 0))),

    # -- PRISMA figure caption ----------------------------------------------
    # The caption exists to stop a reader reconciling the diagram's 47 with the
    # 79 datasets pooled: the two count different things and always will. Every
    # number in it is claimed, and the patterns are anchored to the caption's
    # own wording because several of these figures are stated elsewhere too --
    # a pattern matching two places fails by design.
    Claim("prisma_caption_records",
          r"query above and returned ([\d,]+) records",
          lambda: int(_search_log_last()["after_dedup"]), source="supplementary"),
    Claim("prisma_caption_retained",
          r"hand or through their GEO deposits, (\d+) in all",
          lambda: _literature_included(), source="supplementary"),
    # The section heading states the same count. It said 80 against a corpus
    # of 78 and then 82, having never been read by anything.
    Claim("results_heading_units",
          r"Transcriptomic meta-analysis of (\d+) study units",
          lambda: _pooled_study_count()),
    # Methods states how many accessions the independence registry holds. It
    # said seven while the registry held nine, the two multi-region splits
    # having been added without the sentence being touched.
    Claim("superseded_accessions",
          r"([A-Za-z]+) transcriptomic accessions were excluded to preserve",
          lambda: sum(1 for sid in _read(SUPERSEDED)["study_id"].astype(str)
                      if sid.startswith("GSE")), source="supplementary"),
    # The contrast audit's own summary, in Methods. All three move whenever
    # scripts/audit_contrasts.py is re-run over a changed corpus, which is how
    # the sheet came to describe pre-repair arms for four studies.
    Claim("audit_flagged",
          r"and (\d+) of\s*\nthe \d+ audited accessions carry at least one flag",
          lambda: int(_contrast_audit()["flags"].notna().sum()), source="supplementary"),
    Claim("audit_total",
          r"and \d+ of\s*\nthe (\d+) audited accessions carry at least one flag",
          lambda: len(_contrast_audit()), source="supplementary"),
    Claim("audit_title_only",
          r"being the only thing that separates the arms in (\d+) accessions",
          lambda: int((_contrast_audit()["separating_fields"].fillna("")
                       .str.strip() == "sample title").sum()), source="supplementary"),
    # The Discussion's opening restates the corpus size and the ortholog count.
    # Both were stale by a wide margin -- 58 studies against 82, and an
    # ortholog count from two regenerations earlier -- because the sentence
    # reads as prose and nothing checked it.
    Claim("discussion_tx_units",
          r"encompassing ([\d,]+) transcriptomic\s*\nstudy units",
          lambda: _pooled_study_count()),
    Claim("discussion_pairs",
          r"([\d,]+) cross-species\s*\northolog pairs",
          lambda: int(_cs_summary()["n_ortholog_pairs"].sum())),
    # The corpus size is stated in six places in the manuscript and the
    # ortholog count in three. Each occurrence gets its own claim rather than
    # being reworded away: they are all load-bearing sentences, and on
    # 2026-08-31 the Summary said 79, the Introduction 78 and the forest
    # caption 78 against a corpus of 82.
    #
    # summary_tx_units and summary_pairs are gone with the "Summary."
    # paragraph they read, removed on 2026-09-11. Both figures are still
    # claimed elsewhere -- the study count by intro_tx_units and four others,
    # the ortholog pairs by pairs_total -- so neither lost its only binding.
    # The Introduction's copy of the ortholog count. It is the fourth place the
    # figure appears and was the only one no claim read, so it alone stayed at
    # 33,033 while the other three tracked the analysis.
    Claim("forest_caption_units",
          r"measured in \$\\geq\$3 of the ([\d,]+) study units can be pooled",
          lambda: _pooled_study_count()),

    # -- Supplementary Table S1 ---------------------------------------------
    # S1 is regenerated from the tracked results by scripts/build_table_s1.py.
    # Its counts sat at the 58-study corpus until 2026-08-29 because the prose
    # describing it was read by nothing.
    Claim("s1_tissue_annotated",
          r"Tissue is recorded only for the (\d+) SNI studies",
          lambda: int((_table_s1()["Tissue"].astype(str) != "not annotated").sum()),
          source="supplementary"),
    Claim("s1_counts_recoverable",
          r"recoverable for (\d+) of the [\d,]+ studies",
          lambda: int(pd.to_numeric(_table_s1()["N_Analysed"],
                                    errors="coerce").notna().sum()), source="supplementary"),
    Claim("s1_studies",
          r"recoverable for \d+ of the ([\d,]+) studies",
          lambda: len(_table_s1()), source="supplementary"),
    Claim("s1_samples",
          r"of the [\d,]+ studies \(([\d,]+) samples in total\)",
          lambda: int(pd.to_numeric(_table_s1()["N_Analysed"],
                                    errors="coerce").sum()), source="supplementary"),

    # -- genomics -----------------------------------------------------------
    # -- Figure: GWAS replication against the transcriptomic pool -----------
    # The figure this replaced was drawn on the May 2026 pool and never
    # checked: its genes carried k=3 where the August pool gives k=15-23, and
    # every effect size in it was superseded, ASTN2 with the opposite sign.
    # Nothing caught that because the figure stated no number the registry
    # read. These two claims are what stop it recurring.
    Claim("gwas_fig_top_replication",
          r"reporting them; \\gene\{LRP1\} leads at (\d+)",
          lambda: int(_gwas_replication()["n_studies"].max())),
    # The depth range in the same sentence. Nothing read it, so it still said
    # "6 to 23" against a pool whose range had moved to 7 to 33 -- the figure
    # was regenerated, the caption's other number was, and this one was not.
    Claim("gwas_fig_min_k",
          r"pooled from (\d+) to \d+ studies",
          lambda: int(_gwas_fig_effects()["k"].min())),
    Claim("gwas_fig_max_k",
          r"pooled from \d+ to (\d+) studies",
          lambda: int(_gwas_fig_effects()["k"].max())),
    Claim("gwas_fig_min_padj",
          r"smallest adjusted \\pval-value among\s*\n\s*them being ([\d.]+)",
          lambda: float(_gwas_fig_effects()["padj_pooled"].min()), tol=0.005,
          note="the panel's point is that none of the genetically strongest "
               "genes comes near significance; a value that fell below 0.05 "
               "would change what the figure shows"),
    Claim("magma_genes",
          r"identified ([\d,]+) Bonferroni-significant genes",
          lambda: len(_read(GENOMICS / "magma" / "significant_genes.csv"))),
    Claim("magma_genes_convergence",
          r"concordant transcripts for the ([\d,]+)\s+MAGMA-significant genes",
          lambda: len(_read(GENOMICS / "magma" / "significant_genes.csv"))),

    # -- publication bias ---------------------------------------------------
    Claim("egger_features",
          r"Egger\'s regression test \\cite\{egger1997bias\}, applied to the ([\d,]+) features",
          lambda: len(_read(TX_META / "egger_test.csv", usecols=["k"]))),
    Claim("effect_estimates_total",
          r"study units encompassing ([\d,]+) effect-size estimates",
          lambda: _effect_estimate_count()),
    # The median depth of the Egger set, which nothing read: it stayed at 22
    # across regenerations that took it to 31.
    Claim("egger_median_k",
          r"studies \(median \$k\$\\,=\\,(\d+)\), is significant",
          lambda: int(_egger()["k"].median())),
    Claim("egger_flagged",
          r"is significant for ([\d,]+) of them",
          lambda: int(_egger()["bias_flagged"].sum())),
    # The share, which nothing read: it stayed at 39.7% while its own numerator
    # and denominator were both claimed and both moved.
    Claim("egger_flagged_pct",
          r"is significant for [\d,]+ of them \(([\d.]+)\\% at",
          lambda: 100.0 * float(_egger()["bias_flagged"].sum()) / len(_egger()),
          tol=0.05),
    # The direction split among flagged features. It read "12,405 positive
    # versus 10,927 negative", numbers matching neither the flagged set nor the
    # tested set of any recent run, and supported a symmetry argument the
    # current split (11,141 v 6,577) does not.
    Claim("egger_z_positive",
          r"leans positive \(([\d,]+) with \$z\$\\,\$>\$\\,0",
          lambda: int((_egger_z()["egger_z"] > 0).sum())),
    Claim("egger_z_negative",
          r"against ([\d,]+) with \$z\$\\,\$<\$\\,0\)",
          lambda: int((_egger_z()["egger_z"] < 0).sum())),

    # -- proteomics ---------------------------------------------------------
    Claim("proteomic_combinations",
          r"Of ([\d,]+) peptide--region--condition\s+combinations",
          lambda: len(_read(PROTEOMICS / "PXD013362_effects_normalized.csv",
                            usecols=["feature_id"]))),
    Claim("proteomic_significant",
          r"only (\w+) peptides reached significance",
          lambda: int((_read(PROTEOMICS / "PXD013362_effects_normalized.csv",
                             usecols=["padj"])["padj"] < 0.05).sum())),

    # -- proteomics: screening and triage -----------------------------------
    # These read the PRISMA tables, which are tracked, so they never skip.
    Claim("proteomic_datasets",
          r"A keyword search of PRIDE returned (\d+) pain-targeted",
          lambda: len(_triage())),
    Claim("proteomic_raw_only",
          r"these, (\d+) are deposited as instrument or peak files only",
          lambda: int((_triage()["verdict"] == "raw_only").sum())),
    Claim("proteomic_processed",
          r"and (\w+) carry a\s+processed quantification table",
          lambda: int((_triage()["verdict"] == "processed_table").sum())),
    Claim("proteomic_usable",
          r"of which (\w+) have a usable pain-versus-control",
          lambda: int((_read(PRISMA / "proteomics_manual_review.csv",
                             usecols=["verdict"])["verdict"] == "include").sum())),

    # -- proteomics: the two studies added to the arm -----------------------
    Claim("pxd054342_features",
          r"yielded\s+([\d,]+) protein-level effect sizes",
          lambda: len(_prot_study("PXD054342"))),
    Claim("pxd055816_features",
          r"animals per group\) yielded ([\d,]+), of which",
          lambda: len(_prot_study("PXD055816"))),
    Claim("pxd055816_significant",
          r"yielded [\d,]+, of which ([\d,]+)\s+were significant",
          lambda: int((_prot_study("PXD055816")["padj"] < 0.05).sum())),

    # -- proteomics: the peptide-to-precursor mapping -----------------------
    Claim("peptides_mapped",
          r"containment mapped (\d+) of [\d,]+\s+peptides",
          lambda: int((_read(PEPTIDE_MAP, dtype=str)["peptide"].fillna("")
                       .str.strip() != "").sum())),
    Claim("peptides_total",
          r"containment mapped \d+ of ([\d,]+)\s+peptides",
          lambda: int(_read(PROTEOMICS / "PXD013362_effects_collapsed.csv",
                            usecols=["feature_id"])["feature_id"].nunique())),
    Claim("peptides_curated",
          r"recovered (\d+) of UniProt\'s\s+independently",
          lambda: int(_read(PEPTIDE_MAP, dtype=str)["evidence"]
                      .str.startswith("named_cleavage_product").sum())),

    # -- proteomics: the pooled arm -----------------------------------------
    Claim("proteomic_pooled_features",
          r"produced estimates for\s+([\d,]+) features, of which",
          lambda: len(_prot_pooled())),
    Claim("proteomic_pooled_family",
          r"features, of which (\w+) are measured in at least three units",
          lambda: len(_prot_family())),
    Claim("proteomic_pooled_min_padj",
          r"smallest adjusted \\pval-value is\s+([\d.]+),",
          lambda: float(_prot_family()["padj_pooled"].min()),
          tol=0.005,
          note="quoted to two decimals; nothing in the family is significant"),
    Claim("proteomic_pooled_downward",
          r"(\w+) of\s+the six point downward",
          lambda: int((_prot_family()["yi_pooled"] < 0).sum())),
    Claim("proteomic_i2_min",
          r"I\\textsuperscript\{2\} runs from (\d+) to \d+",
          lambda: float(_prot_family()["I2"].min()),
          tol=0.5),
    Claim("proteomic_i2_max",
          r"I\\textsuperscript\{2\} runs from \d+ to (\d+)",
          lambda: float(_prot_family()["I2"].max()),
          tol=0.5),
    Claim("calca_i2",
          r"disagrees in direction across\s+its four units \(I\\textsuperscript\{2\}\\,=\\,(\d+)\)",
          lambda: _prot_gene("Calca", "I2"),
          tol=0.5),

    # Each pooled effect size the sentence names, at the precision it uses.
    Claim("scg2_g", r"Scg2 \(\$g\$\\,=\\,\$-\$([\d.]+)\)",
          lambda: -_prot_gene("Scg2", "yi_pooled"), tol=0.005),
    Claim("pcsk1n_g", r"Pcsk1n \(\$-\$([\d.]+)\)",
          lambda: -_prot_gene("Pcsk1n", "yi_pooled"), tol=0.005),
    Claim("penk_g", r"Penk \(\$-\$([\d.]+)\)",
          lambda: -_prot_gene("Penk", "yi_pooled"), tol=0.005),
    Claim("scg3_g", r"Scg3\s*\n?\(\$-\$([\d.]+)\)",
          lambda: -_prot_gene("Scg3", "yi_pooled"), tol=0.005),
    Claim("calca_g", r"Calca \(\$\+\$([\d.]+)\)",
          lambda: _prot_gene("Calca", "yi_pooled"), tol=0.005),
    Claim("adcyap1_g", r"Adcyap1 \(\$-\$([\d.]+)\)",
          lambda: -_prot_gene("Adcyap1", "yi_pooled"), tol=0.005),

    # -- metabolomics: screening and triage ---------------------------------
    Claim("metab_datasets",
          r"MetaboLights alone returned (\d+) pain-targeted metabolomics",
          lambda: len(_read(PRISMA / "metabolomics_manual_review.csv",
                            usecols=["accession"]))),
    Claim("metab_screened",
          r"of\s+which (\w+) passed screening on species and pain phenotype",
          lambda: int((_read(PRISMA / "metabolomics_manual_review.csv",
                             usecols=["verdict"])["verdict"] == "include").sum())),
    Claim("metab_empty_mafs",
          r"own right: (\w+) deposit a\s+metabolite assignment file",
          lambda: int((_triage_metab()["reason_code"]
                       == "maf_columns_declared_but_empty").sum())),

    # -- metabolomics: Metabolomics Workbench retrieval ---------------------
    Claim("metab_wb_enumerated",
          r"retrieving all ([\d,]+) studies in one",
          lambda: int(_read(PRISMA / "metabolomics_workbench_counts.csv")
                      ["n_enumerated"].iloc[0]), source="supplementary"),
    Claim("metab_wb_candidates",
          r"screen returned (\d+) metabolomic candidates and one lipidomic",
          lambda: int(_read(PRISMA / "metabolomics_workbench_counts.csv")
                      .query("modality == 'metabolomics'")["n_candidates"].iloc[0]),
                      source="supplementary"),
    Claim("metab_wb_included",
          r"(\w+) studies were\s*\n?included\.",
          lambda: int((_read(PRISMA / "metabolomics_workbench_candidates.csv",
                             usecols=["verdict"])["verdict"] == "include").sum()),
                             source="supplementary"),

    # -- metabolomics: the identifier map -----------------------------------
    Claim("metab_map_corpus_tier",
          r"The (\d+)\s*\ncorpus-stated assignments were re-resolved",
          lambda: int(_read(MET_MAP, dtype=str)["evidence"]
                      .str.startswith("metabolights_maf").sum()), source="supplementary"),
    Claim("metab_map_assignments",
          r"This yielded\s*\n([\d,]+) assignments over",
          lambda: int((_read(MET_MAP, dtype=str)["chebi_id"].fillna("") != "").sum()),
          source="supplementary"),
    Claim("metab_map_distinct_ids",
          r"assignments over ([\d,]+) distinct ChEBI identifiers",
          lambda: _read(MET_MAP, dtype=str).query("chebi_id != ''")["chebi_id"].nunique(),
          source="supplementary"),
    Claim("metab_map_ambiguous",
          r"(\w+) names to which two studies assigned different accessions were",
          lambda: int((_read(MET_MAP, dtype=str)["evidence"]
                       .str.startswith("ambiguous_in_corpus")).sum()), source="supplementary"),
    # The cross-check drifted silently once because it was computed by hand and
    # never registered: the manuscript claimed 255 corpus-derived assignments
    # long after full-recall resolution had taken that tier to 724.
    Claim("metab_crosscheck_resolved",
          r"agreed in \d+ of the (\d+) cases where both resolve",
          lambda: int((_read(MET_CROSSCHECK, dtype=str)["ols4_chebi_id"]
                       .fillna("") != "").sum()), source="supplementary"),
    Claim("metab_crosscheck_agree",
          r"agreed in (\d+) of the \d+ cases where both resolve",
          lambda: int((_read(MET_CROSSCHECK, dtype=str)["agrees"] == "True").sum()),
          source="supplementary"),

    # -- metabolomics: the pooled arm ---------------------------------------
    Claim("metab_pooled_features",
          r"gave\s*\nestimates for ([\d,]+) features",
          lambda: len(_metab_pooled())),
    Claim("metab_pooled_family",
          r"of which (\d+) reach three or more",
          lambda: int(_metab_pooled()["padj_pooled"].notna().sum())),
    Claim("metab_pooled_k4",
          r"studies and enter the correction family; (\d+) reach four studies",
          lambda: int((_metab_pooled()["k"] == 4).sum())),
    Claim("metab_pooled_k5",
          r"reach four studies and (\d+) reach\s*\nfive",
          lambda: int((_metab_pooled()["k"] == 5).sum())),
    # The four significant metabolites' adjusted p-values were stated in the
    # manuscript but not registered, so a re-pool moved three of them
    # (2.6->2.7, 3.6->3.7, 2.0->2.1) with nothing to catch it. Compared at the
    # precision the sentence quotes, as the named-gene claims are.
    Claim("metab_padj_cdca",
          r"\\padj\\,=\\,([\d.]+)\$\\times\$10\\textsuperscript\{-10\}",
          lambda: _metab_padj("NAME:chenodeoxycholicacid") * 1e10, tol=0.05),
    Claim("metab_padj_carnitine",
          r"tetradecenoylcarnitine\s*\n?\(\$\+\$0\.59,\s*([\d.]+)\$\\times\$10",
          lambda: _metab_padj("CHEBI:73060") * 1e4, tol=0.05),
    Claim("metab_padj_linoleic",
          r"linoleic acid \(\$\+\$0\.46,\s*\n([\d.]+)\$\\times\$10",
          lambda: _metab_padj("CHEBI:17351") * 1e2, tol=0.05),
    Claim("metab_padj_riboflavin",
          r"riboflavin \(\$-\$0\.34,\s*\n([\d.]+)\$\\times\$10",
          lambda: _metab_padj("CHEBI:17015") * 1e2, tol=0.05),
    # Both counts are stated in Limitations as evidence that name-keyed
    # features remain; they move whenever 05j re-runs.
    Claim("metab_name_keyed_st000676",
          r"contribute (\d+) and \d+ features that key on names",
          lambda: len(_read(MET_PER_STUDY / "ST000676_Dorsal_root_ganglia"
                            / "ST000676_Dorsal_root_ganglia_effects.csv",
                            usecols=["feature_id"])), source="supplementary"),
    Claim("metab_name_keyed_st000780",
          r"contribute \d+ and (\d+) features that key on names",
          lambda: len(_read(MET_PER_STUDY / "ST000780" / "ST000780_effects.csv",
                            usecols=["feature_id"])), source="supplementary"),
    Claim("metab_significant",
          r"(\w+) features are significant, headed by heptadecanoic",
          lambda: int((_metab_pooled()["padj_pooled"] < 0.05).sum())),
    # The three-study pool let one study set the result. Recording both the old
    # and the new share is the point of the sentence, so both are tested.
    Claim("metab_weight_share",
          r"It now carries\s*\n?([\d.]+)\\%",
          lambda: _metab_weight_share("MTBLS2774"), tol=0.1,
          note="median-SE approximation of the inverse-variance weight share; "
               "the tolerance matches the one decimal the sentence quotes, so a "
               "wrong first decimal fails rather than passing inside a wide band",
               source="supplementary"),
    Claim("metab_dominant_se",
          r"median standard error of ([\d.]+) against 0\.55 and 0\.61",
          lambda: float(_metab_study("MTBLS2774")["se"].median()), tol=0.0005,
          source="supplementary"),
]


def _expected_of(name: str) -> Callable[[], float]:
    """The expectation another claim already computes, for a restatement."""
    return next(c.expected for c in CLAIMS if c.name == name)


# The supplementary describes Tables 7 and 9 by their correction families.
# Neither sentence was claimed, and Table 9's said 74 features for three months
# after the metabolomic family had grown to 236 -- found while moving the text
# for Brain Communications (2026-09-21). Each is bound to the same expectation
# as the paper's own statement of that family, so the two cannot disagree.
CLAIMS += [
    Claim("s7_prot_family",
          r"Correction is applied within the\s+(\w+) features measured in at least three units",
          _expected_of("proteomic_pooled_family"), source="supplementary"),
    Claim("s9_metab_family",
          r"correction is applied within the\s+([\d,]+) features reaching three or more units",
          _expected_of("metab_pooled_family"), source="supplementary"),
]

# The manuscript spells small counts as words. The map ran to ten and stopped,
# so "Nineteen features are significant" failed with a bare float() ValueError
# rather than a claim mismatch -- the count had simply outgrown the vocabulary.
_WORDS = {"one": 1, "two": 2, "three": 3, "four": 4, "five": 5, "six": 6,
          "seven": 7, "eight": 8, "nine": 9, "ten": 10, "eleven": 11,
          "twelve": 12, "thirteen": 13, "fourteen": 14, "fifteen": 15,
          "sixteen": 16, "seventeen": 17, "eighteen": 18, "nineteen": 19,
          "twenty": 20}


@cache
def _document(source: str) -> str:
    path = DOCUMENTS[source]
    if not path.exists():
        pytest.skip(f"no {source} at {path}")
    return path.read_text()


def _manuscript() -> str:
    """The manuscript text. Uncached, so that clearing _document is enough.

    Two caches over the same file would let the mutation harness restore one
    and keep serving the other, which is a stale read that looks like a
    passing claim.
    """
    return _document("manuscript")


def _as_number(token: str) -> float:
    token = token.strip()
    if token.lower() in _WORDS:
        return float(_WORDS[token.lower()])
    # The manuscript writes a negative number as `$-$0.002`, LaTeX's math minus
    # rather than a hyphen. Without this the registry cannot bind a negative at
    # all, which is how rat's Pearson r sat unclaimed while it was negative.
    token = token.replace("$-$", "-")
    try:
        return float(token.replace(",", ""))
    except ValueError:
        raise AssertionError(
            f"the manuscript states {token!r} where a number is expected; if "
            "that is a number word, add it to _WORDS"
        ) from None


def test_claim_registry_has_no_duplicate_names():
    names = [c.name for c in CLAIMS]
    assert len(names) == len(set(names)), "duplicate claim names in the registry"


@pytest.mark.parametrize("claim", CLAIMS, ids=lambda c: c.name)
def test_manuscript_claim_matches_results(claim: Claim):
    text = _document(claim.source)

    matches = re.findall(claim.pattern, text)
    assert matches, (
        f"claim '{claim.name}' no longer appears in the {claim.source}.\n"
        f"  pattern: {claim.pattern}\n"
        "  If the sentence was rewritten, update this claim; if the result was "
        "dropped, delete it. Do not leave it unmatched."
    )
    assert len(matches) == 1, (
        f"claim '{claim.name}' matches {len(matches)} places in the "
        f"{claim.source}: {matches}. The same figure stated twice can drift "
        "apart; give each its own claim."
    )

    try:
        expected = float(claim.expected())
    except MissingInput as exc:
        pytest.skip(f"results not generated: {exc}")

    stated = _as_number(matches[0])
    delta = abs(stated - expected)
    assert delta <= claim.tol, (
        f"{claim.name}: the {claim.source} says {stated:,.4g}, results give "
        f"{expected:,.4g} (difference {delta:,.4g}, tolerance {claim.tol})."
        + (f"\n  {claim.note}" if claim.note else "")
    )


#: Brain Communications' abbreviated summary: at most 50 words, third person
#: ("Zippo reports that ..."), accessible to a non-specialist, and free of
#: abbreviations other than approved gene symbols. It is entered in the
#: submission system and reused on social media, so nothing downstream would
#: catch an over-long or jargon-laden one.
_SUMMARY_MAX_WORDS = 50


def test_abbreviated_summary_meets_the_journal_format():
    if not ABBREVIATED_SUMMARY.exists():
        pytest.skip(f"no abbreviated summary at {ABBREVIATED_SUMMARY}")
    text = ABBREVIATED_SUMMARY.read_text().strip()
    words = text.split()
    assert len(words) <= _SUMMARY_MAX_WORDS, (
        f"abbreviated_summary.txt has {len(words)} words; the journal allows "
        f"{_SUMMARY_MAX_WORDS}")
    assert text.startswith("Zippo reports"), (
        "the journal asks for the third person: 'Author et al. report that ...'")
    abbreviations = re.findall(r"\b[A-Z]{2,}\b", text)
    assert not abbreviations, f"abbreviations the journal disallows: {abbreviations}"
    # A number here would be a copy of a manuscript figure that no claim reads,
    # which is how the Highlights had to be bound one by one. Stating none keeps
    # the summary outside the registry honestly rather than by omission.
    assert not re.search(r"\d", text), (
        "the abbreviated summary states a number; bind it with a claim or drop it")


SUPPLEMENTARY_DIR = REPO_ROOT / "manuscript" / "supplementary"


#: "Supplementary Table~5", "Supplementary Tables~6 and~8", "Supplementary
#: Tables~1, 2, 6 and~8", "Supplementary Tables~1--18": the journal's form,
#: which dropped the S prefix the Elsevier version carried.
_SUPP_TABLE_REF = re.compile(
    r"Supplementary\s+Tables?[~\s]+"
    r"((?:\d+(?:--\d+)?)(?:(?:,\s*|,?\s+and[~\s]+)\d+(?:--\d+)?)*)")


def _cited_supplementary_numbers() -> set[int]:
    """Every Supplementary Table number the two documents reference."""
    text = _manuscript() + _document("supplementary")
    cited: set[int] = set()
    for group in _SUPP_TABLE_REF.findall(text):
        for lo, hi in re.findall(r"(\d+)(?:--(\d+))?", group):
            cited.update(range(int(lo), int(hi or lo) + 1))
    return cited


def test_the_title_still_describes_the_results():
    """The title asserts non-overlapping gene sets. Nothing read it.

    Every other assertion in this paper is bound to a number the registry
    recomputes; the title is the one sentence an editor reads first and the
    only load-bearing claim that stated no figure at all. It does state a
    testable condition: that the evidence layers implicate *non-overlapping*
    sets. Two overlaps would falsify it, and both are computable.

    This is deliberately a floor, not a paraphrase of the title. It fails when
    the word "non-overlapping" stops being true, not when the wording changes.
    """
    title = re.search(r"\\title\{(.+?)\}\}", _manuscript(), re.S)
    assert title, "the title is no longer recognisable; update this test"
    # "Disjoint evidence layers" since the Brain Communications title; the
    # claim is the same one "non-overlapping gene sets" made.
    if not re.search(r"non-overlapping|disjoint", title.group(1), re.I):
        pytest.skip("the title no longer claims disjoint gene sets")

    # (1) Species: no ortholog pair significant on both the rodent and the
    # human side, in either species.
    both = 0
    for slug in ("musculus", "norvegicus"):
        try:
            d = _concordance(slug)
        except MissingInput as exc:
            pytest.skip(f"results not generated: {exc}")
        both += int(d["both_significant"].sum())
    assert both == 0, (
        f"{both} ortholog pair(s) are now significant in both species, so the "
        "title's 'non-overlapping' no longer holds for the species axis")

    # (2) Genomic risk: no genome-wide-significant MAGMA gene is significant in
    # the human transcriptomic pool.
    try:
        magma = set(_read(GENOMICS / "magma" / "significant_genes.csv",
                          usecols=["gene"])["gene"].astype(str))
        human = _species("Homo_sapiens")
    except MissingInput as exc:
        pytest.skip(f"results not generated: {exc}")
    sig_human = set(human.loc[human["padj"] < 0.05, "feature_id"].astype(str))
    overlap = sorted(magma & sig_human)
    assert not overlap, (
        f"MAGMA gene(s) {', '.join(overlap)} are now significant in the human "
        "transcriptomic pool, so the title's 'disjoint' no longer holds "
        "for the genomic axis")


def test_the_data_availability_range_covers_every_table():
    """The paper states its supplementary range as S1--SN in one place.

    It said S1--S10 after seven more tables had been added, and nothing read
    it: the per-number checks below bind each table individually, and a range
    is not a number any of them sees.
    """
    built = {int(p.name.split("_")[1][1:])
             for p in SUPPLEMENTARY_DIR.glob("Table_S*_*.xlsx")}
    if not built:
        pytest.skip("no workbooks built; run scripts/build_supplementary_xlsx.py")
    m = re.search(r"Supplementary Tables~(\d+)--(\d+)", _manuscript())
    assert m, "the data-availability sentence no longer states a table range"
    lo, hi = int(m.group(1)), int(m.group(2))
    assert (lo, hi) == (min(built), max(built)), (
        f"the paper states {lo}--{hi} but S{min(built)}--S{max(built)} ship")


def test_every_cited_supplementary_table_ships_as_a_workbook():
    """A number in the text with no file beside it is a dead reference.

    The paper used to cite its underlying data by repository path, which a
    reader with the PDF cannot open. Those citations are now table numbers, and
    a number is only useful if the journal receives the file. Generate them
    with scripts/build_supplementary_xlsx.py.
    """
    cited = _cited_supplementary_numbers()
    assert cited, "no Supplementary Table references found; did the wording change?"
    missing = sorted(n for n in cited
                     if not list(SUPPLEMENTARY_DIR.glob(f"Table_S{n}_*.xlsx")))
    assert not missing, (
        "cited in the text but not shipped as .xlsx: "
        + ", ".join(f"S{n}" for n in missing)
        + "\n  run: python scripts/build_supplementary_xlsx.py")


def test_every_supplementary_workbook_is_cited():
    """And the converse: a workbook nobody references is one nobody will open."""
    built = {int(p.name.split("_")[1][1:])
             for p in SUPPLEMENTARY_DIR.glob("Table_S*_*.xlsx")}
    if not built:
        pytest.skip("no workbooks built; run scripts/build_supplementary_xlsx.py")
    orphaned = sorted(built - _cited_supplementary_numbers())
    assert not orphaned, (
        "shipped but referenced nowhere in the text: "
        + ", ".join(f"S{n}" for n in orphaned))


# ---------------------------------------------------------------------------
# Named-gene statistics
# ---------------------------------------------------------------------------
# The manuscript quotes a p-value for each gene it names. These are rounded for
# print, so the assertion is that the computed value, rounded to the precision
# the sentence uses, is what the sentence says -- 7.09876e-41 written as
# 7.1e-41 passes, written as 7.2e-41 does not.


@cache
def _magma() -> pd.DataFrame:
    return _read(GENOMICS / "magma" / "gene_results.csv",
                 usecols=["group", "gene", "pval"])


# significant_genes.csv keys genes by Entrez id, so a symbol the manuscript
# names has to be resolved. Only the genes actually quoted are listed.
ENTREZ = {"MLLT10": "8028"}


def _magma_pval(symbol: str, group: str) -> float:
    df = _magma()
    eid = ENTREZ[symbol]
    row = df[(df["gene"].astype(str) == eid) & (df["group"] == group)]
    if row.empty:
        raise MissingInput(f"MAGMA has no {symbol} (Entrez {eid}) in {group}")
    return float(row["pval"].iloc[0])


def _gene_stat(loader: Callable[[], pd.DataFrame], gene: str, col: str) -> float:
    df = loader()
    row = df[df["feature_id"].astype(str) == gene]
    if row.empty:
        raise MissingInput(f"{gene} absent from the table")
    return float(row[col].iloc[0])


@dataclass(frozen=True)
class GeneClaim:
    """A p-value the manuscript attributes to a named gene."""

    name: str
    pattern: str                      # captures (mantissa, exponent) or (decimal,)
    expected: Callable[[], float]
    bound: bool = False               # pattern states an upper bound, not a value


_SCI = r"([\d.]+)\$\\times\$10\\textsuperscript\{(-?\d+)\}"

GENE_CLAIMS: list[GeneClaim] = [
    # SNL stratum
    # Cd36 is no longer the top hit and is quoted as a decimal; the earlier
    # 7.1e-41 was probe-selection bias, not a result.
    GeneClaim("snl_cd36", r"\\gene\{Cd36\} \(k=5, \\padj\\,=\\,([\d.]+)\)",
              lambda: _gene_stat(lambda: _stratum("SNL"), "Cd36", "padj")),
    GeneClaim("snl_h19", r"long non-coding RNA \(\\padj\\,=\\,([\d.]+)\)",
              lambda: _gene_stat(lambda: _stratum("SNL"), "H19", "padj")),
    GeneClaim("snl_bdnf", r"\(\\gene\{Bdnf\}; \\padj\\,=\\,([\d.]+)\)",
              lambda: _gene_stat(lambda: _stratum("SNL"), "Bdnf", "padj")),
    # CCI stratum
    # cci_sting1 deleted: Sting1 pools at padj 0.33 and the STING-pathway
    # sentence it verified has been removed from the manuscript.
    GeneClaim("cci_ctsb", r"\\gene\{Ctsb\}\s*\n?\(\\padj\\,=" + _SCI,
              lambda: _gene_stat(lambda: _stratum("CCI"), "Ctsb", "padj")),
    GeneClaim("cci_lgmn", r"\\gene\{Lgmn\}\s*\n?\(\\padj\\,=\\,([\d.]+)\)",
              lambda: _gene_stat(lambda: _stratum("CCI"), "Lgmn", "padj")),
    GeneClaim("cci_rt1dmb", r"\\gene\{RT1-DMb\}\s*\n?\(\\padj\\,=" + _SCI,
              lambda: _gene_stat(lambda: _stratum("CCI"), "RT1-DMb", "padj")),
    GeneClaim("cci_nfkbid", r"\\gene\{Nfkbid\}\s*\n?\(\\padj\\,=" + _SCI,
              lambda: _gene_stat(lambda: _stratum("CCI"), "Nfkbid", "padj")),
    GeneClaim("cci_jun", r"\\gene\{Jun\}\s*\n?\(\\padj\\,=" + _SCI,
              lambda: _gene_stat(lambda: _stratum("CCI"), "Jun", "padj")),
    # The retired microglial reading. Quoted so that if C1qa ever becomes
    # tissue-dependent again, the sentence denying it fails.
    GeneClaim("sni_tissue_c1qa",
              r"those [\d,]+ and are null \(\\padj\\,=\\,([\d.]+)\)",
              lambda: _gene_stat(_sni_tissue, "C1qa", "padj")),
    # CIPN stratum
    GeneClaim("cipn_bpifb5", r"headed by \\gene\{Bpifb5\} \(\\padj\\,=" + _SCI,
              lambda: _gene_stat(lambda: _stratum("CIPN"), "Bpifb5", "padj")),
    GeneClaim("cipn_itga5", r"\\gene\{Itga5\} \(\\padj\\,=" + _SCI,
              lambda: _gene_stat(lambda: _stratum("CIPN"), "Itga5", "padj")),
    # human pool
    GeneClaim("human_closest_magma",
              r"the closest, \\gene\{TMEM134\}, stands at \\padj\\,=\\,(\d+\.\d+)",
              lambda: _gene_stat(lambda: _species("Homo_sapiens"), "TMEM134", "padj")),
    # genomics
    GeneClaim("magma_mllt10", r"\\gene\{MLLT10\} \(the strongest gene, \\pval\\,=" + _SCI,
              lambda: _magma_pval("MLLT10", "back_pain")),
    GeneClaim("gwas_lrp1", r"\\gene\{LRP1\} \(10 cohorts, minimum\s*\n?\\pval\\,=" + _SCI,
              lambda: float(_read(GENOMICS / "gwas_gene_meta_replication.csv")
                            .set_index("gene").loc["LRP1", "min_p"])),
    # stratified table
    # The study count is not pinned in these patterns. It was, and when the
    # human pool went from 24 studies to 23 the same mistake in _sens_claims
    # unmatched six claims at once -- reported as "no longer appears in the
    # manuscript" rather than as a wrong number. The count has its own claim.
    GeneClaim("table_rt1dmb",
              r"CCI & \d+ & [\d,]+ & [\d,]+ & \d+ & \\gene\{RT1-DMb\} \(([\d.]+)e(-\d+)\)",
              lambda: _gene_stat(lambda: _stratum("CCI"), "RT1-DMb", "padj")),
    GeneClaim("table_exoc8",
              r"SNL & \d+ & [\d,]+ & [\d,]+ & [\d,]+ & \\gene\{Exoc8\} \(([\d.]+)e(-\d+)\)",
              lambda: _gene_stat(lambda: _stratum("SNL"), "Exoc8", "padj")),
    GeneClaim("table_clec2d2l1", r"\\gene\{Clec2d2l1\} \(([\d.]+)e(-\d+)\)",
              lambda: _gene_stat(lambda: _stratum("SNI"), "Clec2d2l1", "padj")),
    # Stated as a value, not as the "$<$1e-300" bound it replaced: a bound
    # claim on a padj that underflows to exactly 0 is satisfied by anything,
    # which is why table_gm35167 was the one claim the mutation harness could
    # never kill.
    GeneClaim("table_bpifb5",
              r"CIPN & \d+ & [\d,]+ & [\d,]+ & \d+ & \\gene\{Bpifb5\} \(([\d.]+)e(-\d+)\)",
              lambda: _gene_stat(lambda: _stratum("CIPN"), "Bpifb5", "padj")),
]


def _significant_digits(mantissa: str) -> int:
    digits = mantissa.replace(".", "").lstrip("0")
    return len(digits) or 1


def _round_sig(x: float, digits: int) -> float:
    # NaN reaches here when the named feature has left the corrected family --
    # it fell below meta.min_studies, so the results carry padj = NA rather
    # than a worse number. That is a real outcome and the claim must fail
    # saying so; before this it raised inside math.floor and took the whole
    # mutation harness down with it (2026-08-29, table_fam150b/table_gm35167).
    if math.isnan(x):
        raise AssertionError(
            "the results give NA for this feature: it is no longer measured in "
            "enough studies to enter the BH family, so the manuscript cannot "
            "quote a corrected p-value for it")
    if x == 0:
        return 0.0
    return round(x, -int(math.floor(math.log10(abs(x)))) + (digits - 1))


def test_gene_claim_registry_has_no_duplicate_names():
    names = [c.name for c in GENE_CLAIMS]
    assert len(names) == len(set(names)), "duplicate gene-claim names"


@pytest.mark.parametrize("claim", GENE_CLAIMS, ids=lambda c: c.name)
def test_named_gene_statistic_matches_results(claim: GeneClaim):
    text = _manuscript()

    matches = re.findall(claim.pattern, text)
    assert matches, (
        f"gene claim '{claim.name}' no longer appears in the manuscript.\n"
        f"  pattern: {claim.pattern}\n"
        "  Update the claim if the sentence was rewritten, or delete it if the "
        "gene is no longer named."
    )
    assert len(matches) == 1, (
        f"gene claim '{claim.name}' matches {len(matches)} places: {matches}"
    )

    try:
        expected = float(claim.expected())
    except MissingInput as exc:
        pytest.skip(f"results not generated: {exc}")

    m = matches[0]
    if claim.bound:
        exponent = abs(int(m if isinstance(m, str) else m[0]))
        # 10**-400 underflows to 0.0, so compare in log space. A padj that
        # underflowed to exactly 0 satisfies any bound and cannot be used to
        # check the exponent -- the test confirms the direction, not the digit.
        if expected == 0:
            return
        assert not math.isnan(expected), (
            f"{claim.name}: the results give NA for this feature -- it no longer "
            f"enters the BH family, so no bound can be asserted for it."
        )
        assert math.log10(expected) < -exponent, (
            f"{claim.name}: manuscript claims padj < 1e-{exponent}, "
            f"results give {expected:.3g}."
        )
        return

    if isinstance(m, str):                       # plain decimal, e.g. 0.074
        mantissa, stated = m, float(m)
    else:                                        # LaTeX scientific notation
        mantissa, exponent = m
        stated = float(f"{mantissa}e{exponent}")

    digits = _significant_digits(mantissa)
    rounded = _round_sig(expected, digits)
    # abs=0 is essential. pytest.approx keeps a default absolute tolerance of
    # 1e-12 and passes when *either* tolerance is met, so with the default
    # every p-value below 1e-12 compares equal to every other and the
    # assertion becomes vacuous for exactly the values this test exists to
    # check.
    assert rounded == pytest.approx(stated, rel=1e-9, abs=0), (
        f"{claim.name}: manuscript states {stated:.3g}, results give "
        f"{expected:.6g}, which to {digits} significant digits is "
        f"{rounded:.3g}."
    )


def test_round_sig_reports_a_feature_that_left_the_corrected_family():
    # A named gene can fall below meta.min_studies between runs, leaving padj
    # as NA. That must fail the claim with an explanation rather than raising
    # ValueError out of math.floor, which is what stopped the mutation harness
    # from running at all on 2026-08-29.
    with pytest.raises(AssertionError, match="no longer measured in enough studies"):
        _round_sig(float("nan"), 2)


def test_no_ortholog_pair_is_significant_in_both_species():
    """The Results assert that not one pair clears padj<0.05 on both sides.

    A count of zero cannot be written as a captured number the way every other
    claim is, and a bound ("fewer than one") is satisfied by anything, which is
    the flaw that made table_gm35167 unkillable. So it is asserted directly.
    Cfap68 was the single rat pair before the 2026-08-30 re-run; if any pair
    comes back, this fails and the sentence has to change.
    """
    try:
        counts = {slug: int(_concordance(slug)["both_significant"].sum())
                  for slug in ("musculus", "norvegicus")}
    except MissingInput as exc:
        pytest.skip(f"results not generated: {exc}")

    for slug, n in counts.items():
        assert n == 0, (
            f"{slug}: {n} ortholog pair(s) are now significant in both species, "
            "so the Results sentence claiming none is wrong"
        )


def test_forest_caption_names_the_top_pooled_feature():
    """Figure 2's caption calls its subject "the most significant pooled feature".

    That is a ranking, not a number, so no Claim pattern covers it -- and it
    drifted: the caption still named Olr830 (padj 0.649 in the current pool)
    after 08_figures.py had regenerated the figure around SERTM2.
    """
    try:
        pooled = _read(TX_META / "pooled_effects.csv",
                       usecols=["feature_id", "padj_pooled"])
    except MissingInput as exc:
        pytest.skip(f"results not generated: {exc}")
    fam = pooled[pooled["padj_pooled"].notna()]
    top = str(fam.sort_values("padj_pooled")["feature_id"].iloc[0])
    caption = re.search(r"\\caption\{(?:\\textbf\{)?Random-effects forest plot.*?\}",
                        _manuscript(), re.S)
    assert caption, "the forest figure's caption is no longer recognisable"
    assert top in caption.group(0), (
        f"the caption names a different feature than the top of the pool ({top})"
    )


def test_between_model_sentence_names_the_top_q_features():
    """The Results name the five features leading the between-model Q test.

    A list of gene names is a ranking, not a number, so no Claim covers it --
    the same gap the forest caption fell through. It had drifted: the sentence
    named Alox12e, Ambp, Apoa2, Asic5 and Clec4m, of which only two are still
    in the top five, and Alox12e stood at the head of a table row whose own
    counts were two regenerations stale.

    Ranked on padj then raw p, as Table 2's rows are: the leaders all underflow
    to padj exactly 0, so padj alone leaves the order arbitrary.
    """
    try:
        q = _between_model_q()
    except MissingInput as exc:
        pytest.skip(f"results not generated: {exc}")
    top5 = list(q.sort_values(["padj_QM", "QM_pval"])["feature_id"].head(5))
    sentence = re.search(
        r"The between-model Q test identified.*?heterogeneity", _manuscript(), re.S)
    assert sentence, "the between-model Q sentence is no longer recognisable"
    missing = [g for g in top5 if g not in sentence.group(0)]
    assert not missing, (
        f"the sentence does not name {missing}, which now lead the Q test "
        f"(top five: {top5})"
    )


# ---------------------------------------------------------------------------
# Identifier space: is every pooled study actually keyed by gene symbols?
#
# This is the defect that keeps returning by a different door. The 2026-08-16
# repair removed unmapped array addresses from the pools; the 2026-08-29
# ensure_gpl_soft() fix removed fifteen studies written out in probe space.
# Both were found by audit, because "could not obtain ID map -- keeping native
# IDs" is also what a study whose ids are already symbols prints.
#
# A study in an unmappable id space is not merely untidy. Its features join
# nothing, so they inflate the feature total and the below-threshold count
# while contributing nothing to the correction family: on 2026-08-29 three such
# studies supplied 215,320 of the 371,856 pooled features (58%) and 0 of the
# 78,836 in the BH family.
# ---------------------------------------------------------------------------

#: Identifier forms that are addresses or accessions, never gene symbols.
_UNMAPPED_ID = re.compile(
    r"^(ENS[A-Z]*\d{6,}"        # Ensembl gene/transcript accessions
    r"|ASHGV\d+"                # Arraystar array addresses
    r"|A_\d+_P\d+"              # Agilent probe names
    r"|\(\+\)|\(-\)"           # Agilent spike-in controls, e.g. (+)E1A_r60_1
    r"|[NX][MR]_\d"             # RefSeq accessions
    r"|\d+_[sx]?_?at$|\d+_st$)"  # Affymetrix probeset addresses
)


def test_no_pooled_study_sits_in_an_unmapped_identifier_space():
    """Every study the meta-analysis pools must be keyed by gene symbols.

    Judged per study rather than per feature: a symbol table carrying a few
    stray accessions is normal (GSE221921 is 99.8% symbols with 42 ENSG rows),
    whereas a study that is *predominantly* accessions was never mapped at all.
    """
    per_study = REPO_ROOT / "results" / "per_study" / "transcriptomics"
    tables = sorted(per_study.glob("*_effects_normalized.csv"))
    if not tables:
        pytest.skip("no normalized per-study tables; run pipeline/05b")

    offenders = {}
    for path in tables:
        study = path.name.replace("_effects_normalized.csv", "")
        ids = _read(path, usecols=["feature_id"], dtype=str)["feature_id"].astype(str)
        frac = float(ids.str.match(_UNMAPPED_ID).mean())
        if frac > 0.5:
            offenders[study] = (frac, len(ids))

    assert not offenders, (
        "these studies are pooled in an identifier space that is not gene "
        "symbols, so their features can never join another study's:\n"
        + "\n".join(f"  {s}: {f:.1%} of {n:,} rows"
                     for s, (f, n) in sorted(offenders.items()))
        + "\nEither resolve them in pipeline/05b or hold them in "
          "literature/prisma/transcriptomics_candidates.csv and delete their "
          "effect tables -- 06 globs the directory, so a table on disk pools "
          "regardless of its PRISMA verdict."
    )


# ---------------------------------------------------------------------------
# Freshness: do the results files describe the corpus that is screened in?
#
# Every claim above compares the manuscript against a results file. None of
# them asks whether that file still describes the current corpus, so a file
# that is a whole re-run behind passes every claim that reads it, indefinitely.
# That is not hypothetical: 06b was never re-run during the 2026-08-28 chain
# (05 -> 05b -> 06 -> 06g -> 07), so the stratified tables described a corpus
# that had stopped existing, the manuscript was verified against them, and the
# suite reported green. SNL's Features column kept its pre-2026-08-16 value of
# 287,084 through two such runs.
#
# These tests close that gap from the other side: the results must agree with
# literature/prisma about which studies exist.
# ---------------------------------------------------------------------------

TX_CANDIDATES = PRISMA / "transcriptomics_candidates.csv"
SUPERSEDED = REPO_ROOT / "conf" / "analysis" / "superseded_studies.csv"


@cache
def _included_studies() -> set[str]:
    d = _read(TX_CANDIDATES)
    return set(d.loc[d["verdict"] == "include", "accession"].astype(str))


@cache
def _superseded_studies() -> set[str]:
    d = _read(SUPERSEDED)
    return set(d["study_id"].astype(str))


def _parent(study: str) -> str:
    """Derived ids (GSE241361_DRG) inherit their parent's screening verdict.

    The token may carry a digit: `GSE180627_S1` is the primary somatosensory
    cortex, not an accession. Requiring it to be purely alphabetic left that
    one unit of the four unrecognised as derived, so it matched neither the
    screening record nor the supersede replacements and reported as a study
    the corpus had dropped. The token must still *begin* with a letter, so a
    bare accession is never truncated.
    """
    return re.sub(r"_[A-Za-z][A-Za-z0-9_]*$", "", study)


def _studies_in(path: Path) -> set[str]:
    d = _read(path, usecols=["study_ids"])
    return set(d["study_ids"].astype(str).str.split(";").explode().str.strip()) - {""}


@cache
def _supersede_replacements() -> set[str]:
    """Ids the registry promotes in place of a superseded study.

    GSE241361 is superseded *by* GSE241361_DRG and GSE241361_Spinal_cord, which
    are therefore poolable even though they never appear in the screening
    record -- only their parent accession does.
    """
    d = _read(SUPERSEDED)
    named = set(
        d["superseded_by"].astype(str).str.split(";").explode().str.strip()
    ) - {""}
    # Only derived ids qualify. A replacement that appears in the screening
    # record has a verdict of its own and must satisfy it -- GSE57733 replaces
    # GSE57832 and is nonetheless held, for an unrelated reason. And the
    # registry spans every modality, so metabolomic replacements
    # (ST000676_Dorsal_root_ganglia) are not expected in a transcriptomic pool:
    # require the parent accession to be in this modality's record.
    known = set(_read(TX_CANDIDATES)["accession"].astype(str))
    return {r for r in named if r not in known and _parent(r) in known}


def _poolable() -> set[str]:
    return (_included_studies() - _superseded_studies()) | _supersede_replacements()


def test_pooled_table_contains_no_study_that_screening_excludes():
    """A study dropped from the corpus must not survive in a results file.

    This is the direction that caught GSE102937: superseded on 2026-08-29 for
    duplicating GSE102721, it would otherwise have kept contributing to any
    results file written before that decision.
    """
    try:
        found = _studies_in(TX_META / "pooled_effects.csv")
        poolable = _poolable()
    except MissingInput as exc:
        pytest.skip(f"results not generated: {exc}")
    stale = {s for s in found
             if s not in poolable and _parent(s) not in poolable}
    assert not stale, (
        f"pooled_effects.csv still contains {len(stale)} study(ies) that are not "
        f"screened in, or that are superseded: {sorted(stale)[:8]}. The file "
        f"predates a screening decision and every claim reading it is verifying "
        f"the manuscript against a corpus that no longer exists."
    )


def test_pooled_table_contains_every_study_screening_admits():
    """And the other direction: a study admitted since the last run must appear.

    This is what would have caught the stratified tables in August, and what
    caught four studies here -- GSE63442, GSE89224, GSE91396 and GSE92718
    produced no normalized table on 2026-08-28 because their annotation
    packages were absent, and nothing noticed.
    """
    try:
        found = _studies_in(TX_META / "pooled_effects.csv")
        found |= {_parent(s) for s in found}
        poolable = _poolable()
    except MissingInput as exc:
        pytest.skip(f"results not generated: {exc}")
    missing = poolable - found
    assert not missing, (
        f"{len(missing)} study(ies) are screened in but absent from "
        f"pooled_effects.csv: {sorted(missing)[:8]}. Either the meta-analysis "
        f"has not been re-run since they were admitted, or they failed "
        f"normalisation silently."
    )


@pytest.mark.parametrize("stratum_file", sorted(STRAT.glob("*_pooled.csv")) or [None],
                         ids=lambda p: p.stem if p is not None else "no-strata")
def test_stratified_tables_contain_no_excluded_study(stratum_file):
    """Each stratum is a subset of the corpus, so the same rule applies to it.

    Run per file rather than over their union: a single stale stratum is
    exactly the failure this exists to catch, and a union test would let it
    hide behind the fresh ones.
    """
    if stratum_file is None:
        pytest.skip("no stratified tables generated")
    try:
        found = _studies_in(stratum_file)
        poolable = _poolable()
    except MissingInput as exc:
        pytest.skip(f"results not generated: {exc}")
    stale = {s for s in found
             if s not in poolable and _parent(s) not in poolable}
    assert not stale, (
        f"{stratum_file.name} contains {len(stale)} study(ies) not screened in "
        f"or superseded: {sorted(stale)[:8]}. This file predates a screening "
        f"decision; 06b was last re-run before it."
    )


@pytest.mark.parametrize("species", ["Mus musculus", "Rattus norvegicus"])
def test_sensitivity_arms_share_the_main_arm_rodent_pool(species):
    """Table 4's three arms must differ only in which human studies are pooled.

    That is the whole claim the table makes: the rodent pools are identical
    across the arms, so a difference between arms is about tissue compartment
    and not about pooling depth. Nothing enforced it, and the arms are written
    by a separate driver from the main analysis -- so when the corpus grew from
    78 study units to 82 on 2026-08-30, the main arm was regenerated and the
    two sensitivity arms were not. Their claims stayed green because the
    manuscript was stale in exactly the same way the files were, which is the
    one failure mode a claim registry cannot catch on its own.

    Compared on the rodent feature count rather than the pooled estimates: it
    is the cheapest field that moves whenever the rodent corpus does.
    """
    try:
        main = int(_sens_row("all_human", species)["n_animal_features"])
        arms = {slug: int(_sens_row(slug, species)["n_animal_features"])
                for slug in ("tissue_matched", "published_pool")}
    except MissingInput as exc:
        pytest.skip(f"results not generated: {exc}")
    stale = {slug: n for slug, n in arms.items() if n != main}
    assert not stale, (
        f"the main arm pools {main:,} {species} features but "
        + ", ".join(f"{slug} pools {n:,}" for slug, n in sorted(stale.items()))
        + ". The sensitivity arms were not re-run with the corpus; re-run "
          "06g/07 per arm (temp/driver_sens.sh) before trusting Table 4."
    )


CONTRAST_AUDIT = PRISMA / "transcriptomics_contrast_audit.csv"


def test_every_pooled_study_has_been_contrast_audited():
    """A study cannot enter the pool without its arms having been looked at.

    `groups.csv` declares which samples are cases and nothing downstream
    questions it, so a study split on drug treatment, on microarray channel or
    on timepoint yields a well-formed effect table and reports that as a pain
    effect. Two audits found eleven such studies between them
    (plan/2026-08-30-contrast-audit.md), and both were run by hand over the
    corpus as it stood on the day.

    This is the standing form of that check: not pass/fail on the flags --
    whether a split is wrong is a reading of the study, which is why
    scripts/audit_contrasts.py emits a review sheet and not a verdict -- but
    fail when a study reaches the pool that no sheet covers. Derived units
    inherit their parent's row, the arms being curated per accession.
    """
    try:
        pooled = _studies_in(TX_META / "pooled_effects.csv")
        audited = set(_read(CONTRAST_AUDIT, usecols=["accession"])["accession"]
                      .astype(str))
    except MissingInput as exc:
        pytest.skip(f"results not generated: {exc}")

    unreviewed = {s for s in pooled
                  if s not in audited and _parent(s) not in audited}
    assert not unreviewed, (
        f"{len(unreviewed)} pooled study(ies) have no row in "
        f"{CONTRAST_AUDIT.name}: {sorted(unreviewed)[:8]}. Run "
        "`python scripts/audit_contrasts.py` and review the flags before "
        "these contribute to a pooled estimate."
    )


@cache
def _literature_screening() -> pd.DataFrame:
    return _read(PRISMA / "screening.csv",
                 usecols=["source", "screening_decision", "exclusion_reason"])


def _lit(source: str | None = None, decision: str | None = None,
         reason: str | None = None) -> int:
    d = _literature_screening()
    if source is not None:
        d = d[d["source"] == source]
    if decision == "none":
        d = d[d["screening_decision"].isna()]
    elif decision is not None:
        d = d[d["screening_decision"] == decision]
    if reason is not None:
        d = d[d["exclusion_reason"] == reason]
    return len(d)


# The bibliographic search, reported in the supplementary since the flow
# diagram moved to the repositories (2026-09-21). The 921 are the records the
# automated step passed and no one screened; stating them is the point of the
# paragraph, so they are bound like any other result.
CLAIMS += [
    Claim("lit_auto_excluded", r"vocabularies, which excluded ([\d,]+)\s+records",
          lambda: _lit(decision="exclude"), source="supplementary"),
    Claim("lit_no_omics", r"records \(([\d,]+) naming no omics technology",
          lambda: _lit(reason="no_omics_term"), source="supplementary"),
    Claim("lit_not_pain", r"technology and (\d+) no pain phenotype\)",
          lambda: _lit(reason="not_pain_phenotype"), source="supplementary"),
    Claim("lit_pubmed_retained", r"then retained (\d+) PubMed records",
          lambda: _lit(source="pubmed", decision="include"), source="supplementary"),
    Claim("lit_added", r"on reading them, and (\d+) more publications",
          lambda: _lit(decision="include") - _lit(source="pubmed", decision="include"),
          source="supplementary"),
    Claim("lit_unscreened", r"The remaining (\d+)\s+records that passed",
          lambda: _lit(decision="none"), source="supplementary"),
]


def _restate(name: str, of: str, pattern: str, source: str) -> Claim:
    """A claim restating another's number elsewhere: same expectation, same tol."""
    base = next(c for c in CLAIMS if c.name == of)
    return Claim(name, pattern, base.expected, tol=base.tol, source=source)


# The cover letter's numbers, each bound to the manuscript claim it restates.
# The Journal of Pain letter was not bound, and its "354 checks" was already
# wrong by the time the paper moved journals.
CLAIMS += [
    _restate("cl_units", "studies_total_abstract",
             r"It pools (\d+) transcriptomic study units", "cover_letter"),
    _restate("cl_q_testable", "q_testable", r"Of ([\d,]+) features\s+testable", "cover_letter"),
    _restate("cl_q_significant", "q_significant",
             r"at least two models, ([\d,]+) differ significantly", "cover_letter"),
    _restate("cl_null_pct", "q_null_pct", r"([\d.]+)% of features still appear", "cover_letter"),
    _restate("cl_null_p", "q_null_p", r"\(P = ([\d.]+)\)", "cover_letter"),
    _restate("cl_conc_mouse", "concordance_mouse", r"effects is ([\d.]+)% for\s+mouse",
             "cover_letter"),
    _restate("cl_conc_rat", "concordance_rat", r"mouse and ([\d.]+)% for rat, and falls",
             "cover_letter"),
    _restate("cl_tm_mouse", "tm_pct_mouse", r"falls to ([\d.]+)% and [\d.]+% once",
             "cover_letter"),
    _restate("cl_tm_rat", "tm_pct_rat", r"falls to [\d.]+% and ([\d.]+)% once", "cover_letter"),
    _restate("cl_hvh", "hvh_replicate_pct_abstract", r"agree ([\d.]+)% of the time",
             "cover_letter"),
    _restate("cl_magma", "magma_genes", r"None of\s+the (\d+) genome-wide-significant genes",
             "cover_letter"),
    # The registry's own size, counted when the tests run, so it includes
    # these claims and moves whenever one is added or retired.
    Claim("cl_registry_size", r"each of those\s+([\d,]+) tests",
          lambda: len(CLAIMS) + len(GENE_CLAIMS), source="cover_letter"),
]


def test_every_checklist_location_names_a_real_heading():
    """The PRISMA checklist cites sections by name, in italics.

    A renamed heading would leave the checklist pointing at nothing, which is
    how the Journal of Pain checklist's page numbers went stale. Every
    italicised location must be a heading of the paper or the supplementary.
    """
    if not PRISMA_CHECKLIST.exists():
        pytest.skip("no PRISMA checklist")
    headings = {"Abstract"}
    for doc in (_manuscript(), _document("supplementary")):
        headings |= set(re.findall(r"\\(?:sub)*section\*?\{([^}]*)\}", doc))
        headings |= set(re.findall(r"\\paragraph\{([^}]*?)\.?\}", doc))
    rows = [ln.split("|") for ln in PRISMA_CHECKLIST.read_text().splitlines()
            if ln.startswith("|") and not ln.startswith("|---")]
    cited = {m for row in rows if len(row) > 4 for m in re.findall(r"\*([^*|]+)\*", row[4])}
    missing = sorted(c for c in cited if c not in headings)
    assert cited, "no italicised locations found; did the checklist format change?"
    assert not missing, f"the checklist cites headings that do not exist: {missing}"

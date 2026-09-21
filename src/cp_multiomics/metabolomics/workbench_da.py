"""Contrast declarations and arm selection for the Workbench studies.

Contrasts are declared here rather than inferred, on the same reasoning as
`pipeline/05j_metabolomics_da.py`: choosing a comparator is a design
judgement, and a judgement that is reviewable belongs in a table someone can
read. Each entry below records what it compares and why the other arms are not
folded into either side.

Three shapes are needed that the MetaboLights specs did not require.

**Patterns, not only level lists.** ST003954 groups its controls (`HC`, 10
samples) but labels each case individually — `IBS-C-1` through `IBS-C-61`. A
list of exact levels cannot express that, so a level spec may instead be
`"re:<regex>"`, matched against the normalised level.

**Splitting by tissue.** ST000676 measures dorsal root ganglia, sciatic nerve
and plasma; ST003984 measures plasma and peritoneal fluid. Pooling tissues
inside one study would average unrelated compartments, and entering the study
once per tissue without splitting would count its samples several times. Both
therefore split into derived study ids — `ST000676_Dorsal_root_ganglia` and so
on — and the bare accession is registered in
`conf/analysis/superseded_studies.csv`, the mechanism GSE241361 already uses.

**Arms excluded rather than assigned.** ST000676 carries diet-reversal arms
(`HF DR`, `HF STZ DR`) and a high-fat-without-diabetes arm; ST001412 carries a
lean arm. Putting a treated arm with the cases measures the treatment and
putting it with the controls measures protection — the trap MTBLS5667's
herbal-formula arm set. Any level named in neither `case` nor `control` is
dropped, so exclusion is the default and inclusion is explicit.
"""

from __future__ import annotations

import logging
import re

import pandas as pd

from .workbench_mwtab import Analysis

logger = logging.getLogger(__name__)

# A level spec is a list of exact normalised levels, or "re:<regex>".
LevelSpec = list[str] | str

# Factor keys as deposited differ in case between studies (DIET, diet, Group,
# Group1), so keys are matched case-insensitively, exactly as levels are.
WORKBENCH_STUDIES: dict[str, dict] = {
    "ST000603": {
        "case": {"Treatment": ["ic"]},
        "control": {"Treatment": ["control"]},
        "region": "urine", "condition": "interstitial cystitis",
    },
    # Case is `HF STZ`: high fat plus streptozotocin, the diabetic-neuropathy
    # model the study set out to phenotype. `HF` alone is diet-induced obesity
    # without diabetes — a different exposure, not a milder case — and the two
    # `DR` arms reverse the diet, so all three are excluded rather than
    # assigned. `Control` and `control` are one arm; level matching is
    # case-insensitive, so the deposited inconsistency merges on its own.
    "ST000676": {
        "case": {"DIET": ["hf stz"]},
        "control": {"DIET": ["control"]},
        "split_by": "SAMPLE_TYPE",
        "condition": "diabetic peripheral neuropathy",
    },
    # Genotype (wt, db+, ob+) and strain vary alongside diet. Both carriers are
    # heterozygous and phenotypically near-normal, so the diet contrast is
    # taken across them and the resulting strain/genotype spread is left as
    # unmodelled heterogeneity, which the random-effects pool is there to
    # absorb.
    "ST000780": {
        "case": {"diet": ["high fat"]},
        "control": {"diet": ["control"]},
        "region": "sciatic nerve", "condition": "high-fat-diet neuropathy",
    },
    # The matched control is `Obese non neuropathy`, not `Lean`. Comparing
    # obese-with-neuropathy against lean would confound neuropathy with
    # obesity and report the sum as a neuropathy effect; the lean arm is
    # therefore excluded.
    "ST001412": {
        "case": {"Group": ["obese neuropathy"]},
        "control": {"Group": ["obese non neuropathy"]},
        "region": "plasma", "condition": "obesity-associated neuropathy",
    },
    "ST003954": {
        "case": {"Group1": r"re:^ibs"},
        "control": {"Group1": ["hc"]},
        "region": "serum", "condition": "irritable bowel syndrome",
    },
    # Rheumatoid arthritis and juvenile idiopathic arthritis entered scope on
    # 2026-08-30 (conf/analysis/phenotype_scope.csv). All three deposits carry
    # several arms beyond the contrast, and in every case the extra arm is the
    # same trap MTBLS5667 and ST000676 document: folding a treatment arm into
    # the cases measures the drug, and folding it into the controls measures
    # protection. Level matching is exact after normalisation, so `ra` does not
    # capture `ra+mtx`, `ra-follow`, `pre-ra`, `rf-ra` or `ccp-ra`.
    #
    # `RA+MTX` is methotrexate-treated: the study's own subject is how the
    # plasma metabolome *normalises* under treatment, so that arm is the
    # endpoint of a different question.
    "ST001949": {
        "case": {"Condition": ["ra"]},
        "control": {"Condition": ["control"]},
        "region": "plasma", "condition": "rheumatoid arthritis",
    },
    # 210 RA subjects against 219 health. Four analysis rows per sample, which
    # run_workbench_study collapses; they are not repeated subjects. Five other
    # arms are excluded and each for its own reason: `RA-follow` is
    # longitudinal follow-up of 387 subjects and would enter the same people
    # twice; `pre-RA` is pre-clinical and not yet the disease; `OA` is a
    # different condition; `RF-RA` and `CCP-RA` are serology-defined subsets
    # whose subjects may overlap the RA arm, and the deposit gives no way to
    # tell; `-` is unlabelled.
    "ST003177": {
        "case": {"Clinical group": ["ra"]},
        "control": {"Clinical group": ["health"]},
        "region": "plasma", "condition": "rheumatoid arthritis",
    },
    # 112 JIA subjects against 72 healthy. The Crohn's disease arm is a third
    # condition, not a second control, and is excluded.
    "ST003520": {
        "case": {"Group": ["jia"]},
        "control": {"Group": ["healthy_control"]},
        "region": "plasma", "condition": "juvenile idiopathic arthritis",
    },
    # ST000218 is deliberately absent though its NORA arm is new-onset RA and
    # so within scope. It is a faecal short-chain fatty-acid panel of eight
    # metabolites; pooling faecal SCFAs with the plasma, serum and urine
    # metabolome of every other unit would compare compartments, not
    # conditions. Recorded in the screening file with that reason rather than
    # as a scope decision.
    # Stage and cycle phase vary within the cases only, so they cannot be
    # pinned without discarding most of the cohort; they stay unmodelled.
    "ST003984": {
        "case": {"Group": ["endometriosis"]},
        "control": {"Group": ["control"]},
        "split_by": "Sample source",
        "condition": "endometriosis",
    },
}


def norm_level(value: str) -> str:
    """Match `05j`'s normalisation, so the two agree on what a level is."""
    return re.sub(r"\s+", " ", (value or "").strip().lower())


def parse_factor_cell(cell: str) -> dict[str, str]:
    """`Key:Value | Key:Value` into a dict, keys and values normalised."""
    out: dict[str, str] = {}
    for part in (cell or "").split("|"):
        if ":" not in part:
            continue
        key, _, value = part.partition(":")
        key = norm_level(key)
        if key:
            out[key] = norm_level(value)
    return out


def level_matches(value: str, spec: LevelSpec) -> bool:
    """Whether one factor value satisfies a level spec."""
    if isinstance(spec, str):
        if not spec.startswith("re:"):
            raise ValueError(f"string level spec must start with 're:': {spec!r}")
        return re.search(spec[3:], value, re.I) is not None
    return value in {norm_level(s) for s in spec}


def in_arm(factors: dict[str, str], arm: dict[str, LevelSpec]) -> bool:
    """Whether a sample joins an arm: every declared factor must match.

    A sample missing a declared factor does not join. Treating an absent
    factor as a match would sweep unlabelled samples into whichever arm was
    tested first.
    """
    for key, spec in arm.items():
        value = factors.get(norm_level(key))
        if value is None or not level_matches(value, spec):
            return False
    return True


def split_values(analysis: Analysis, key: str) -> list[str]:
    """Distinct values of a splitting factor, in deposition order."""
    seen: dict[str, None] = {}
    for cell in analysis.factors:
        value = parse_factor_cell(cell).get(norm_level(key))
        if value:
            seen.setdefault(value, None)
    return list(seen)


def arm_columns(analysis: Analysis, spec: dict,
                split_value: str | None = None) -> tuple[list[str], list[str]]:
    """Sample ids belonging to the case and control arms.

    `split_value` restricts both arms to one level of `spec["split_by"]`, so a
    tissue-split study compares like with like.
    """
    split_by = spec.get("split_by")
    case, control = [], []
    for sample, cell in zip(analysis.samples, analysis.factors or
                            [""] * analysis.n_samples, strict=False):
        factors = parse_factor_cell(cell)
        if split_by and split_value is not None:
            if factors.get(norm_level(split_by)) != split_value:
                continue
        if in_arm(factors, spec["case"]):
            case.append(sample)
        elif in_arm(factors, spec["control"]):
            control.append(sample)
    return case, control


def analysis_frame(analysis: Analysis) -> tuple[pd.DataFrame, pd.Series]:
    """Abundances as a (metabolite x sample) frame, plus the metabolite names.

    Names are the study's own. Where the study supplies a RefMet name it is
    preferred, because RefMet is a standardised space and resolves to ChEBI far
    more reliably than a vendor's free text — the same reason the arm keys on
    accessions rather than spellings at all.
    """
    names, rows = [], []
    for raw_name, values in analysis.values.items():
        names.append(analysis.refmet.get(raw_name, raw_name))
        rows.append(values)
    frame = pd.DataFrame(rows, columns=analysis.samples)
    return frame, pd.Series(names, dtype=object)


def derived_study_id(accession: str, split_value: str) -> str:
    """`ST000676` + `dorsal root ganglia` -> `ST000676_Dorsal_root_ganglia`.

    Mirrors GSE241361_DRG: species and parentage are recovered by stripping
    the suffix, so the separator must not appear inside the accession.
    """
    slug = re.sub(r"[^A-Za-z0-9]+", "_", split_value.strip()).strip("_")
    return f"{accession}_{slug[:1].upper()}{slug[1:]}"

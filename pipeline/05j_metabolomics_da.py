#!/usr/bin/env python
"""Per-study differential abundance for the usable MetaboLights studies.

Emits the same schema as `05g_proteomics_pride_da.py`, so `06` consumes these
without special-casing: Hedges' g with its standard error, a Welch t-test
p-value, and BH adjustment within the study.

Three of the nine included studies survived file triage with populated
abundances, identified metabolites and a control arm
(`literature/prisma/metabolomics_file_triage.csv`). Their contrasts are declared
in `STUDIES` below rather than inferred, because choosing the comparator is a
design judgement that should be reviewable:

- **MTBLS13869** rat serum, recurrent pelvic pain. `animal model` against
  `control`, 5 v 5. Quality-control injections are apparatus and are dropped.
- **MTBLS2774** human serum, irritable bowel syndrome. `Irritable Bowel
  Syndrome` against `control`.
- **MTBLS5667** rat plasma, knee osteoarthritis. `osteoarthritis` against
  `normal control`, 6 v 6. **The `high-dose YQYXF` arm is excluded**, not
  folded into either side: it is a herbal-formula treatment arm, so putting it
  with the cases would measure the drug and putting it with the controls would
  measure protection. This is the same trap as PXD054342's prior-exercise arm.

**Ionisation modes are not separate studies.** Each of these deposits two to
four MAFs — LC-MS positive and negative, and for MTBLS2774 two of each — over
the *same* samples. Treating them as independent units is precisely the
pseudo-replication that let PXD013362 enter the proteomic pool fourteen times.
DA therefore runs per MAF, where the sample-to-group mapping is unambiguous,
and the study is then collapsed to one row per metabolite by smallest p-value,
the rule `05b`'s `aggregate_to_gene` and `05i` already apply.

**Features are keyed on ChEBI where possible.** `conf/analysis/metabolite_id_map.csv`
resolves names to ChEBI; without it the studies share no identifier namespace
and cannot pool at all. A metabolite with no mapping is emitted as
`NAME:<normalised>` so it is still visible and can still join another study
that used the same name — a weaker join than an accession, and marked as such
rather than silently mixed in with the ChEBI-keyed rows.

Linking MAF columns to samples takes three hops: MAF column -> assay row
(`MS Assay Name`, or a data-file basename, or the sample name itself) ->
`Sample Name` -> the sample table's `Factor Value`. Studies differ in which of
those the MAF column matches, so every alias is indexed and the column is
resolved against all of them.

Usage:
    python pipeline/05j_metabolomics_da.py
    python pipeline/05j_metabolomics_da.py --accessions MTBLS5667
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import logging
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from scipy import stats
from statsmodels.stats.multitest import multipletests

from cp_multiomics.metabolomics import workbench_da as wda
from cp_multiomics.metabolomics.workbench_mwtab import parse_mwtab

REPO_ROOT = Path(__file__).resolve().parent.parent
MAP = REPO_ROOT / "conf" / "analysis" / "metabolite_id_map.csv"
OUT_DIR = REPO_ROOT / "results" / "per_study" / "metabolomics"

FILES_API = "https://www.ebi.ac.uk/metabolights/ws/studies/{acc}/files"
FTP = "https://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/{acc}/{name}"
MWTAB = "https://www.metabolomicsworkbench.org/rest/study/study_id/{acc}/mwtab/txt"

logging.basicConfig(level=logging.INFO, format="%(message)s")
logger = logging.getLogger(__name__)

MIN_PER_GROUP = 3

SCHEMA = ["feature_id", "effect_size", "se", "pval", "control_mean",
          "treatment_mean", "n_control", "n_treatment", "study_id",
          "region", "condition", "cohort", "padj"]

MAF_FIXED = {
    "database_identifier", "chemical_formula", "smiles", "inchi",
    "metabolite_identification", "mass_to_charge", "fragmentation",
    "modifications", "charge", "retention_time", "taxid", "species",
    "database", "database_version", "reliability", "uri", "search_engine",
    "search_engine_score", "smallmolecule_abundance_sub",
    "smallmolecule_abundance_stdev_sub", "smallmolecule_abundance_std_error_sub",
    "chemical_shift", "multiplicity",
}

# Contrast per study. `case` and `control` are {Factor Value name: allowed
# levels}; a sample joins an arm only if it matches every entry. Levels are
# compared case-insensitively after whitespace collapse.
#
# The dict form supports factorial designs, where an arm is defined by several
# factors at once — needed to pin a treatment factor to its untreated level so
# the contrast measures disease rather than drug.
STUDIES: dict[str, dict] = {
    "MTBLS13869": {
        "case": {"Modeling status": ["animal model"]},
        "control": {"Modeling status": ["control"]},
        "region": "serum", "condition": "recurrent pelvic pain",
    },
    "MTBLS2774": {
        "case": {"Disease": ["irritable bowel syndrome"]},
        "control": {"Disease": ["control"]},
        "region": "serum", "condition": "irritable bowel syndrome",
    },
    "MTBLS5667": {
        "case": {"Cohort": ["osteoarthritis"]},
        "control": {"Cohort": ["normal control"]},
        "region": "plasma", "condition": "knee osteoarthritis",
    },
}

# Four included studies are deliberately absent, for two distinct reasons.
#
# MTBLS9662 has data and identifiers but its only factor is Cohort =
# Mild/Moderate/Severe, 30 each: a severity gradient among cases with no
# unaffected group. A severity trend is a different estimand from a
# case-versus-control effect size, and pooling the two would compare severe
# disease against mild disease and report it as a disease effect.
#
# MTBLS13513, MTBLS14059 and MTBLS14295 **deposit no abundances**. Each
# declares sample columns beside named metabolites — 91 columns and 1,225
# metabolites for MTBLS13513 — and every cell is empty. Their group labels are
# perfectly good and were extracted successfully (50 v 30, 12 v 12, 4 v 4);
# there is simply nothing to compute an effect size from. The file triage
# originally called them usable because it counted declared columns rather than
# populated ones; `scripts/triage_metabolomics.py` now counts values.

# Quality-control and blank injections are apparatus, not biology, and must
# never land in an arm. Two mechanisms are needed because the metadata cannot
# be trusted alone: MTBLS13513, MTBLS13869 and MTBLS5667 annotate their QC
# samples properly, but MTBLS14059 and MTBLS14295 label *every* sample
# "experimental blank", so keying on that string would exclude the entire
# study. Their QC injections carry the control arm's factor levels — seven of
# MTBLS14059's thirteen (0,0,0) samples are QC — so without the name fallback
# they would silently pad the control group.
QC_SAMPLE_TYPES = {
    "pooled quality control sample", "sample preparation blank",
    "solvent blank", "quality control",
}
QC_NAME = re.compile(r"^(qc|blank|pool(ed)?)[\W_]*\d*$", re.I)

UNNAMED = re.compile(
    r"^(cluster[_\- ]?\d+|unknown.*|unidentified.*|peak[_\- ]?\d+|m\d+t\d+"
    r"|feature[_\- ]?\d+|metabolite[_\- ]?\d+|compound[_\- ]?\d+|\d+(\.\d+)?)$",
    re.I,
)


def norm_name(name: str) -> str:
    """Match `build_metabolite_id_map.py`: keeps l-/d- stereo prefixes."""
    return re.sub(r"[^a-z0-9]+", "", name.strip().lower())


def norm_level(value: str) -> str:
    return re.sub(r"\s+", " ", value.strip().lower())


def fetch(url: str) -> bytes:
    r = requests.get(url, timeout=180)
    r.raise_for_status()
    return r.content


def read_tsv(acc: str, name: str) -> tuple[list[str], list[list[str]]]:
    text = fetch(FTP.format(acc=acc, name=name)).decode("utf-8", "replace")
    rows = list(csv.reader(io.StringIO(text), delimiter="\t"))
    if not rows:
        return [], []
    return [c.strip() for c in rows[0]], rows[1:]


def list_files(acc: str) -> list[str]:
    data = json.loads(fetch(FILES_API.format(acc=acc)))
    return [f["file"] for f in data.get("study", []) or [] if not f.get("directory")]


def sample_groups(acc: str, sfile: str, spec: dict) -> tuple[dict[str, str], int]:
    """Map Sample Name -> 'case' | 'control'; return also the QC count dropped.

    A sample joins an arm only if it matches every factor the arm names, which
    is what keeps the factorial studies' treatment arms out of the contrast.
    """
    header, rows = read_tsv(acc, sfile)
    if "Sample Name" not in header:
        raise SystemExit(f"{acc}: no 'Sample Name' in {sfile}")
    si = header.index("Sample Name")
    ti = header.index("Characteristics[Sample type]") \
        if "Characteristics[Sample type]" in header else -1

    wanted = {f for arm in ("case", "control") for f in spec[arm]}
    idx = {}
    for f in wanted:
        col = f"Factor Value[{f}]"
        if col not in header:
            raise SystemExit(f"{acc}: no {col!r} in {sfile}")
        idx[f] = header.index(col)

    def matches(row: list[str], arm: dict[str, list[str]]) -> bool:
        for f, levels in arm.items():
            i = idx[f]
            if len(row) <= i or norm_level(row[i]) not in {norm_level(x) for x in levels}:
                return False
        return True

    groups: dict[str, str] = {}
    dropped_qc = 0
    for r in rows:
        if len(r) <= si or not r[si].strip():
            continue
        sample = r[si].strip()
        stype = norm_level(r[ti]) if ti >= 0 and len(r) > ti else ""
        if stype in QC_SAMPLE_TYPES or QC_NAME.match(sample):
            dropped_qc += 1
            continue
        if matches(r, spec["case"]):
            groups[sample] = "case"
        elif matches(r, spec["control"]):
            groups[sample] = "control"
    return groups, dropped_qc


def column_aliases(acc: str, assay_files: list[str]) -> dict[str, str]:
    """Every identifier an assay row can be known by -> its Sample Name."""
    alias: dict[str, str] = {}
    for af in assay_files:
        header, rows = read_tsv(acc, af)
        if "Sample Name" not in header:
            continue
        si = header.index("Sample Name")
        idx = [header.index(c) for c in
               ("MS Assay Name", "NMR Assay Name", "Raw Spectral Data File",
                "Derived Spectral Data File", "Free Induction Decay Data File")
               if c in header]
        for r in rows:
            if len(r) <= si or not r[si].strip():
                continue
            sample = r[si].strip()
            alias.setdefault(sample.lower(), sample)
            for i in idx:
                if len(r) > i and r[i].strip():
                    v = r[i].strip()
                    alias.setdefault(v.lower(), sample)
                    base = Path(v).name
                    alias.setdefault(base.lower(), sample)
                    alias.setdefault(re.sub(r"\.(raw|mzxml|mzml|mzdata|xml|zip|d)$", "",
                                            base, flags=re.I).lower(), sample)
    return alias


def resolve_column(col: str, alias: dict[str, str]) -> str | None:
    """MAF columns carry suffixes the assay files do not (e.g. SAMP001_AUC)."""
    candidates = [col, re.sub(r"_(auc|area|height|intensity|abundance)$", "", col, flags=re.I)]
    for c in candidates:
        hit = alias.get(c.lower())
        if hit:
            return hit
    return None


def hedges_g(case: np.ndarray, ctrl: np.ndarray) -> tuple[float, float]:
    n1, n2 = len(case), len(ctrl)
    df = n1 + n2 - 2
    if df <= 0:
        return np.nan, np.nan
    v1, v2 = case.var(ddof=1), ctrl.var(ddof=1)
    sp2 = ((n1 - 1) * v1 + (n2 - 1) * v2) / df
    if not np.isfinite(sp2) or sp2 <= 0:
        return np.nan, np.nan
    d = (case.mean() - ctrl.mean()) / np.sqrt(sp2)
    j = 1.0 - 3.0 / (4.0 * df - 1.0)
    g = j * d
    var_g = (j ** 2) * ((n1 + n2) / (n1 * n2) + d ** 2 / (2.0 * df))
    return g, float(np.sqrt(var_g))


def maf_frame(acc: str, mfile: str) -> tuple[pd.DataFrame, pd.Series]:
    """Return (values indexed by row, metabolite names) for one MAF."""
    header, rows = read_tsv(acc, mfile)
    if "metabolite_identification" not in header:
        return pd.DataFrame(), pd.Series(dtype=str)
    mi = header.index("metabolite_identification")
    scols = [(i, c) for i, c in enumerate(header) if c and c not in MAF_FIXED]
    names, data = [], []
    for r in rows:
        if len(r) <= mi:
            continue
        raw = r[mi].strip()
        if not raw or UNNAMED.match(raw):
            continue
        names.append(raw)
        data.append([r[i].strip() if len(r) > i else "" for i, _ in scols])
    frame = pd.DataFrame(data, columns=[c for _, c in scols])
    frame = frame.apply(pd.to_numeric, errors="coerce")
    return frame, pd.Series(names)


def differential(values: pd.DataFrame, names: pd.Series, case_cols: list[str],
                 ctrl_cols: list[str], assay: str) -> pd.DataFrame:
    """Hedges' g per metabolite for one MAF."""
    # Intensities are concentrations on a raw scale; log2 stabilises variance
    # and zeros are absent measurements rather than true zeros.
    vals = values.where(values > 0)
    if np.nanmax(vals.to_numpy(dtype=float)) > 100:
        vals = np.log2(vals)
    case_mat = vals[case_cols].to_numpy(dtype=float)
    ctrl_mat = vals[ctrl_cols].to_numpy(dtype=float)
    out = []
    for i in range(len(vals)):
        case = case_mat[i][np.isfinite(case_mat[i])]
        ctrl = ctrl_mat[i][np.isfinite(ctrl_mat[i])]
        if len(case) < MIN_PER_GROUP or len(ctrl) < MIN_PER_GROUP:
            continue
        g, se = hedges_g(case, ctrl)
        if not np.isfinite(g) or not np.isfinite(se) or se <= 0:
            continue
        out.append({
            "name": names.iloc[i], "effect_size": g, "se": se,
            "pval": float(stats.ttest_ind(case, ctrl, equal_var=False).pvalue),
            "control_mean": ctrl.mean(), "treatment_mean": case.mean(),
            "n_control": len(ctrl), "n_treatment": len(case), "assay": assay,
        })
    return pd.DataFrame(out)


def load_map() -> dict[str, str]:
    if not MAP.exists():
        raise SystemExit(f"no metabolite map at {MAP}; run "
                         "scripts/build_metabolite_id_map.py --apply first")
    m = pd.read_csv(MAP, dtype=str).fillna("")
    m = m[m["chebi_id"] != ""]
    return dict(zip(m["normalised_name"], m["chebi_id"]))


def run_study(acc: str, spec: dict, id_map: dict[str, str]) -> pd.DataFrame:
    names_ = list_files(acc)
    mafs = sorted(n for n in names_ if n.startswith("m_") and n.endswith(".tsv"))
    sfiles = sorted(n for n in names_ if n.startswith("s_") and n.endswith(".txt"))
    assays = sorted(n for n in names_ if n.startswith("a_") and n.endswith(".txt"))
    if not mafs or not sfiles:
        raise SystemExit(f"{acc}: missing MAF or sample table")

    groups, dropped_qc = sample_groups(acc, sfiles[0], spec)
    alias = column_aliases(acc, assays)
    logger.info("  %d samples in the contrast (%d case, %d control); "
                "%d QC/blank dropped", len(groups),
                sum(v == "case" for v in groups.values()),
                sum(v == "control" for v in groups.values()), dropped_qc)

    parts = []
    for mfile in mafs:
        values, names = maf_frame(acc, mfile)
        if values.empty:
            continue
        # One column per sample. MTBLS13869's MAF carries both `SAMP001_AUC`
        # and `SAMP001`, which resolve to the same animal; taking both would
        # count every animal twice, halving the standard error. Here the plain
        # columns happen to be empty so the finite-value filter would have
        # hidden it, which is exactly why the rule is enforced rather than
        # relied upon. The column with the most measurements wins.
        by_sample: dict[str, list[str]] = {}
        unmapped = 0
        for col in values.columns:
            sample = resolve_column(col, alias)
            if sample is None:
                unmapped += 1
                continue
            by_sample.setdefault(sample, []).append(col)

        case_cols, ctrl_cols, collisions = [], [], 0
        for sample, cols in by_sample.items():
            g = groups.get(sample)
            if g not in ("case", "control"):
                continue
            if len(cols) > 1:
                collisions += 1
                cols = sorted(cols, key=lambda c: (-int(values[c].notna().sum()), c))
            (case_cols if g == "case" else ctrl_cols).append(cols[0])
        assay = mfile.replace(f"m_{acc}_", "").replace("_v2_maf.tsv", "")
        logger.info("    %-50s case=%d ctrl=%d samples (unresolved cols %d, "
                    "collapsed duplicates %d)",
                    assay[:50], len(case_cols), len(ctrl_cols), unmapped, collisions)
        if len(case_cols) < MIN_PER_GROUP or len(ctrl_cols) < MIN_PER_GROUP:
            continue
        part = differential(values, names, case_cols, ctrl_cols, assay)
        if not part.empty:
            parts.append(part)

    if not parts:
        return pd.DataFrame(columns=SCHEMA)
    df = pd.concat(parts, ignore_index=True)

    # One row per metabolite per study. Ionisation modes measure the same
    # samples, so keeping both would be within-study replication; the most
    # significant represents the metabolite.
    return finalise(df, acc, spec["region"], spec["condition"], id_map)



def finalise(df: pd.DataFrame, study_id: str, region: str, condition: str,
             id_map: dict[str, str]) -> pd.DataFrame:
    """Shared tail of both repositories' per-study analysis.

    Collapse by name, key on ChEBI, collapse again on the identifier, then
    adjust within the study. Factored out of `run_study` so the Workbench path
    cannot drift from the MetaboLights one — a copy of this would be a copy
    that stops matching.

    The second collapse is the sample-independence guard. The map deliberately
    merges synonyms — Erucamide and 13Z-Docosenamide are both CHEBI:142245 —
    so two distinct names in one study can become one identifier, and the
    by-name collapse cannot see it. Left in, MTBLS13869 contributed
    CHEBI:142245 twice and the feature reached the pool at k=3 with only two
    real studies behind it: a sample-independence violation arriving through
    the identifier map rather than through the assay files.
    """
    df["normalised_name"] = df["name"].map(norm_name)
    df = df.sort_values("pval").drop_duplicates(subset=["normalised_name"], keep="first")

    df["feature_id"] = df["normalised_name"].map(id_map)
    unmapped = df["feature_id"].isna()
    df.loc[unmapped, "feature_id"] = "NAME:" + df.loc[unmapped, "normalised_name"]

    before = len(df)
    df = df.sort_values("pval").drop_duplicates(subset=["feature_id"], keep="first")
    merged = before - len(df)

    n_chebi = int(df["feature_id"].str.startswith("CHEBI").sum())
    logger.info("  %d metabolites, %d keyed on ChEBI, %d on name only%s",
                len(df), n_chebi, len(df) - n_chebi,
                f" ({merged} synonym rows merged)" if merged else "")

    df["study_id"] = study_id
    df["region"] = region
    df["condition"] = condition
    df["cohort"] = ""
    df["padj"] = multipletests(df["pval"], method="fdr_bh")[1]
    return df[SCHEMA + ["assay"]].sort_values("pval")


def run_workbench_study(acc: str, spec: dict,
                        id_map: dict[str, str]) -> dict[str, pd.DataFrame]:
    """Per-study differential abundance for one Workbench accession.

    Returns one frame per emitted study id: a single entry keyed on the
    accession, or one per tissue where the contrast declares `split_by`.

    An mwtab file can hold several analyses over the *same* samples — ST003984
    deposits three. They are analysed separately and concatenated here, then
    collapsed to one row per metabolite by `finalise`, exactly as the
    MetaboLights path collapses its several MAFs. Treating them as independent
    units is the pseudo-replication that let PXD013362 into the proteomic pool
    fourteen times.
    """
    text = fetch(MWTAB.format(acc=acc)).decode("utf-8", "replace")
    analyses = parse_mwtab(text)
    if not analyses:
        logger.warning("  no abundance matrix in mwtab")
        return {}
    logger.info("  %d analysis document(s) over %s samples",
                len(analyses), "/".join(str(a.n_samples) for a in analyses))

    split_by = spec.get("split_by")
    splits: list[tuple[str, str | None]]
    if split_by:
        values: list[str] = []
        for a in analyses:
            for v in wda.split_values(a, split_by):
                if v not in values:
                    values.append(v)
        splits = [(wda.derived_study_id(acc, v), v) for v in values]
    else:
        splits = [(acc, None)]

    out: dict[str, pd.DataFrame] = {}
    for study_id, split_value in splits:
        frames = []
        for a in analyses:
            case_cols, ctrl_cols = wda.arm_columns(a, spec, split_value)
            if len(case_cols) < MIN_PER_GROUP or len(ctrl_cols) < MIN_PER_GROUP:
                continue
            values, names = wda.analysis_frame(a)
            frames.append(differential(values, names, case_cols, ctrl_cols,
                                       a.analysis_id))
            logger.info("    %-12s %-26s case=%d ctrl=%d",
                        a.analysis_id, split_value or "", len(case_cols), len(ctrl_cols))
        frames = [f for f in frames if not f.empty]
        if not frames:
            logger.warning("    %s: no arm reached %d per group", study_id, MIN_PER_GROUP)
            continue
        region = split_value or spec.get("region", "")
        out[study_id] = finalise(pd.concat(frames, ignore_index=True), study_id,
                                 region, spec["condition"], id_map)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--accessions", help="comma-separated subset")
    args = ap.parse_args()

    accs = ([a.strip() for a in args.accessions.split(",")]
            if args.accessions else list(STUDIES) + list(wda.WORKBENCH_STUDIES))
    id_map = load_map()
    logger.info("metabolite map: %d names -> ChEBI\n", len(id_map))

    total, written = 0, 0
    for acc in accs:
        if acc in STUDIES:
            spec = STUDIES[acc]
            logger.info("=== %s (%s) ===", acc, spec["condition"])
            emitted = {acc: run_study(acc, spec, id_map)}
        elif acc in wda.WORKBENCH_STUDIES:
            spec = wda.WORKBENCH_STUDIES[acc]
            logger.info("=== %s (%s) [Workbench] ===", acc, spec["condition"])
            emitted = run_workbench_study(acc, spec, id_map)
        else:
            raise SystemExit(f"{acc}: no contrast declared for this accession")

        for study_id, df in emitted.items():
            if df.empty:
                logger.warning("  %s: no effect sizes; nothing written", study_id)
                continue
            out = OUT_DIR / study_id / f"{study_id}_effects.csv"
            out.parent.mkdir(parents=True, exist_ok=True)
            df.to_csv(out, index=False)
            sig = int((df["padj"] < 0.05).sum())
            logger.info("  -> %s  %d features, %d at padj<0.05\n",
                        out.relative_to(REPO_ROOT), len(df), sig)
            total += len(df)
            written += 1

    logger.info("%d effect sizes across %d study units (%d accessions)",
                total, written, len(accs))
    return 0


if __name__ == "__main__":
    sys.exit(main())

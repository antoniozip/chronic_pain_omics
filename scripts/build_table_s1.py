"""Regenerate Supplementary Table S1 (transcriptomic studies).

The table shipped with five columns while the manuscript advertised six
fields, three of which it did not contain. This script derives every column
from a named source so the table cannot drift from the analysis again:

    Study_ID, Pain_Model      pipeline/06b_stratified_meta.R (STUDY_MODEL),
                              merged with conf/analysis/study_strata.csv the
                              way 06b merges them - the file wins, and a
                              derived unit inherits from its parent accession
    Species, Platform, Year   data/interim/transcriptomics/harmonized.jsonl
    Tissue                    data/interim/sni_tissue_annotations.csv
    N_Case, N_Control, N_Analysed
                              results/per_study/transcriptomics/*_effects.csv
    N_Features, Mean_Abs_log2FC
                              results/per_study/transcriptomics/
                              effects_normalized_all.csv - the features that
                              actually enter the meta-analysis, which is fewer
                              than the raw per-study tables carry

Tissue is annotated only for the SNI studies that step 06d stratifies, so it
is reported as "not annotated" elsewhere rather than guessed.

    uv run python scripts/build_table_s1.py
"""

from __future__ import annotations

import argparse
import json
import logging
import re
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
logger = logging.getLogger(__name__)

REPO_ROOT = Path(__file__).resolve().parents[1]
PER_STUDY_DIR = REPO_ROOT / "results" / "per_study" / "transcriptomics"
MANIFEST = REPO_ROOT / "data" / "interim" / "transcriptomics" / "harmonized.jsonl"
STRATIFIED_R = REPO_ROOT / "pipeline" / "06b_stratified_meta.R"
STRATA_CSV = REPO_ROOT / "conf" / "analysis" / "study_strata.csv"
TISSUE_CSV = REPO_ROOT / "data" / "interim" / "sni_tissue_annotations.csv"
SUPERSEDED_CSV = REPO_ROOT / "conf" / "analysis" / "superseded_studies.csv"
DEFAULT_OUT = (REPO_ROOT / "manuscript" / "supplementary"
               / "Table_S1_transcriptomic_studies.csv")

UNKNOWN = "not annotated"
_STUDY_MODEL_RE = re.compile(r"(GSE[A-Za-z0-9_]+)\s*=\s*\"([^\"]+)\"")
_COMMON_NAME = {
    "Mus musculus": "Mouse",
    "Rattus norvegicus": "Rat",
    "Homo sapiens": "Human",
}


def load_pain_models(r_path: Path, strata_path: Path = STRATA_CSV) -> dict[str, str]:
    """The stratum of every study, exactly as 06b resolves it.

    06b reads two sources: the STUDY_MODEL vector in its own source, and
    conf/analysis/study_strata.csv, which the 2026-08-28 amendment introduced
    so that a retrieval no longer edits analysis source. The file wins where
    both name a study. Reading only the vector left the 23 amendment studies
    in this table reported as "not annotated" while 06b was stratifying them.
    """
    text = r_path.read_text()
    start = text.index("STUDY_MODEL <- c(")
    block = text[start:text.index("\n)", start)]
    models = dict(_STUDY_MODEL_RE.findall(block))
    if not models:
        raise ValueError(f"no STUDY_MODEL entries parsed from {r_path}")
    if strata_path.exists():
        extra = pd.read_csv(strata_path).dropna(subset=["study_id", "stratum"])
        models.update(zip(extra["study_id"].astype(str),
                          extra["stratum"].astype(str)))
    return models


def resolve_pain_model(study_id: str, models: dict[str, str]) -> str:
    """Look a study up, falling back to its parent accession.

    A derived unit (GSE180627_S1, GSE241361_DRG) is in no accession-keyed map
    by construction, so without this it reports as unstratified while 06b
    pools it. Mirrors resolve_by_parent() in R/meta_utils.R; an exact hit wins,
    so a split that changes the stratum can still be recorded per unit.
    """
    candidate = study_id
    while True:
        if candidate in models:
            return models[candidate]
        parent = candidate.rsplit("_", 1)[0]
        if parent == candidate:
            return UNKNOWN
        candidate = parent


def load_manifest(path: Path) -> dict[str, dict]:
    """Index the harmonized manifest by accession."""
    records: dict[str, dict] = {}
    with open(path) as f:
        for line in f:
            if line.strip():
                rec = json.loads(line)
                records[rec["accession"]] = rec
    return records


def load_tissues(path: Path) -> dict[str, str]:
    if not path.exists():
        logger.warning("no tissue annotations at %s", path)
        return {}
    df = pd.read_csv(path)
    return dict(zip(df["study_id"], df["tissue_category"]))


def _species(rec: dict | None) -> str:
    if not rec:
        return UNKNOWN
    names = rec.get("species_canonical") or []
    if isinstance(names, str):
        names = [names]
    return "; ".join(_COMMON_NAME.get(n, n) for n in names) or UNKNOWN


def _platform(rec: dict | None) -> str:
    """Manifest stores a bare GEO platform id; restore the GPL prefix."""
    if not rec:
        return UNKNOWN
    raw = str(rec.get("platform_canonical") or "").strip()
    if not raw:
        return UNKNOWN
    # Multi-platform studies are stored as "87;86;85"; every id needs the
    # prefix, not just the first.
    ids = [p.strip() for p in raw.split(";") if p.strip()]
    return ";".join(p if p.upper().startswith("GPL") else f"GPL{p}" for p in ids)


def _first_int(df: pd.DataFrame, column: str) -> int | None:
    """First value of `column` as an int, or None when absent or missing.

    Some per-study tables carry the column with NaN values, so presence of the
    column is not sufficient.
    """
    if column not in df.columns or df.empty:
        return None
    value = df[column].iloc[0]
    return None if pd.isna(value) else int(value)


def build(per_study_dir: Path, out_path: Path) -> pd.DataFrame:
    pain_models = load_pain_models(STRATIFIED_R)
    # Feature counts and effect magnitudes must describe the analysed set, not
    # the raw per-study tables: ID normalization drops features that never
    # reach the meta-analysis.
    normalized = pd.read_csv(per_study_dir / "effects_normalized_all.csv",
                             usecols=["study_id", "effect_size"], low_memory=False)
    analysed = normalized.groupby("study_id")["effect_size"].agg(
        N_Features="size", Mean_Abs_log2FC=lambda s: s.abs().mean())
    manifest = load_manifest(MANIFEST)
    tissues = load_tissues(TISSUE_CSV)
    # The table must list the studies that entered the analysis, so it honours
    # the same exclusion registry as steps 05b and 06.
    superseded = (set(pd.read_csv(SUPERSEDED_CSV)["study_id"].astype(str))
                  if SUPERSEDED_CSV.exists() else set())

    rows: list[dict] = []
    for effects_path in sorted(per_study_dir.glob("*_effects.csv")):
        study_id = effects_path.name.removesuffix("_effects.csv")
        if study_id in superseded:
            logger.info("excluding %s: superseded (%s)", study_id, SUPERSEDED_CSV.name)
            continue
        df = pd.read_csv(effects_path)
        # Split studies (GSE241361_DRG) carry a suffix the manifest lacks.
        rec = manifest.get(study_id) or manifest.get(study_id.split("_")[0])

        n_case = _first_int(df, "n_case")
        n_control = _first_int(df, "n_control")

        rows.append({
            "Study_ID": study_id,
            "Species": _species(rec),
            "Pain_Model": resolve_pain_model(study_id, pain_models),
            "Tissue": tissues.get(study_id, UNKNOWN),
            "Platform": _platform(rec),
            "Year": (rec or {}).get("year", UNKNOWN),
            "N_Case": n_case,
            "N_Control": n_control,
            "N_Analysed": None if n_case is None or n_control is None else n_case + n_control,
            "N_Features": (int(analysed.loc[study_id, "N_Features"])
                           if study_id in analysed.index else None),
            "Mean_Abs_log2FC": (round(float(analysed.loc[study_id, "Mean_Abs_log2FC"]), 4)
                                if study_id in analysed.index else None),
        })

    table = pd.DataFrame(rows).sort_values("Study_ID").reset_index(drop=True)
    # Nullable integers: several studies have no recorded sample size, and a
    # plain int column would render the rest as floats ("21.0").
    for col in ("N_Case", "N_Control", "N_Analysed", "N_Features"):
        table[col] = table[col].astype("Int64")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(out_path, index=False)

    missing_model = int((table["Pain_Model"] == UNKNOWN).sum())
    missing_tissue = int((table["Tissue"] == UNKNOWN).sum())
    missing_n = int(table["N_Analysed"].isna().sum())
    logger.info("%d studies -> %s", len(table), out_path)
    logger.info("unannotated: pain model %d, tissue %d, sample size %d",
                missing_model, missing_tissue, missing_n)
    return table


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--per-study-dir", type=Path, default=PER_STUDY_DIR)
    p.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = p.parse_args()
    build(args.per_study_dir, args.out)


if __name__ == "__main__":
    main()

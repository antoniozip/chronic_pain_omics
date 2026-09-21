"""Pseudobulk each single-cell study to one counts column per sample.

Summing across the cells within a sample is what makes this arm comparable
with the bulk arm. A per-sample pseudobulk profile *is* a bulk profile, so the
existing per-study DA machinery applies unchanged and the sample rather than
the cell is the unit of replication. Testing across cells instead treats them
as replicates, inflates p-values by orders of magnitude, and would make every
number the arm produces incomparable with the pooled bulk estimates.

The ten included studies deposit their matrices in four shapes and the
`LAYOUT` table below says which each uses, where its files are, and how a GSM
maps onto them. `src/cp_multiomics/single_cell/readers.py` holds the readers.

Two things this step is the authority on, and step 18 cannot be:

  * **Which samples actually carry data.** A GSM can be described in the
    series metadata and absent from the deposit. GSE155622 describes 18 10x
    libraries and deposits 10. Those are recorded as `absent_from_deposit`,
    not dropped, and the study's arm counts are re-checked afterwards.
  * **How many cells each sample contributed.** Recorded per GSM in the
    manifest, because a pseudobulk profile summed over 40 cells and one summed
    over 12,000 are not equally trustworthy and nothing downstream can tell
    them apart.

Where step 18 assigned two libraries of one animal the same subject id, their
counts are summed here into one column, so the animal enters once.

Outputs:
    data/interim/single_cell/pseudobulk/{accession}_counts.tsv.gz
    data/interim/single_cell/pseudobulk_manifest.csv

Usage:
    python pipeline/19_pseudobulk_single_cell.py
    python pipeline/19_pseudobulk_single_cell.py --accession GSE328175
    python pipeline/19_pseudobulk_single_cell.py --min-umi-per-cell 200
"""

from __future__ import annotations

import argparse
import logging
import re
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "src"))

from cp_multiomics.single_cell import (  # noqa: E402
    Pseudobulk,
    read_10x_h5,
    read_dense_split,
    read_lines,
    read_mtx,
    strip_version,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                    datefmt="%H:%M:%S")
logger = logging.getLogger(__name__)

RAW = REPO_ROOT / "data" / "raw" / "single_cell"
INTERIM = REPO_ROOT / "data" / "interim" / "single_cell"
GROUPS_DIR = INTERIM / "groups"
OUT_DIR = INTERIM / "pseudobulk"
MANIFEST = INTERIM / "pseudobulk_manifest.csv"
SHARD_DIR = INTERIM / "pseudobulk_manifest.d"

#: Barcodes below this many UMIs are not counted as cells. See the readers
#: module: this is a cell call for the deposits that ship a raw droplet matrix,
#: and a no-op for the ones that ship CellRanger's filtered output.
DEFAULT_MIN_UMI = 500

#: GSE155622's merged tables name each cell by its library and an index. The
#: map is written out rather than derived because the deposit's library names
#: ("SNI 6h_2") and the series titles ("SNI 6h_rep2") do not share a form.
#:
#: The library names themselves are spelled two ways *within the deposit*: the
#: metadata files carry "SNI 6h" and the counts tables carry "SNI.6h", because
#: R's make.names() replaced the space when the matrix columns were written.
#: Both spellings are keys. Keying only the metadata's spelling mapped the
#: three control libraries and none of the seven SNI ones, which left the study
#: with an empty case arm -- a shape indistinguishable from a deposit that
#: simply omits its case samples, and eight of GSE155622's samples really are
#: omitted, so the warning it produced looked entirely ordinary.
GSE155622_PREFIXES = {
    prefix: gsm
    for name, gsm in {
        "Ctrl_1": "GSM4708669", "Ctrl_2": "GSM4708672", "Ctrl_3": "GSM4708676",
        "SNI 6h": "GSM4708670", "SNI 6h_2": "GSM4708677",
        "SNI 24h": "GSM4708671", "SNI 24h_2": "GSM4708678",
        "SNI 2d": "GSM4708673", "SNI 7d": "GSM4708674", "SNI 14d": "GSM4708675",
    }.items()
    for prefix in {name, name.replace(" ", ".")}
}

#: Per-study deposit layout. `kind` selects the reader; `id_space` says what
#: the emitted feature ids are, which step 05b needs in order to leave the
#: already-symbol studies alone.
LAYOUT: dict[str, dict] = {
    "GSE134003": {
        "kind": "dense_merged", "id_space": "symbol",
        "file": "GSE134003_SI_SNI_rawcounts.txt.gz",
        "prefix": r"^([A-Z]+\d+)_",       # SI01_AAAATCGCACCG -> SI01
        "prefix_is": "title",
    },
    "GSE155622": {
        "kind": "dense_merged", "id_space": "symbol",
        "file": ["GSE155622_raw_UMI_counts_1.txt.gz",
                 "GSE155622_raw_UMI_counts_2.txt.gz"],
        "prefix": r"^(.*)_\d+$",          # Ctrl_1_24 -> Ctrl_1
        "prefix_is": "map", "map": GSE155622_PREFIXES,
    },
    "GSE162807": {"kind": "h5", "id_space": "symbol",
                  "glob": "unpacked/{gsm}_*_bc_matri*.h5"},
    "GSE179640": {"kind": "h5", "id_space": "symbol",
                  "glob": "unpacked/{gsm}_*_bc_matrix.h5"},
    "GSE198608": {
        "kind": "dense_per_gsm", "id_space": "symbol",
        "glob": "unpacked/{gsm}_*.umi_counts.xls.gz",
        # Row ids are R-mangled Gene_id_unique; the deposit ships the map.
        "feature_map": ("GSE198608_Rnor_6.0.102.csv.gz", "Gene_id_unique", "Gene_id_raw"),
    },
    "GSE214411": {"kind": "mtx_10x", "id_space": "symbol",
                  "glob": "unpacked/{gsm}_*_matrix.mtx.gz"},
    "GSE254360": {"kind": "mtx_10x", "id_space": "symbol",
                  "glob": "unpacked/{gsm}_*_matrix.mtx.gz"},
    "GSE283177": {
        "kind": "mtx_sidecar", "id_space": "ensembl_gene", "cells_are_rows": True,
        "glob": "unpacked/{gsm}_*_cells_x_genes.mtx.gz",
        "features": "_cells_x_genes.genes.txt.gz", "feature_column": 0,
        "strip_version": True,
    },
    "GSE289659": {
        "kind": "mtx_sidecar", "id_space": "symbol", "cells_are_rows": True,
        "glob": "unpacked/{gsm}_*_DGE_filtered_DGE.mtx.gz",
        "features": "_DGE_filtered_all_genes.csv.gz", "feature_column": 1,
        "feature_sep": ",", "feature_header": True,
    },
    "GSE328175": {"kind": "mtx_10x", "id_space": "symbol",
                  "glob": "unpacked/{gsm}_*_matrix.mtx.gz"},
}


def one_file(directory: Path, pattern: str) -> Path | None:
    hits = sorted(directory.glob(pattern))
    return hits[0] if hits else None


def read_gsm(accession: str, gsm: str, spec: dict, min_umi: int) -> Pseudobulk | None:
    """Pseudobulk one GSM of a per-sample deposit."""
    root = RAW / accession
    mtx = one_file(root, spec["glob"].format(gsm=gsm))
    if mtx is None:
        return None

    if spec["kind"] == "h5":
        return read_10x_h5(mtx, gsm, min_umi)

    if spec["kind"] == "dense_per_gsm":
        return read_dense_split(mtx, lambda _c: gsm, min_umi)

    if spec["kind"] == "mtx_10x":
        features_path = Path(str(mtx).replace("matrix.mtx.gz", "features.tsv.gz"))
        if not features_path.exists():
            features_path = Path(str(mtx).replace("matrix.mtx.gz", "genes.tsv.gz"))
        # column 1 is the gene symbol in a 10x features file, column 0 the id
        features = read_lines(features_path, column=1)
        return read_mtx(mtx, features, gsm, min_umi, cells_are_rows=False)

    features_path = one_file(root / "unpacked", f"{gsm}_*{spec['features']}")
    if features_path is None:
        return None
    features = read_lines(features_path, column=spec.get("feature_column", 0),
                          sep=spec.get("feature_sep", "\t"),
                          skip_header=spec.get("feature_header", False))
    if spec.get("strip_version"):
        features = strip_version(features)
    return read_mtx(mtx, features, gsm, min_umi,
                    cells_are_rows=spec.get("cells_are_rows", False))


def read_merged(accession: str, spec: dict, groups: pd.DataFrame,
                min_umi: int) -> Pseudobulk:
    """Pseudobulk a study deposited as one table carrying several samples."""
    prefix_re = re.compile(spec["prefix"])
    if spec["prefix_is"] == "title":
        lookup = dict(zip(groups["title"].astype(str), groups["sample_id"]))
    else:
        lookup = spec["map"]

    unmapped: Counter[str] = Counter()

    def owner(column: str) -> str | None:
        found = prefix_re.match(column.strip('"'))
        prefix = found.group(1) if found else column.strip('"')
        sample = lookup.get(prefix) if found else None
        if sample is None:
            unmapped[prefix] += 1
        return sample

    # A study deposited as several tables is several submissions, and two
    # submissions agree on their feature order only by luck: GSE155622's two
    # tables do not, any more than GSE198608's two batches do. They go through
    # the same aligner as the per-sample deposits.
    files = spec["file"] if isinstance(spec["file"], list) else [spec["file"]]
    parts = [read_dense_split(RAW / accession / name, owner, min_umi)
             for name in files]

    assert_every_prefix_maps(accession, unmapped, spec.get("unmapped_prefixes", ()))
    return align_features(accession, parts) if len(parts) > 1 else parts[0]


def assert_every_prefix_maps(accession: str, unmapped: Counter[str],
                             allowed: tuple[str, ...] | list[str]) -> None:
    """Refuse a merged table whose columns do not all belong to a sample.

    A merged deposit holds that study's cells and nothing else, so a column
    naming a library the map does not know is a broken map, not a missing
    sample. The distinction matters because the two look identical downstream:
    GSE155622's counts tables spell a library `SNI.6h` where its metadata
    spells it `SNI 6h`, so seven of ten libraries mapped to nothing and the
    study emerged with three control samples and an empty case arm -- alongside
    eight samples that genuinely are absent from the deposit, which made the
    warning look like the ordinary finding it was not.

    A study that really does deposit cells from samples outside the arm names
    those prefixes in `unmapped_prefixes`, so the exclusion is written down.
    """
    surprising = {p: n for p, n in unmapped.items() if p not in allowed}
    if not surprising:
        return
    detail = ", ".join(f"{p!r} ({n} cells)"
                       for p, n in sorted(surprising.items(), key=lambda kv: -kv[1]))
    raise ValueError(
        f"{accession}: {sum(surprising.values())} cell columns name a library "
        f"the sample map does not know: {detail}. Either the map is wrong or "
        f"the deposit carries samples outside the arm, in which case name the "
        f"prefixes in the layout's `unmapped_prefixes`.")


def feature_lookup(accession: str, spec: dict) -> dict[str, str]:
    """The id map a deposit ships alongside its matrices, keyed both ways.

    GSE198608's two submission batches spell the same 32,883 genes two ways:
    the earlier one carries R's syntactic mangling of the id (`5S_rRNA_1`) and
    the later one the id itself (`5-8S-rRNA-1`). Keying on both is what lets
    the two batches meet on one feature index; without it 2,816 genes appear
    in one batch only and are dropped as unmeasured.
    """
    name, from_col, to_col = spec["feature_map"]
    table = pd.read_csv(RAW / accession / name)
    source = table[from_col].astype(str)
    target = table[to_col].astype(str)
    mangled = source.str.replace(r"[^0-9A-Za-z_.]", "_", regex=True)
    return dict(zip(source, target)) | dict(zip(mangled, target))


def apply_feature_map(lookup: dict[str, str], block: Pseudobulk) -> Pseudobulk:
    """Rewrite feature ids through a map the deposit ships alongside."""
    return Pseudobulk([lookup.get(f, f) for f in block.features],
                      block.columns, block.cells)


def align_features(accession: str, parts: list[Pseudobulk]) -> Pseudobulk:
    """Put every sample of a study on one feature index.

    Sample order is not a guarantee even within a study: GSE198608's two
    submission batches deposit the same 32,884 genes in two different orders,
    so a positional merge would silently transpose the whole annotation.
    Features not measured in every sample are dropped rather than filled with
    zero, because an absent row is an absent measurement and a fabricated zero
    is an effect.
    """
    # Duplicates are summed per sample before anything is aligned. A 10x
    # features file routinely gives two Ensembl ids the same symbol, and
    # keeping one occurrence of each would discard the other's counts.
    frames = [aggregate_duplicate_features(
        pd.DataFrame(part.columns, index=pd.Index(part.features, name="feature_id")),
        quiet=True) for part in parts]

    shared: set[str] | None = None
    for frame in frames:
        names = set(frame.index)
        shared = names if shared is None else (shared & names)
    common = sorted(shared or set())

    widest = max(len(f) for f in frames)
    if len(common) < widest:
        logger.warning("  feature sets differ across samples: keeping the %d "
                       "measured in every sample, of %d at the widest",
                       len(common), widest)
    if not common:
        raise ValueError(f"{accession}: no feature is present in every sample")

    columns: dict[str, np.ndarray] = {}
    for frame in frames:
        aligned = frame.loc[common]
        for sample in aligned.columns:
            columns[sample] = aligned[sample].to_numpy()
    cells = {k: v for part in parts for k, v in part.cells.items()}
    return Pseudobulk(common, columns, cells)


def collapse_subjects(block: Pseudobulk, groups: pd.DataFrame) -> Pseudobulk:
    """Sum the columns of GSMs step 18 assigned the same biological unit."""
    subject_of = dict(zip(groups["sample_id"], groups["subject"]))
    out: dict[str, np.ndarray] = {}
    cells: dict[str, tuple[int, int]] = {}
    for gsm, values in block.columns.items():
        subject = subject_of.get(gsm, gsm)
        if subject in out:
            out[subject] = out[subject] + values
            seen, kept = cells[subject]
            cells[subject] = (seen + block.cells[gsm][0], kept + block.cells[gsm][1])
        else:
            out[subject] = values.copy()
            cells[subject] = block.cells[gsm]
    return Pseudobulk(block.features, out, cells)


def aggregate_duplicate_features(frame: pd.DataFrame,
                                 quiet: bool = False) -> pd.DataFrame:
    """Sum rows sharing a feature id, which a symbol map can create."""
    if frame.index.has_duplicates:
        before = len(frame)
        frame = frame.groupby(level=0, sort=False).sum()
        if not quiet:
            logger.info("  %d rows collapsed onto %d distinct features",
                        before, len(frame))
    return frame


def run_study(accession: str, min_umi: int) -> list[dict]:
    spec = LAYOUT[accession]
    groups = pd.read_csv(GROUPS_DIR / f"{accession}_groups.csv")
    logger.info("%s: %d assigned GSM(s), %s deposit",
                accession, len(groups), spec["kind"])

    # The map has to be applied before samples are aligned, not after: the ids
    # it resolves are exactly the ones that otherwise fail to match.
    lookup = feature_lookup(accession, spec) if "feature_map" in spec else None

    if spec["kind"] == "dense_merged":
        block = read_merged(accession, spec, groups, min_umi)
        if lookup:
            block = apply_feature_map(lookup, block)
        missing = [g for g in groups["sample_id"] if g not in block.columns]
    else:
        parts: list[Pseudobulk] = []
        missing = []
        for gsm in groups["sample_id"]:
            part = read_gsm(accession, gsm, spec, min_umi)
            if part is None:
                missing.append(gsm)
                continue
            parts.append(apply_feature_map(lookup, part) if lookup else part)
        if not parts:
            logger.error("  no sample of %s could be read", accession)
            return [{"accession": accession, "sample_id": g, "status": "absent_from_deposit"}
                    for g in groups["sample_id"]]
        block = align_features(accession, parts)

    for gsm in missing:
        logger.warning("  %s: assigned to an arm but absent from the deposit", gsm)

    # Per-GSM cell counts have to be taken before the collapse: it re-keys the
    # dictionary by subject, so a later lookup by GSM answers "no match" for
    # exactly the studies that have a repeated animal, and records them as
    # zero cells. GSE162807 read as 0 of 0 barcodes across all 28 samples.
    cells_per_gsm = dict(block.cells)
    block = collapse_subjects(block, groups)

    frame = pd.DataFrame(block.columns, index=pd.Index(block.features, name="feature_id"))
    frame = aggregate_duplicate_features(frame)
    frame = frame.loc[frame.sum(axis=1) > 0]
    frame = frame.round().astype("int64")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out = OUT_DIR / f"{accession}_counts.tsv.gz"
    frame.to_csv(out, sep="\t", compression="gzip")
    logger.info("  %d features x %d samples -> %s",
                len(frame), frame.shape[1], out.name)

    return manifest_rows(accession, groups, cells_per_gsm, frame, missing)


def manifest_rows(accession: str, groups: pd.DataFrame,
                  cells_per_gsm: dict[str, tuple[int, int]],
                  frame: pd.DataFrame, missing: list[str]) -> list[dict]:
    """One manifest row per assigned GSM, placed or absent.

    `cells_per_gsm` is keyed by GSM and `frame` by subject, and the two differ
    exactly where an animal was sequenced twice. Reading either through the
    other's key returns nothing and records a zero, which is indistinguishable
    from a sample that genuinely contributed no cells.
    """
    subject_of = dict(zip(groups["sample_id"], groups["subject"]))
    group_of = dict(zip(groups["sample_id"], groups["group"]))
    rows: list[dict] = []
    for gsm in groups["sample_id"]:
        subject = subject_of.get(gsm, gsm)
        if gsm in missing:
            rows.append({"accession": accession, "sample_id": gsm,
                         "subject": subject, "group": group_of.get(gsm, ""),
                         "cells_seen": 0, "cells_kept": 0, "total_counts": 0,
                         "status": "absent_from_deposit"})
            continue
        seen, kept = cells_per_gsm.get(gsm, (0, 0))
        rows.append({"accession": accession, "sample_id": gsm, "subject": subject,
                     "group": group_of.get(gsm, ""), "cells_seen": seen,
                     "cells_kept": kept,
                     "total_counts": int(frame[subject].sum()) if subject in frame else 0,
                     "status": "pseudobulked"})
    return rows


def merged_manifest() -> pd.DataFrame:
    """The manifest, assembled from every run's shard.

    A later run of the same accession supersedes an earlier one, so the shards
    are read newest last and duplicates resolved on the sample.
    """
    shards = sorted(SHARD_DIR.glob("*.csv"), key=lambda p: p.stat().st_mtime)
    frames = [pd.read_csv(s) for s in shards]
    if not frames:
        return pd.DataFrame()
    frame = pd.concat(frames, ignore_index=True)
    frame = frame.drop_duplicates(subset=["accession", "sample_id"], keep="last")
    return frame.sort_values(["accession", "sample_id"], ignore_index=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--accession", action="append")
    parser.add_argument("--min-umi-per-cell", type=int, default=DEFAULT_MIN_UMI)
    args = parser.parse_args()

    wanted = args.accession or sorted(LAYOUT)
    rows: list[dict] = []
    for accession in wanted:
        if not (GROUPS_DIR / f"{accession}_groups.csv").exists():
            logger.warning("%s: no group assignment, run step 18 first", accession)
            continue
        rows.extend(run_study(accession, args.min_umi_per_cell))

    if rows:
        frame = pd.DataFrame(rows)
        INTERIM.mkdir(parents=True, exist_ok=True)
        # One shard file per run, merged on read. The manifest is the only
        # record of how many cells each sample carried, and a read-modify-write
        # of a single file loses a whole study's worth of that whenever two
        # runs overlap -- which they will, because the studies differ by an
        # order of magnitude in cost and running the slow one alongside the
        # rest is the obvious thing to do.
        SHARD_DIR.mkdir(parents=True, exist_ok=True)
        shard = SHARD_DIR / (",".join(sorted(frame["accession"].unique())) + ".csv")
        frame.to_csv(shard, index=False)
        frame = merged_manifest()
        frame.to_csv(MANIFEST, index=False)
        absent = frame[frame["status"] == "absent_from_deposit"]
        logger.info("%d sample(s) across %d stud%s -> %s",
                    len(frame), frame["accession"].nunique(),
                    "y" if frame["accession"].nunique() == 1 else "ies", MANIFEST)
        if len(absent):
            logger.warning("%d assigned sample(s) absent from their deposit: %s",
                           len(absent), ", ".join(absent["sample_id"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

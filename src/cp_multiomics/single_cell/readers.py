"""Readers that turn one deposited single-cell matrix into per-gene totals.

Every reader returns the same thing: gene identifiers, and one column of
summed counts per sample. Summing across the cells of a sample is the whole
point -- a per-sample pseudobulk profile is a bulk profile, so the sample
rather than the cell is the unit of replication and the effect sizes stay
comparable with the bulk arm.

Four deposit shapes are covered, because the ten studies of the arm use four:

    ``read_10x_h5``      CellRanger ``filtered_feature_bc_matrix.h5``
    ``read_mtx``         Matrix Market plus feature and barcode sidecars, in
                         either orientation
    ``read_dense``       a delimited genes-by-cells table, read in row chunks
                         because two of them are ~19k genes by ~36k cells
    ``read_dense_split`` the same, split into several samples by a prefix on
                         each cell's column name

None of them holds a full matrix in memory when it can avoid it. The
pseudobulk of a sample is a row sum, so the cell axis is reduced as it is
read.

A UMI floor is applied uniformly. Some of these deposits are CellRanger's
filtered output and some are the raw droplet matrix -- GSE254360 deposits
2.4 million barcodes for 2,400 cells -- and summing every droplet of a raw
matrix adds the ambient profile to the sample. The floor is a cell call, not a
quality filter: on an already-filtered matrix it removes nothing.
"""

from __future__ import annotations

import gzip
import logging
import re
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class Pseudobulk:
    """Per-gene summed counts for one or more samples of one study."""

    features: list[str]
    #: sample id -> summed counts, aligned to `features`
    columns: dict[str, np.ndarray]
    #: sample id -> (barcodes seen, barcodes kept)
    cells: dict[str, tuple[int, int]]


def _open(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") \
        if path.suffix == ".gz" else open(path, encoding="utf-8", errors="replace")


def _reduce(rows: np.ndarray, cols: np.ndarray, data: np.ndarray,
            n_genes: int, n_cells: int, min_umi: int) -> tuple[np.ndarray, int]:
    """Sum `data` over the cells whose total reaches `min_umi`."""
    per_cell = np.bincount(cols, weights=data, minlength=n_cells)
    keep = per_cell >= min_umi
    if not keep.all():
        mask = keep[cols]
        rows, data = rows[mask], data[mask]
    totals = np.bincount(rows, weights=data, minlength=n_genes)
    return totals, int(keep.sum())


def read_10x_h5(path: Path, sample_id: str, min_umi: int) -> Pseudobulk:
    """CellRanger HDF5. Genes are rows; the CSC pointer runs over barcodes."""
    import h5py

    with h5py.File(path, "r") as handle:
        group = handle["matrix"] if "matrix" in handle else handle[list(handle)[0]]
        shape = group["shape"][:]
        data = group["data"][:].astype(np.float64)
        indices = group["indices"][:]
        indptr = group["indptr"][:]
        feats = group["features"]
        names = feats["name"][:] if "name" in feats else feats["id"][:]
        features = [n.decode() for n in names]
        # A multiome deposit carries Peaks alongside Gene Expression in one
        # file, and a peak is not a gene. GSE253345 is the reason this is here.
        if "feature_type" in feats:
            kinds = np.array([t.decode() for t in feats["feature_type"][:]])
            gene_rows = np.flatnonzero(kinds == "Gene Expression")
        else:
            gene_rows = np.arange(len(features))

    n_genes, n_cells = int(shape[0]), int(shape[1])
    cols = np.repeat(np.arange(n_cells), np.diff(indptr))
    totals, kept = _reduce(indices, cols, data, n_genes, n_cells, min_umi)
    return Pseudobulk([features[i] for i in gene_rows],
                      {sample_id: totals[gene_rows]},
                      {sample_id: (n_cells, kept)})


def read_mtx(mtx: Path, features: list[str], sample_id: str, min_umi: int,
             *, cells_are_rows: bool) -> Pseudobulk:
    """Matrix Market counts, streamed as coordinate triplets.

    `cells_are_rows` because the orientation is not a convention: 10x writes
    genes by cells, while the Parse and kallisto outputs in this arm write
    cells by genes and say so only in a comment or a filename.
    """
    with _open(mtx) as handle:
        for line in handle:
            if line.startswith("%"):
                continue
            dim_a, dim_b, _ = (int(v) for v in line.split())
            break
        triplets = np.loadtxt(handle, dtype=np.float64)

    if triplets.ndim == 1:                       # a matrix with one entry
        triplets = triplets.reshape(1, 3)
    a = triplets[:, 0].astype(np.int64) - 1      # Matrix Market is 1-based
    b = triplets[:, 1].astype(np.int64) - 1
    values = triplets[:, 2]
    genes, cells = (b, a) if cells_are_rows else (a, b)
    n_genes, n_cells = (dim_b, dim_a) if cells_are_rows else (dim_a, dim_b)

    if n_genes != len(features):
        raise ValueError(f"{mtx.name}: {n_genes} gene rows against "
                         f"{len(features)} feature names")
    totals, kept = _reduce(genes, cells, values, n_genes, n_cells, min_umi)
    return Pseudobulk(features, {sample_id: totals}, {sample_id: (n_cells, kept)})


def read_dense(path: Path, sample_id: str, min_umi: int,
               sep: str = "\t") -> Pseudobulk:
    """A genes-by-cells table, whose every column belongs to one sample."""
    return read_dense_split(path, lambda _column: sample_id, min_umi, sep=sep)


def read_dense_split(path: Path, sample_of: Callable[[str], str | None],
                     min_umi: int,
                     sep: str = "\t") -> Pseudobulk:
    """A genes-by-cells table carrying several samples, split by `sample_of`.

    Read in row chunks: the two tables this serves are ~19k genes tall by 19k
    and 63k cells, which is 350 to 1,200 million values, and materialising one
    of those as a frame costs several gigabytes to produce a few columns.

    Column positions are resolved once against the header, never per chunk.
    `Index.get_loc` on a 63,000-entry object index is a scan, so calling it
    inside the chunk loop is quadratic in the cell count and turns a two-minute
    read into an unfinished one.
    """
    header = pd.read_csv(path, sep=sep, nrows=0)
    columns = list(header.columns)[1:]
    owners = [(i, c, sample_of(c)) for i, c in enumerate(columns)]
    known = [(i, c, o) for i, c, o in owners if o]
    if not known:
        raise ValueError(f"{path.name}: no column mapped to a sample")

    samples = sorted({o for _i, _c, o in known})
    take = {s: np.array([i for i, _c, o in known if o == s], dtype=np.int64)
            for s in samples}
    names = {s: [c for _i, c, o in known if o == s] for s in samples}

    features: list[str] = []
    per_cell = {s: np.zeros(len(take[s])) for s in samples}
    blocks: dict[str, list[np.ndarray]] = {s: [] for s in samples}
    for chunk in pd.read_csv(path, sep=sep, index_col=0, chunksize=2000):
        features.extend(str(i) for i in chunk.index)
        values = chunk.to_numpy(dtype=np.float64, na_value=0.0)
        for sample in samples:
            mine = values[:, take[sample]]
            per_cell[sample] += mine.sum(axis=0)
            blocks[sample].append(mine)

    result: dict[str, np.ndarray] = {}
    counts: dict[str, tuple[int, int]] = {}
    for sample in samples:
        keep = per_cell[sample] >= min_umi
        result[sample] = np.vstack(blocks[sample])[:, keep].sum(axis=1)
        counts[sample] = (len(names[sample]), int(keep.sum()))
    return Pseudobulk(features, result, counts)


def read_lines(path: Path, column: int = 0, sep: str = "\t",
               skip_header: bool = False) -> list[str]:
    """One field per line, for the feature and barcode sidecars of an mtx."""
    out: list[str] = []
    with _open(path) as handle:
        if skip_header:
            next(handle, None)
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split(sep)
            out.append(parts[column] if column < len(parts) else parts[0])
    return out


def strip_version(identifiers: list[str]) -> list[str]:
    """`ENSMUSG00000102693.1` -> `ENSMUSG00000102693`."""
    return [re.sub(r"\.\d+$", "", i) for i in identifiers]

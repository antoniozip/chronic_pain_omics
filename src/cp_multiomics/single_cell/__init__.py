"""Single-cell arm: pseudobulk aggregation of deposited cell-level matrices."""

from cp_multiomics.single_cell.readers import (
    Pseudobulk,
    read_10x_h5,
    read_dense,
    read_dense_split,
    read_lines,
    read_mtx,
    strip_version,
)

__all__ = [
    "Pseudobulk",
    "read_10x_h5",
    "read_dense",
    "read_dense_split",
    "read_lines",
    "read_mtx",
    "strip_version",
]

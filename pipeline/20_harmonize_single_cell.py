"""Write the harmonized manifest for the single-cell arm.

Steps 05, 05b and 06 all key on `data/interim/{modality}/harmonized.jsonl`,
so the arm needs one under its own modality name to travel the existing path.
Writing it here rather than by hand keeps the species, sample counts and id
space derived from what steps 18 and 19 actually produced, instead of from the
registry's description of what the deposits claim.

`single_cell` is deliberately not one of step 05's or step 06's
`ALL_MODALITIES`. That is the mechanism keeping the arm out of the bulk pool:
`--modality all` cannot reach it, so merging the two requires someone to type
the name, and there is no path by which it happens silently.

Usage:
    python pipeline/20_harmonize_single_cell.py
"""

from __future__ import annotations

import argparse
import gzip
import json
import logging
import re
from datetime import UTC, datetime
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                    datefmt="%H:%M:%S")
logger = logging.getLogger(__name__)

REPO_ROOT = Path(__file__).resolve().parent.parent
REGISTRY = REPO_ROOT / "conf" / "analysis" / "single_cell_studies.csv"
CANDIDATES = REPO_ROOT / "literature" / "prisma" / "single_cell_candidates.csv"
INTERIM = REPO_ROOT / "data" / "interim" / "single_cell"
PSEUDOBULK = INTERIM / "pseudobulk"
MANIFEST = INTERIM / "harmonized.jsonl"

#: Species as the registry spells them, canonicalised. GSE162807 is deposited
#: as "Homo sapiens; Mus musculus" and enters mouse-only: its three human
#: samples are all naive, so step 18 excluded them, and recording the study as
#: multi-species here would make step 05b refuse to resolve its Ensembl ids and
#: step 06g drop it from the species pools for a human half it no longer has.
SPECIES_OVERRIDE = {"GSE162807": ["Mus musculus"]}

ENSEMBL_PREFIX = re.compile(r"^ENS(MUS|RNO)?[GT]\d")


def id_space_of(counts: Path) -> str:
    """What the pseudobulk table's feature ids are, read from the table."""
    with gzip.open(counts, "rt") as handle:
        next(handle)
        sample = [next(handle, "").split("\t", 1)[0] for _ in range(200)]
    sample = [s for s in sample if s]
    if not sample:
        return "unknown"
    if sum(bool(ENSEMBL_PREFIX.match(s)) for s in sample) / len(sample) >= 0.95:
        return "ensembl"
    return "symbol"


def main() -> int:
    argparse.ArgumentParser(description=__doc__).parse_args()

    registry = pd.read_csv(REGISTRY).set_index("accession")
    candidates = pd.read_csv(CANDIDATES)
    included = candidates[candidates["verdict"] == "include"]["accession"]

    records = []
    for accession in sorted(included):
        counts = PSEUDOBULK / f"{accession}_counts.tsv.gz"
        if not counts.exists():
            logger.warning("%s: no pseudobulk table, run step 19 first", accession)
            continue
        meta = registry.loc[accession]
        header = pd.read_csv(counts, sep="\t", nrows=0)
        samples = [c for c in header.columns if c != "feature_id"]
        species = SPECIES_OVERRIDE.get(
            accession, [s.strip() for s in str(meta["species"]).split(";")])

        records.append({
            "modality": "single_cell",
            "source": "GEO",
            "accession": accession,
            "title": meta["title"],
            "species_canonical": species,
            "has_human": "Homo sapiens" in species,
            "has_animal": any(s != "Homo sapiens" for s in species),
            # Pseudobulk counts are integer counts of reads, so the per-study
            # DA runs the RNA-seq branch and not the microarray one.
            "platform_canonical": "RNA-seq",
            "year": str(meta["year"]),
            "n_samples": len(samples),
            "assay": meta["assay"],
            "id_space": id_space_of(counts),
            "effect_size_unit": "log2FC",
            "url": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}",
            "ortholog_needed": species != ["Homo sapiens"],
            "quality_flags": [],
            "harmonized_at": datetime.now(UTC).isoformat(),
        })
        logger.info("%-10s %-18s %2d samples, ids in %s space", accession,
                    "/".join(species), len(samples), records[-1]["id_space"])

    INTERIM.mkdir(parents=True, exist_ok=True)
    with open(MANIFEST, "w") as handle:
        for record in records:
            handle.write(json.dumps(record) + "\n")
    logger.info("%d studies -> %s", len(records), MANIFEST)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

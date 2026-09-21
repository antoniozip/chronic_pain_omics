"""Guard the SPARC screening record and the identifier space around it.

SPARC is the first repository in this project whose accessions are not GEO
series, and a mixed identifier space is how several silent failures here have
started: an id that matches nothing answers "no match", and every caller reads
that as "nothing to do". The record therefore lives in its own sheet, and these
tests hold the two properties that keep it harmless.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
PRISMA = REPO_ROOT / "literature" / "prisma"
SPARC_SHEET = PRISMA / "sparc_candidates.csv"

sys.path.insert(0, str(REPO_ROOT / "scripts"))

import fetch_sparc_dataset as fsd  # noqa: E402


def rows() -> list[dict[str, str]]:
    with SPARC_SHEET.open() as fh:
        return list(csv.DictReader(fh))


def test_every_row_records_a_verdict_and_a_reason():
    # PRISMA discipline: never drop a study without saying why.
    for row in rows():
        assert row["verdict"], row
        assert row["reason_code"], row
        assert row["reason_source"], row


def test_no_sparc_id_leaks_into_a_geo_keyed_sheet():
    # load_included_accessions() reads {modality}_candidates.csv and matches
    # against harmonized GEO accessions. A SPARC id there would match nothing,
    # which is indistinguishable from a study with no data.
    ids = {row["accession"] for row in rows()}
    for sheet in sorted(PRISMA.glob("*_candidates.csv")) + \
            sorted(PRISMA.glob("*_manual_review.csv")):
        if sheet == SPARC_SHEET:
            continue
        text = sheet.read_text()
        for accession in ids:
            assert accession not in text, f"{accession} leaked into {sheet.name}"


def test_the_sheet_is_not_named_for_a_modality():
    # load_included_accessions(prisma_dir, modality) globs for
    # "{modality}_candidates.csv". Naming this sheet after a modality would
    # silently substitute it for that modality's real screening record.
    modalities = {"genomics", "transcriptomics", "proteomics", "metabolomics",
                  "lipidomics", "single_cell"}
    assert SPARC_SHEET.stem.replace("_candidates", "") not in modalities


def test_an_included_row_fails_until_something_can_read_it():
    # Nothing reads this sheet yet. Flipping a verdict to "include" would look
    # like an inclusion and change nothing, which is the failure mode this
    # whole file exists to prevent. When a SPARC study is genuinely admitted,
    # wire it into ingestion first and then relax this test.
    included = [row["accession"] for row in rows() if row["verdict"] == "include"]
    assert not included, (
        f"{included} are marked include, but no pipeline step reads "
        f"{SPARC_SHEET.name}. Add an ingestion route before including a study, "
        "or the inclusion is silently a no-op.")


def test_object_key_prefers_the_listing_uri():
    # The rebuilt key is a guess that happens to be right for the datasets seen
    # so far; the listing's own URI is authoritative.
    uri = "s3://sparc-prod-aod-discover-publish50-use1/478/files/primary/a.csv"
    assert fsd.object_key(uri, 478, "files/primary/a.csv") == \
        "478/files/primary/a.csv"


def test_object_key_falls_back_when_the_listing_carries_no_uri():
    assert fsd.object_key("", 478, "files/primary/a.csv") == \
        "478/files/primary/a.csv"


@pytest.mark.parametrize("column", ["accession", "repository", "doi", "url"])
def test_every_row_can_be_traced_back_to_its_deposit(column: str):
    for row in rows():
        assert row[column].strip(), (column, row["accession"])

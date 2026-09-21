#!/usr/bin/env python
"""Write the reviewer's verdicts into the Workbench PRISMA candidate record.

`validate_workbench_candidates.py` assembles evidence and proposes; it does not
decide. This records what was decided on 2026-08-27, so that
`metabolomics_workbench_candidates.csv` carries a verdict and a reason code for
every retrieved study, as `metabolomics_manual_review.csv` does for
MetaboLights.

Twelve of the 23 candidates clear every mechanical gate (a control arm of at
least three against three, not in vitro, abundances actually deposited). Six
were included; the exclusions divide into two kinds, and the distinction
matters for anyone re-reading this record.

**Evidence-based exclusions** are properties of the data. ST001327 and ST001273
model sonication-induced traumatic optic neuropathy, whose abstracts describe
vision loss and retinal ganglion cell death with no pain measure. The rest fail
a mechanical gate and carry that gate's code.

**One scope decision** is not a property of the data. ST003177, ST003520,
ST001949 and ST000218 are rheumatoid, juvenile idiopathic and psoriatic
arthritis cohorts with control arms and, in three cases, excellent annotation
depth — ST003177 alone is 2,492 samples. They were excluded to keep the arm's
phenotypes centred on conditions defined by pain.

**That decision sits awkwardly beside MTBLS5667**, a knee-osteoarthritis study
already among the three analysable MetaboLights studies. Arthritis is therefore
in the corpus on one side of the line and out on the other. The code
`phenotype_scope_decision` records this as a choice rather than a finding, so
it can be revisited without re-deriving the evidence; the evidence itself is in
`metabolomics_workbench_validation.csv` and does not change either way.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
CANDIDATES = REPO_ROOT / "literature" / "prisma" / "metabolomics_workbench_candidates.csv"

INCLUDE = {
    "ST000676": "usable_pending_file_check",   # HF/STZ neuropathy; DRG + sciatic nerve
    "ST001412": "usable_pending_file_check",   # obese neuropathy v obese non-neuropathy
    "ST003954": "usable_pending_file_check",   # IBS v HC; matches MTBLS2774's phenotype
    "ST000603": "usable_pending_file_check",   # interstitial cystitis / bladder pain
    "ST003984": "usable_pending_file_check",   # endometriosis; chronic pelvic pain
    "ST000780": "usable_pending_file_check",   # HF diet, sciatic nerve function
}

EXCLUDE = {
    # Mechanical gates, from metabolomics_workbench_validation.csv.
    "ST000876": "no_control_arm_in_factors",
    "ST001121": "no_control_arm_in_factors",
    "ST001122": "no_control_arm_in_factors",
    "ST001815": "no_control_arm_in_factors",
    "ST001940": "intervention_pre_post_not_case_control",
    "ST004518": "no_pain_versus_control_contrast",
    "ST002540": "no_pain_versus_control_contrast",
    "ST002974": "no_pain_versus_control_contrast",
    "ST003523": "no_pain_versus_control_contrast",
    "ST000585": "no_abundances_deposited",
    "ST003579": "no_abundances_deposited",
    # Phenotype, from the abstracts: optic neuropathy is vision loss.
    "ST001327": "no_pain_phenotype",
    "ST001273": "no_pain_phenotype",
    # Scope decision, not a property of the data. See module docstring.
    "ST003177": "phenotype_scope_decision",
    "ST003520": "phenotype_scope_decision",
    "ST001949": "phenotype_scope_decision",
    "ST000218": "phenotype_scope_decision",
}


def main() -> int:
    with open(CANDIDATES) as f:
        rows = list(csv.DictReader(f))
        fields = list(rows[0])

    unknown = {r["accession"] for r in rows} - set(INCLUDE) - set(EXCLUDE)
    if unknown:
        raise SystemExit(f"no verdict recorded for {sorted(unknown)}")

    for r in rows:
        acc = r["accession"]
        r["verdict"] = "include" if acc in INCLUDE else "exclude"
        r["reason_code"] = INCLUDE.get(acc) or EXCLUDE[acc]
        r["reason_source"] = "reviewer+rest_validation"
        r["reason_code_proposed"] = ""

    with open(CANDIDATES, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    n_inc = sum(1 for r in rows if r["verdict"] == "include")
    print(f"-> {CANDIDATES.relative_to(REPO_ROOT)}  {len(rows)} rows, {n_inc} include")
    return 0


if __name__ == "__main__":
    sys.exit(main())

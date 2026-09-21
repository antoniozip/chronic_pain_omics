"""
Analyze the 500 PRIDE proteomics datasets for chronic pain relevance.
Searches titles, descriptions, species, and keywords for pain vocabulary.
Checks PRIDE API for quantification data availability.
"""

import json
import ssl
import time
import urllib.error
import urllib.request
from collections import Counter
from pathlib import Path

from cp_multiomics.vocabulary import load_pain_vocabulary

WORKDIR = Path("/media/antonio/data/CP_multi-omics")
HARMONIZED_FILE = WORKDIR / "data/interim/proteomics/harmonized.jsonl"
RAW_FILE = WORKDIR / "data/raw/proteomics/datasets.jsonl"
PLAN_DIR = WORKDIR / "plan"
PLAN_DIR.mkdir(parents=True, exist_ok=True)

# Pain-related keywords, loaded from the single source of truth in
# conf/search/pain_vocabulary.yaml (see src/cp_multiomics/vocabulary.py).
PAIN_KEYWORDS = list(load_pain_vocabulary().all_terms)

# Species we care about
RELEVANT_SPECIES = {"Homo sapiens (human)", "Mus musculus (mouse)", "Rattus norvegicus (rat)"}

# Harmful species (irrelevant)
HARMFUL_SPECIES = [
    "Listeria", "Escherichia", "Salmonella", "Staphylococcus", "Streptococcus",
    "Arabidopsis", "Drosophila", "Caenorhabditis", "Saccharomyces",
    "Zea mays", "Oryza sativa", "Triticum", "Glycine max",
    "Sus scrofa", "Bos taurus", "Gallus gallus",
]


def load_jsonl(path):
    records = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line:
                records.append(json.loads(line))
    return records


def search_text(text, keywords):
    """Return set of matched keywords in text (case-insensitive)."""
    text_lower = text.lower()
    matches = set()
    for kw in keywords:
        if kw.lower() in text_lower:
            matches.add(kw)
    return matches


def score_pain_relevance(record):
    """
    Score a dataset for chronic pain relevance based on title, description, keywords.
    Returns (score, matched_terms, reasoning).
    """
    title = record.get("title", "")
    description = record.get("description", "")
    species = record.get("species_canonical", record.get("species", []))
    keywords_raw = record.get("keywords", [])

    matched_terms = set()

    # Search title
    title_matches = search_text(title, PAIN_KEYWORDS)
    matched_terms.update(title_matches)

    # Search description
    desc_matches = search_text(description, PAIN_KEYWORDS)
    matched_terms.update(desc_matches)

    # Search keywords
    for kw in keywords_raw:
        kw_matches = search_text(kw, PAIN_KEYWORDS)
        matched_terms.update(kw_matches)

    # Check species relevance
    species_set = set(species)
    relevant_species = species_set & RELEVANT_SPECIES

    # Count unique matched terms
    n_matches = len(matched_terms)

    # Heuristic scoring
    title_boost = sum(1 for kw in PAIN_KEYWORDS if kw.lower() in title.lower()) * 3

    score = n_matches + title_boost

    # Boost for strong pain-specific terms in title
    strong_pain_terms = ["chronic pain", "neuropathic pain", "inflammatory pain", "nerve injury",
                         "dorsal root ganglion", "DRG", "spinal cord injury", "nociceptor",
                         "hyperalgesia", "allodynia", "migraine", "fibromyalgia", "sciatic",
                         "opioid", "PACAP"]
    for term in strong_pain_terms:
        if term.lower() in title.lower():
            score += 5

    # Require relevant species for high relevance
    if not relevant_species:
        if matched_terms:
            score = max(1, score // 2)  # reduce score if no relevant species

    return score, matched_terms, relevant_species


def check_pride_api(accession):
    """Check what files/quant data are available via PRIDE REST API."""
    url = f"https://www.ebi.ac.uk/pride/ws/archive/v2/projects/{accession}"
    try:
        ctx = ssl.create_default_context()
        ctx.check_hostname = False
        ctx.verify_mode = ssl.CERT_NONE
        req = urllib.request.Request(url, headers={"Accept": "application/json"})
        with urllib.request.urlopen(req, context=ctx, timeout=15) as resp:
            data = json.loads(resp.read().decode())
        return data
    except Exception as e:
        return {"error": str(e)}


def check_pride_files(accession):
    """Check available assay files for a PRIDE project."""
    url = f"https://www.ebi.ac.uk/pride/ws/archive/v2/projects/{accession}/files"
    try:
        ctx = ssl.create_default_context()
        ctx.check_hostname = False
        ctx.verify_mode = ssl.CERT_NONE
        req = urllib.request.Request(url, headers={"Accept": "application/json"})
        with urllib.request.urlopen(req, context=ctx, timeout=15) as resp:
            data = json.loads(resp.read().decode())
        return data
    except Exception as e:
        return {"error": str(e)}


def main():
    # Load datasets
    harmonized = load_jsonl(HARMONIZED_FILE)
    raw = load_jsonl(RAW_FILE)

    print(f"Loaded {len(harmonized)} harmonized records, {len(raw)} raw records")
    print()

    # Build raw lookup by accession
    raw_by_accession = {r["accession"]: r for r in raw}

    # --- Screen all datasets for pain relevance ---
    scored = []
    for rec in harmonized:
        accession = rec["accession"]
        score, matched, species = score_pain_relevance(rec)
        raw_rec = raw_by_accession.get(accession, {})
        # Get keywords from raw
        keywords = raw_rec.get("keywords", [])

        if score > 0:
            scored.append({
                "accession": accession,
                "title": rec.get("title", ""),
                "score": score,
                "matched_terms": matched,
                "species": species,
                "species_canonical": rec.get("species_canonical", []),
                "has_human": rec.get("has_human", False),
                "has_animal": rec.get("has_animal", False),
                "n_samples": rec.get("n_samples", -1),
                "year": rec.get("year", ""),
                "quality_flags": rec.get("quality_flags", []),
                "keywords": keywords,
                "description": raw_rec.get("description", "")[:200],
            })

    # Sort by score descending
    scored.sort(key=lambda x: x["score"], reverse=True)

    print("=" * 100)
    print("STUDIES WITH PAIN-RELEVANT MATCHES")
    print("=" * 100)
    for s in scored:
        print(f"\n{s['accession']} | Score: {s['score']} | {s['year']} | {s['n_samples']} samples")
        print(f"  Title: {s['title'][:120]}")
        print(f"  Species: {', '.join(s['species_canonical'])}")
        print(f"  Matched: {', '.join(sorted(s['matched_terms']))}")
        print(f"  Flags: {s['quality_flags']}")

    print()
    print("=" * 100)
    print("SUMMARY")
    print("=" * 100)
    print(f"Total datasets screened: {len(harmonized)}")
    print(f"Datasets with pain-relevant matches: {len(scored)}")

    # Category breakdown
    high_confidence = [s for s in scored if s["score"] >= 10]
    medium_confidence = [s for s in scored if 5 <= s["score"] < 10]
    low_confidence = [s for s in scored if s["score"] < 5]

    print(f"  High confidence (score >= 10): {len(high_confidence)}")
    for s in high_confidence:
        print(f"    {s['accession']}: {s['title'][:100]} "
              f"(matched: {', '.join(sorted(s['matched_terms']))})")
    print(f"  Medium confidence (5-9): {len(medium_confidence)}")
    for s in medium_confidence:
        print(f"    {s['accession']}: {s['title'][:100]} "
              f"(matched: {', '.join(sorted(s['matched_terms']))})")
    print(f"  Low confidence (1-4): {len(low_confidence)}")

    # Species breakdown
    species_counts = Counter()
    for rec in harmonized:
        for s in rec.get("species_canonical", []):
            species_counts[s] += 1
    print("\nSpecies distribution (all 500):")
    for s, c in species_counts.most_common(20):
        print(f"  {s}: {c}")

    # Relevant species breakdown
    human_studies = [rec for rec in harmonized if rec.get("has_human")]
    mouse_studies = [rec for rec in harmonized
                     if rec.get("has_animal")
                     and "Mus musculus" in str(rec.get("species_canonical", []))]
    rat_studies = [rec for rec in harmonized
                   if "Rattus norvegicus" in str(rec.get("species_canonical", []))]
    print("\nStudies with relevant species:")
    print(f"  Human: {len(human_studies)}")
    print(f"  Mouse: {len(mouse_studies)}")
    print(f"  Rat: {len(rat_studies)}")

    # --- Check PRIDE API for top candidates ---
    print("\n" + "=" * 100)
    print("PRIDE API CHECK - What data is available?")
    print("=" * 100)

    # Check API structure first with a known accession
    probe_accessions = ["PXD013362", "PXD015949", "PXD003830", "PXD006512", "PXD011963"]
    for acc in probe_accessions[:3]:  # Check first 3
        print(f"\n--- {acc} ---")
        data = check_pride_api(acc)
        if "error" in data:
            print(f"  API Error: {data['error']}")
        else:
            print(f"  Title: {data.get('title', 'N/A')}")
            print(f"  Organisms: {data.get('organisms', [])}")
            print(f"  Files: {len(data.get('files', []))} files listed")
            print(f"  Sample count: {data.get('samplesCount', '?')}")
            print(f"  Assay count: {data.get('assayCount', '?')}")
            # Check for quant files
            files = data.get("files", [])
            for f in files[:20]:
                if not isinstance(f, dict):
                    continue
                fname = f.get("fileName", f.get("name", ""))
                ft = f.get("fileType")
                if isinstance(ft, dict):
                    ftype = ft.get("name", "")
                else:
                    ftype = str(ft) if ft else ""
                if any(q in fname.lower() for q in
                       ["quant", "intensity", "maxquant", "proteingroups",
                        "peptide", "protein"]):
                    print(f"    -> {fname} ({ftype})")

        # Check files endpoint
        files_data = check_pride_files(acc)
        if "error" not in files_data:
            file_list = (files_data if isinstance(files_data, list)
                         else files_data.get("_embedded", {}).get("files", [])
                         if isinstance(files_data, dict) else [])
            print(f"  Files endpoint returned {len(file_list)} items")
            quant_files = []
            for f in file_list[:30]:
                if isinstance(f, dict):
                    fname = f.get("fileName", f.get("name", ""))
                    ft = f.get("fileType")
                    if isinstance(ft, dict):
                        ftype = ft.get("name", "")
                    else:
                        ftype = str(ft) if ft else ""
                    if any(q in fname.lower() for q in
                           ["quant", "maxquant", "proteingroups", "protein",
                            "intensity", "peptide"]):
                        quant_files.append(fname)
                        print(f"    QUANT: {fname} ({ftype}), size={f.get('fileSizeBytes', 0)}")
            if not quant_files:
                non_raw = [f.get("fileName", f.get("name", "")[:80])
                           for f in file_list[:15] if isinstance(f, dict)]
                print(f"    No quant files found in first 15. Sample files: {non_raw}")

        time.sleep(0.5)  # Rate limiting

    # --- Save results ---
    output = {
        "total_datasets": len(harmonized),
        "pain_relevant_total": len(scored),
        "high_confidence": [
            {"accession": s["accession"], "title": s["title"], "score": s["score"],
             "matched_terms": list(s["matched_terms"]), "species": list(s["species_canonical"]),
             "n_samples": s["n_samples"], "year": s["year"]}
            for s in high_confidence
        ],
        "medium_confidence": [
            {"accession": s["accession"], "title": s["title"], "score": s["score"],
             "matched_terms": list(s["matched_terms"]), "species": list(s["species_canonical"]),
             "n_samples": s["n_samples"], "year": s["year"]}
            for s in medium_confidence
        ],
        "low_confidence": [
            {"accession": s["accession"], "title": s["title"], "score": s["score"],
             "matched_terms": list(s["matched_terms"]), "species": list(s["species_canonical"]),
             "n_samples": s["n_samples"], "year": s["year"]}
            for s in low_confidence
        ],
    }

    # Save full results
    with open(WORKDIR / "data/interim/proteomics/pain_relevance_scores.json", "w") as f:
        json.dump(output, f, indent=2, default=str)
    print("\nResults saved to data/interim/proteomics/pain_relevance_scores.json")

    # Save a filtered subset of pain-relevant harmonized records
    pain_records = []
    for rec in harmonized:
        accession = rec["accession"]
        if any(s["accession"] == accession for s in scored):
            # Add relevance score
            s = next(s for s in scored if s["accession"] == accession)
            rec["pain_relevance_score"] = s["score"]
            rec["pain_matched_terms"] = list(s["matched_terms"])
            pain_records.append(rec)

    with open(WORKDIR / "data/interim/proteomics/pain_relevant.jsonl", "w") as f:
        for rec in pain_records:
            f.write(json.dumps(rec, default=str) + "\n")
    print(f"Saved {len(pain_records)} pain-relevant records to "
          "data/interim/proteomics/pain_relevant.jsonl")

    return scored, output


if __name__ == "__main__":
    scored, output = main()

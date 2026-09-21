"""
Deep-dive investigation of top pain-relevant PRIDE datasets.
Check what quantification data is available and how to access it.
"""
import json
import ssl
import time
import urllib.error
import urllib.request
from pathlib import Path

WORKDIR = Path("/media/antonio/data/CP_multi-omics")

def api_get(url, timeout=20):
    """Fetch JSON from PRIDE API with retry."""
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE
    for attempt in range(3):
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json"})
            with urllib.request.urlopen(req, context=ctx, timeout=timeout) as resp:
                return json.loads(resp.read().decode())
        except Exception as e:
            if attempt == 2:
                return {"error": str(e)}
            time.sleep(1)

# Top pain-relevant datasets
TARGETS = [
    "PXD013362",  # PACAP/chronic migraine/opioid hyperalgesia (mouse) - HIGHEST SCORE
    "PXD015949",  # DRG/dorsal horn neuropeptides (mouse) - HIGH
    "PXD001390",  # CSF proteome mapping (human)
    "PXD002424",  # Human CSF exosome (human)
    "PXD000004",  # Brain cortex proteome (human)
    "PXD002681",  # Rat brain mitochondrial proteins (rat)
    "PXD001195",  # AMPAR proteome across brain regions (rat)
    "PXD001049",  # Nitrated protein profile of rat hypoxic brain (rat)
]

print("=" * 120)
print("PRIDE API DEEP DIVE - Data Availability for Pain-Relevant Studies")
print("=" * 120)

for acc in TARGETS:
    print(f"\n{'─' * 100}")
    print(f"  ACCESSIONSION: {acc}")
    print(f"{'─' * 100}")

    # 1. Project metadata
    data = api_get(f"https://www.ebi.ac.uk/pride/ws/archive/v2/projects/{acc}")
    if "error" in data:
        print(f"  API Error: {data['error']}")
        continue

    print(f"  Title: {data.get('title', 'N/A')}")
    print(f"  Organisms: {[o.get('name','') for o in data.get('organisms', [])]}")
    print(f"  Sample count: {data.get('samplesCount', 'N/A')}")
    print(f"  Assay count: {data.get('assayCount', 'N/A')}")
    print(f"  Project tags: {data.get('projectTags', [])}")
    print(f"  Submission type: {data.get('submissionType', 'N/A')}")

    # 2. Files
    files_data = api_get(f"https://www.ebi.ac.uk/pride/ws/archive/v2/projects/{acc}/files")
    if "error" not in files_data:
        file_list = (files_data if isinstance(files_data, list)
                     else files_data.get("_embedded", {}).get("files", [])
                     if isinstance(files_data, dict) else [])
        print(f"\n  Total files: {len(file_list)}")

        # Categorize files
        raw_files = []
        quant_files = []
        mzid_files = []
        other_files = []

        for f in file_list:
            if not isinstance(f, dict):
                continue
            fname = f.get("fileName", f.get("name", ""))
            public_uri = f.get("publicFileUri", f.get("publicFileUri", ""))
            fsize = f.get("fileSizeBytes", 0)

            if fname.endswith((".raw", ".mzML", ".mzXML", ".wiff", ".d")):
                raw_files.append((fname, fsize, public_uri))
            elif any(q in fname.lower() for q in
                     ["quant", "proteingroups", "maxquant", "peptide",
                      "intensity"]):
                quant_files.append((fname, fsize, public_uri))
            elif fname.endswith(".mzid") or fname.endswith(".mzid.gz"):
                mzid_files.append((fname, fsize, public_uri))
            else:
                other_files.append((fname, fsize, public_uri))

        print(f"  Raw files (.raw/.mzML etc): {len(raw_files)}")
        for fn, sz, uri in raw_files[:5]:
            print(f"    - {fn} ({sz//1024//1024} MB)")
        if len(raw_files) > 5:
            print(f"    ... and {len(raw_files)-5} more")

        print(f"  Quant/result files: {len(quant_files)}")
        for fn, sz, uri in sorted(quant_files, key=lambda x: x[1], reverse=True):
            sz_str = (f"{sz/1024/1024/1024:.1f} GB" if sz > 1e9
                      else f"{sz/1024/1024:.0f} MB" if sz > 1e6
                      else f"{sz/1024:.0f} KB")
            print(f"    - {fn} ({sz_str})")
            if uri:
                print(f"      URI: {uri[:200]}")

        print(f"  MzID files: {len(mzid_files)}")
        for fn, sz, uri in mzid_files[:3]:
            print(f"    - {fn} ({sz/1024:.0f} KB)")

        print(f"  Other files: {len(other_files)}")
        for fn, sz, uri in sorted(other_files, key=lambda x: x[1], reverse=True)[:5]:
            sz_str = (f"{sz/1024/1024/1024:.1f} GB" if sz > 1e9
                      else f"{sz/1024/1024:.0f} MB" if sz > 1e6
                      else f"{sz/1024:.0f} KB")
            print(f"    - {fn} ({sz_str})")
    else:
        print(f"  Files endpoint error: {files_data['error']}")

    # 3. Check how to download
    print("\n  Download info:")
    # Try to get download URL
    dl_data = api_get(f"https://www.ebi.ac.uk/pride/ws/archive/v2/projects/{acc}/download")
    if "error" not in dl_data:
        if isinstance(dl_data, dict):
            print(f"    Download link type: {dl_data.get('type', 'N/A')}")
            for k, v in dl_data.items():
                if isinstance(v, str):
                    print(f"    {k}: {v[:200]}")
    else:
        print(f"    Download endpoint: {dl_data['error']}")

    time.sleep(0.5)

print("\n\nDone.")

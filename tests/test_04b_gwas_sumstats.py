import importlib.util
from pathlib import Path

_SPEC = importlib.util.spec_from_file_location(
    "pipeline_04b", Path(__file__).parent.parent / "pipeline" / "04b_gwas_sumstats.py"
)
pipeline_04b = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(pipeline_04b)


def _cfg(tmp_path):
    return {
        "efo_traits": ["back pain", "fibromyalgia"],
        "ftp_base": "https://ftp.example/ss",
        "gwas_rest": "https://rest.example",
        "output": {
            "raw_dir": str(tmp_path / "raw"),
            "interim_dir": str(tmp_path / "interim"),
            "manifest": str(tmp_path / "manifest.csv"),
            "max_file_mb": 1500,
        },
    }


_REST = {
    "back pain": [{"accessionId": "GCST001", "fullPvalueSet": True,
                   "diseaseTrait": {"trait": "Back pain"},
                   "initialSampleSize": "5,000 UK Biobank cases, 100,000 controls"},
                  {"accessionId": "GCST002", "fullPvalueSet": True,
                   "diseaseTrait": {"trait": "Back pain (Gene-based burden)"},
                   "initialSampleSize": "400,000 UK Biobank"}],
    "fibromyalgia": [{"accessionId": "GCST003", "fullPvalueSet": True,
                      "diseaseTrait": {"trait": "Fibromyalgia"},
                      "initialSampleSize": "1,362 European cases, 44,047 controls"}],
}


def test_dry_run_writes_manifest_without_downloading(tmp_path):
    def http_get(url, params=None):
        return {"_embedded": {"studies": _REST.get(params["efoTrait"], [])}}
    def list_dir(url):        # would be called only on a real fetch
        raise AssertionError("dry-run must not touch the FTP")
    def http_stream(url):
        raise AssertionError("dry-run must not download")

    manifest = pipeline_04b.run(_cfg(tmp_path), http_get, list_dir, http_stream,
                                dry_run=True, limit=None)
    import csv
    rows = list(csv.DictReader(open(manifest)))
    by_acc = {r["accession"]: r for r in rows}
    # burden dropped, both non-burden studies present
    assert by_acc["GCST002"]["disposition"] == "excluded"
    assert by_acc["GCST002"]["reason"] == "burden"
    assert by_acc["GCST001"]["disposition"] == "excluded"   # dry-run: not fetched
    assert by_acc["GCST003"]["disposition"] == "excluded"
    assert {r["reason"] for r in rows if r["accession"] in ("GCST001", "GCST003")} == {"dry_run"}


def test_missing_harmonised_file_is_recorded_not_fatal(tmp_path):
    def http_get(url, params=None):
        return {"_embedded": {"studies": _REST.get(params["efoTrait"], [])}}
    def list_dir(url):
        return []                       # no harmonised file anywhere
    def http_stream(url):
        raise AssertionError("nothing to download")

    manifest = pipeline_04b.run(_cfg(tmp_path), http_get, list_dir, http_stream,
                                dry_run=False, limit=None)
    import csv
    rows = {r["accession"]: r for r in csv.DictReader(open(manifest))}
    assert rows["GCST001"]["reason"] == "no_harmonised_file"
    assert rows["GCST003"]["reason"] == "no_harmonised_file"


def test_successful_fetch_records_fetched_and_variant_count(tmp_path):
    import gzip
    hdr = ("chromosome\tbase_pair_location\teffect_allele\tother_allele\tbeta\t"
           "standard_error\teffect_allele_frequency\tp_value\tvariant_id\t"
           "hm_coordinate_conversion\thm_code\trsid\n")

    def http_get(url, params=None):
        return {"_embedded": {"studies": _REST.get(params["efoTrait"], [])}}
    def list_dir(url):
        return ["file.h.tsv.gz"]
    def http_stream(url):
        body = (hdr + "1\t10177\tAC\tA\t-0.08\t0.20\t0.45\t0.66\tNA\tlo\t11\trs1\n").encode()
        gz = gzip.compress(body)
        return {"Content-Length": str(len(gz))}, [gz]

    manifest = pipeline_04b.run(_cfg(tmp_path), http_get, list_dir, http_stream,
                                dry_run=False, limit=None)
    import csv
    rows = {r["accession"]: r for r in csv.DictReader(open(manifest))}
    assert rows["GCST001"]["disposition"] == "fetched"
    assert rows["GCST001"]["n_variants"] == "1"


def test_cached_raw_is_reharmonized_without_redownloading(tmp_path):
    # A prior run left the raw .h.tsv.gz on disk; re-running must harmonize from
    # cache and never re-download (recovering harmonize_failed after a fixer
    # change costs no ~20 GB re-fetch).
    import gzip
    hdr = ("chromosome\tbase_pair_location\teffect_allele\tother_allele\tbeta\t"
           "standard_error\teffect_allele_frequency\tp_value\tvariant_id\t"
           "hm_coordinate_conversion\thm_code\trsid\n")
    raw_dir = tmp_path / "raw"
    raw_dir.mkdir(parents=True)
    body = (hdr + "1\t10177\tAC\tA\t-0.08\t0.20\t0.45\t0.66\tNA\tlo\t11\trs1\n").encode()
    for acc in ("GCST001", "GCST003"):
        with gzip.open(raw_dir / f"{acc}.h.tsv.gz", "wb") as f:
            f.write(body)
        # B5: a cached raw is only trusted when its verification sidecar exists
        # (an interrupted stream leaves a short file that still decompresses).
        pipeline_04b.write_verification(raw_dir / f"{acc}.h.tsv.gz")

    def http_get(url, params=None):
        return {"_embedded": {"studies": _REST.get(params["efoTrait"], [])}}
    def list_dir(url):
        return ["f.h.tsv.gz"]
    def http_stream(url):
        raise AssertionError("cached raw must not be re-downloaded")

    manifest = pipeline_04b.run(
        _cfg(tmp_path), http_get, list_dir, http_stream, dry_run=False, limit=None)
    import csv
    rows = {r["accession"]: r for r in csv.DictReader(open(manifest))}
    assert rows["GCST001"]["disposition"] == "fetched"
    assert rows["GCST001"]["n_variants"] == "1"
    assert rows["GCST003"]["disposition"] == "fetched"


def test_download_failure_is_recorded_not_fatal(tmp_path):
    # oversized Content-Length -> FileTooLargeError -> reason download_failed, run continues
    def http_get(url, params=None):
        return {"_embedded": {"studies": _REST.get(params["efoTrait"], [])}}
    def list_dir(url):
        return ["f.h.tsv.gz"]
    def http_stream(url):
        return {"Content-Length": str(5_000_000_000)}, [b"x"]   # exceeds max_file_mb cap
    manifest = pipeline_04b.run(
        _cfg(tmp_path), http_get, list_dir, http_stream, dry_run=False, limit=None)
    import csv
    rows = {r["accession"]: r for r in csv.DictReader(open(manifest))}
    assert rows["GCST001"]["reason"] == "download_failed"
    assert rows["GCST003"]["reason"] == "download_failed"   # loop continued to the 2nd kept study


def test_harmonize_failure_is_recorded_not_fatal(tmp_path):
    # a downloaded file with the WRONG columns -> harmonize raises -> reason harmonize_failed,
    # manifest is still written for ALL studies
    import gzip
    def http_get(url, params=None):
        return {"_embedded": {"studies": _REST.get(params["efoTrait"], [])}}
    def list_dir(url):
        return ["f.h.tsv.gz"]
    def http_stream(url):
        bad = gzip.compress(b"wrong\tcolumns\n1\t2\n")   # missing GWAS-SSF columns
        return {"Content-Length": str(len(bad))}, [bad]
    manifest = pipeline_04b.run(
        _cfg(tmp_path), http_get, list_dir, http_stream, dry_run=False, limit=None)
    import csv
    rows = {r["accession"]: r for r in csv.DictReader(open(manifest))}
    assert rows["GCST001"]["reason"] == "harmonize_failed"
    assert rows["GCST003"]["reason"] == "harmonize_failed"
    assert Path(manifest).exists()   # a bad file did NOT destroy the whole manifest


def test_limit_skipped_studies_marked_not_attempted(tmp_path):
    import gzip
    hdr = ("chromosome\tbase_pair_location\teffect_allele\tother_allele\tbeta\tstandard_error\t"
           "effect_allele_frequency\tp_value\tvariant_id\thm_coordinate_conversion\thm_code\trsid\n")
    def http_get(url, params=None):
        return {"_embedded": {"studies": _REST.get(params["efoTrait"], [])}}
    def list_dir(url):
        return ["f.h.tsv.gz"]
    def http_stream(url):
        body = (hdr + "1\t10177\tAC\tA\t-0.08\t0.20\t0.45\t0.66\tNA\tlo\t11\trs1\n").encode()
        gz = gzip.compress(body)
        return {"Content-Length": str(len(gz))}, [gz]
    manifest = pipeline_04b.run(
        _cfg(tmp_path), http_get, list_dir, http_stream, dry_run=False, limit=1)
    import csv
    rows = {r["accession"]: r for r in csv.DictReader(open(manifest))}
    fetched = [a for a, r in rows.items() if r["disposition"] == "fetched"]
    not_attempted = [a for a, r in rows.items() if r["reason"] == "not_attempted"]
    assert len(fetched) == 1                     # only the first kept study fetched
    assert len(not_attempted) == 1               # the other kept study marked not_attempted
    # every resolved study still present exactly once
    assert set(rows) == {"GCST001", "GCST002", "GCST003"}


# --- appended: integrity fixes (B8 dry-run isolation, B5 cache trust) -------

def _good_stream(hdr):
    import gzip

    def http_stream(url):
        body = (hdr + "1\t10177\tAC\tA\t-0.08\t0.20\t0.45\t0.66\tNA\tlo\t11\trs1\n").encode()
        gz = gzip.compress(body)
        return {"Content-Length": str(len(gz))}, [gz]
    return http_stream


_HDR = ("chromosome\tbase_pair_location\teffect_allele\tother_allele\tbeta\t"
        "standard_error\teffect_allele_frequency\tp_value\tvariant_id\t"
        "hm_coordinate_conversion\thm_code\trsid\n")


def test_dry_run_never_overwrites_the_real_manifest(tmp_path):
    """B8: --dry-run rewrote the PRISMA manifest with every study 'excluded'."""
    cfg = _cfg(tmp_path)
    real = Path(cfg["output"]["manifest"])
    real.parent.mkdir(parents=True, exist_ok=True)
    real.write_text("accession,disposition\nGCST999,fetched\n")

    def http_get(url, params=None):
        return {"_embedded": {"studies": _REST.get(params["efoTrait"], [])}}

    written = pipeline_04b.run(cfg, http_get, lambda u: [], _good_stream(_HDR),
                               dry_run=True, limit=None)

    assert Path(written) != real
    assert real.read_text() == "accession,disposition\nGCST999,fetched\n"


def test_a_cached_raw_without_verification_is_redownloaded(tmp_path):
    """B5: a truncated gzip left by an interrupted download must not be trusted."""
    import gzip

    raw_dir = tmp_path / "raw"
    raw_dir.mkdir(parents=True)
    body = (_HDR + "1\t10177\tAC\tA\t-0.08\t0.20\t0.45\t0.66\tNA\tlo\t11\trs1\n").encode()
    complete = gzip.compress(body)
    for acc in ("GCST001", "GCST003"):
        (raw_dir / f"{acc}.h.tsv.gz").write_bytes(complete[:-6])   # truncated, no .ok

    calls = {"n": 0}

    def http_stream(url):
        calls["n"] += 1
        gz = gzip.compress(body)
        return {"Content-Length": str(len(gz))}, [gz]

    def http_get(url, params=None):
        return {"_embedded": {"studies": _REST.get(params["efoTrait"], [])}}

    manifest = pipeline_04b.run(
        _cfg(tmp_path), http_get, lambda u: ["f.h.tsv.gz"], http_stream,
        dry_run=False, limit=None)

    import csv
    rows = {r["accession"]: r for r in csv.DictReader(open(manifest))}
    assert rows["GCST001"]["disposition"] == "fetched"
    assert calls["n"] == 2, "both truncated files must be re-downloaded"
    assert pipeline_04b.raw_is_verified(raw_dir / "GCST001.h.tsv.gz")


#!/usr/bin/env bash
# Fetch the external reference data the genomics arm needs, with checksums.
#
# These files are the second reproducibility gap: MAGMA's 1000 Genomes EUR LD
# panel and NCBI gene locations are ~4 GB, live under a gitignored tools/
# directory, and are captured by neither renv.lock nor uv.lock. A fresh clone
# had the code and no reproducible way to obtain what it needs. This script is
# that way, and it verifies what it downloads: the checksums below are measured
# from the copies that produced the published results, so a changed upstream
# archive fails loudly instead of silently altering the gene-based analysis.
#
# The panel is deliberately not baked into the container image -- 4 GB would
# make it undistributable -- so this populates a directory the image mounts.
#
# Usage:
#   containers/fetch_reference_data.sh [target-dir]      # default: tools/magma
#   containers/fetch_reference_data.sh --verify-only     # check what is there

set -euo pipefail

TARGET="${1:-tools/magma}"
VERIFY_ONLY=0
if [[ "${1:-}" == "--verify-only" ]]; then
    VERIFY_ONLY=1
    TARGET="tools/magma"
fi

# sha256 measured 2026-08-15 on the artefacts behind results/meta/genomics/magma.
G1000_SHA256="83a48fd9dcaa0b9a874b18c63143a4ede93f05505b215b0bd8790130a0d7a954"
NCBI_SHA256="8c2bc9cba581819b2df5a63e741a360d92efef7d819f553a06e2d946ae60a90b"

# MAGMA distributes through SURFsara; the auxiliary-file URLs are versioned by
# an opaque token, so they are recorded here rather than constructed.
G1000_URL="https://vu.data.surfsara.nl/index.php/s/VZNByNwpD8qqINe/download"
NCBI_URL="https://vu.data.surfsara.nl/index.php/s/Pj2orwuVBhYCZWX/download"

log() { printf '%s\n' "$*" >&2; }

verify() {
    local file="$1" want="$2"
    [[ -f "$file" ]] || { log "  MISSING  $file"; return 1; }
    local got
    got="$(sha256sum "$file" | cut -d' ' -f1)"
    if [[ "$got" == "$want" ]]; then
        log "  ok       $file"
        return 0
    fi
    log "  MISMATCH $file"
    log "           expected $want"
    log "           got      $got"
    return 1
}

fetch() {
    local url="$1" out="$2" want="$3"
    if [[ -f "$out" ]] && verify "$out" "$want" >/dev/null 2>&1; then
        log "  cached   $out"
        return 0
    fi
    log "  fetching $out"
    curl -fL --retry 3 --retry-delay 5 "$url" -o "$out.part"
    mv "$out.part" "$out"
    verify "$out" "$want"
}

if [[ "$VERIFY_ONLY" == "1" ]]; then
    log "Verifying reference data in $TARGET"
    rc=0
    verify "$TARGET/g1000_eur.zip" "$G1000_SHA256" || rc=1
    verify "$TARGET/NCBI37.3.zip" "$NCBI_SHA256" || rc=1
    for f in g1000_eur.bed g1000_eur.bim g1000_eur.fam NCBI37.3.gene.loc; do
        [[ -f "$TARGET/$f" ]] && log "  ok       $TARGET/$f" || { log "  MISSING  $TARGET/$f"; rc=1; }
    done
    exit "$rc"
fi

mkdir -p "$TARGET"
log "Fetching MAGMA reference data into $TARGET (~4 GB unpacked)"

fetch "$G1000_URL" "$TARGET/g1000_eur.zip" "$G1000_SHA256"
fetch "$NCBI_URL" "$TARGET/NCBI37.3.zip" "$NCBI_SHA256"

log "Unpacking"
unzip -o -q "$TARGET/g1000_eur.zip" -d "$TARGET"
unzip -o -q "$TARGET/NCBI37.3.zip" -d "$TARGET"

log "Done. Verify at any time with: containers/fetch_reference_data.sh --verify-only"

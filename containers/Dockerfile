# Analysis environment for CP_multi-omics.
#
# Closes the two reproducibility gaps recorded in
# plan/2026-08-15-reproducibility-check.md:
#
#   1. The R environment is described but absent. renv.lock pins R 4.5.3 and
#      Bioconductor 3.22; the workstation runs R 4.3.3 with 56 of the 145
#      locked packages missing, so `renv::restore()` cannot bridge it and
#      05b/05 cannot run. The base image supplies exactly that pinned pair.
#
#   2. MAGMA is present but undescribed. The binary and the 1000G EUR panel
#      live under a gitignored tools/ directory, captured by neither lockfile,
#      so a fresh clone has no reproducible way to obtain them. The binary is
#      baked in here against a recorded checksum; the 4 GB panel is fetched to
#      a mounted volume by containers/fetch_reference_data.sh, which verifies
#      its checksum too.
#
# 05b's failure mode is why this image installs annotation packages explicitly
# rather than trusting the lockfile alone: with one absent, 05b does not error,
# it silently keeps probe ids instead of gene symbols and every downstream
# number moves for the wrong reason.
#
# Build:  docker build -f containers/Dockerfile -t cp-multiomics:latest .
# Verify: docker run --rm cp-multiomics:latest Rscript containers/verify_environment.R

FROM bioconductor/bioconductor_docker:RELEASE_3_22

LABEL org.opencontainers.image.title="CP_multi-omics analysis environment" \
      org.opencontainers.image.description="R 4.5 / Bioconductor 3.22 + Python 3.12 + MAGMA 1.10" \
      org.opencontainers.image.licenses="MIT"

ENV DEBIAN_FRONTEND=noninteractive \
    RENV_CONFIG_REPOS_OVERRIDE="" \
    RENV_PATHS_CACHE=/renv-cache

# unzip and curl are needed by fetch_reference_data.sh; libcurl/libxml2 back
# biomaRt and the Bioconductor annotation packages.
RUN apt-get update && apt-get install -y --no-install-recommends \
        curl unzip ca-certificates git \
        libcurl4-openssl-dev libxml2-dev libssl-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /work

# ---------------------------------------------------------------------------
# R packages
# ---------------------------------------------------------------------------
# Restore first from the lockfile, which is the authority for versions.
COPY renv.lock /work/renv.lock
RUN Rscript -e 'install.packages("renv", repos="https://cloud.r-project.org")' \
 && Rscript -e 'renv::restore(lockfile="/work/renv.lock", prompt=FALSE)'

# Then install every annotation package GPL_TO_ANNO_PKG can select. The lockfile
# held only 5 of these 11 until 2026-08-15 because renv discovers dependencies
# by static analysis and 05b names them through a lookup table. _dependencies.R
# now lists all 11 so future snapshots capture them, but this layer does not
# depend on that having happened -- a missing annotation package is a silent
# wrong answer, not a build failure, so it is pinned down twice.
RUN Rscript -e 'BiocManager::install(c( \
      "rgu34a.db", "rgu34b.db", "rgu34c.db", "rae230a.db", "rat2302.db", \
      "mouse4302.db", "mogene10sttranscriptcluster.db", \
      "mogene11sttranscriptcluster.db", "mta10transcriptcluster.db", \
      "illuminaMousev1.db", "illuminaMousev2.db", "illuminaRatv1.db", \
      "org.Mm.eg.db", "org.Rn.eg.db", \
      "DESeq2", "edgeR", "limma", "biomaRt"), ask=FALSE, update=FALSE)'

# ---------------------------------------------------------------------------
# Python
# ---------------------------------------------------------------------------
COPY --from=ghcr.io/astral-sh/uv:0.5.11 /uv /usr/local/bin/uv
COPY pyproject.toml uv.lock /work/
RUN uv sync --frozen --extra dev --python 3.12

# ---------------------------------------------------------------------------
# MAGMA
# ---------------------------------------------------------------------------
# v1.10 linux, the version recorded in results/meta/genomics/magma/*.log. The
# checksum is the one measured on the copy that produced the published results,
# so a substituted upstream archive fails the build rather than silently
# changing the gene-based analysis.
ARG MAGMA_URL=https://vu.data.surfsara.nl/index.php/s/zkKbNeNOZAhFXZB/download
ARG MAGMA_SHA256=cfaf44469064d3b4b6717cea11666b9754d0250950d35ce290d059327ae5a583
RUN set -eux; \
    mkdir -p /opt/magma; \
    curl -fsSL "$MAGMA_URL" -o /tmp/magma.zip; \
    unzip -o -j /tmp/magma.zip magma -d /opt/magma || unzip -o /tmp/magma.zip -d /opt/magma; \
    echo "${MAGMA_SHA256}  /opt/magma/magma" | sha256sum -c -; \
    chmod +x /opt/magma/magma; \
    rm -f /tmp/magma.zip
ENV PATH="/opt/magma:${PATH}"

COPY . /work

# The repository ships an renv .Rprofile, and copying it in makes every R
# process inside the container autoload renv, which repoints .libPaths() at
# /work/renv/library -- a directory .dockerignore deliberately excludes. The
# result is an image with a complete environment that reports "no package
# called 'jsonlite'", because renv is hiding the library the build populated.
#
# The autoloader's job is to recreate a pinned environment; here the image *is*
# that environment, restored from the same lockfile at build time. So it is
# switched off, and R uses the site library where everything was installed.
# Set RENV_CONFIG_AUTOLOADER_ENABLED=TRUE to opt back in if you mount a
# populated renv/library from the host.
ENV RENV_CONFIG_AUTOLOADER_ENABLED=FALSE

# Reference panels are mounted, not baked: g1000_eur is ~4 GB and would make the
# image unusable to distribute. fetch_reference_data.sh populates /work/tools.
VOLUME ["/work/tools", "/work/data", "/work/results"]

CMD ["bash"]

# Run this once to initialize renv and install core packages.
# After: commit renv.lock and renv/activate.R to version control.

renv::init()

packages <- c(
  # Meta-analysis
  "metafor", "meta",
  # CLI / config / IO
  "optparse", "yaml", "jsonlite",
  # Data wrangling
  "dplyr", "tidyr",
  # HTTP (genomics_da.R)
  "httr",
  # Visualization (step 06 figures)
  "ggplot2",
  # Bioconductor dependency manager
  "BiocManager"
)

install.packages(packages)

BiocManager::install(c(
  "GEOquery",       # step 05: download GEO data
  "DESeq2",         # step 05: RNA-seq DA
  "limma",          # step 05: microarray / proteomics DA
  "edgeR",          # step 05: RNA-seq (alternative normalisation)
  "clusterProfiler",# step 07+: pathway enrichment
  "org.Hs.eg.db",   # human gene annotation
  "org.Mm.eg.db",   # mouse gene annotation
  "AnnotationDbi"
))

# Snapshot after installation
renv::snapshot()

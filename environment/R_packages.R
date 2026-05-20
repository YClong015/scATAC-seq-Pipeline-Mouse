# ============================================================
# R package installation script for scatac-aging-disease
# Tested with R 4.3.1 / 4.4.2
# Run: Rscript environment/R_packages.R
# ============================================================

# ── CRAN packages ──────────────────────────────────────────────
cran_pkgs <- c(
  # Core single-cell
  "Seurat",            # >=5.0
  "harmony",
  "dplyr", "tidyr", "ggplot2", "patchwork", "ggrepel", "scales",
  "Matrix", "data.table",
  "future",
  # Plotting
  "RColorBrewer", "ggforce",
  "UpSetR",
  "magick",   # PNG composition for assemble_Fig6.R / Fig7.R
  "pdftools", # PDF composition
  # Aim 3 helpers
  "tibble",
  # Misc
  "BiocManager"
)

new_cran <- cran_pkgs[!sapply(cran_pkgs, requireNamespace, quietly = TRUE)]
if (length(new_cran) > 0) {
  install.packages(new_cran, repos = "https://cloud.r-project.org", quiet = FALSE)
}

# ── Bioconductor packages ──────────────────────────────────────
bioc_pkgs <- c(
  "Signac",                       # >=1.10
  "DESeq2",                       # >=1.40
  "tximport",                     # only used by sugarcane-wgcna side project
  "ComplexHeatmap", "circlize",
  "EnsDb.Mmusculus.v79",
  "BSgenome.Mmusculus.UCSC.mm10",
  "GenomicRanges", "GenomeInfoDb", "rtracklayer",
  "BiocParallel",
  "scDblFinder",
  # ComplexUpset alternative to UpSetR (for Fig6 v2)
  "ComplexUpset"
)

new_bioc <- bioc_pkgs[!sapply(bioc_pkgs, requireNamespace, quietly = TRUE)]
if (length(new_bioc) > 0) {
  BiocManager::install(new_bioc, update = FALSE, ask = FALSE)
}

# ── Optional: motif tools (only needed by deprecated DAR_general.R) ──
#
# install.packages("https://cran.r-project.org/src/contrib/Archive/TFMPvalue/TFMPvalue_0.0.9.tar.gz",
#                  repos = NULL, type = "source")
# BiocManager::install(c("TFBSTools", "motifmatchr", "JASPAR2022"), update = FALSE)

cat("\n=== Installed ===\n")
installed.packages()[c(cran_pkgs, bioc_pkgs), c("Package", "Version")]

# R package installation for scATAC-ageing-disease-mouse (tested with R 4.4.2).
# Run: Rscript environment/R_packages.R

### CRAN packages ###
cran_pkgs <- c(
  "Seurat",
  "harmony",
  "dplyr", "tidyr", "ggplot2", "patchwork", "ggrepel", "scales",
  "Matrix", "data.table",
  "future",
  "RColorBrewer", "ggforce",
  "UpSetR",
  "magick",   
  "pdftools", 
  "tibble",
  "BiocManager"
)

new_cran <- cran_pkgs[!sapply(cran_pkgs, requireNamespace, quietly = TRUE)]
if (length(new_cran) > 0) {
  install.packages(new_cran, repos = "https://cloud.r-project.org", quiet = FALSE)
}

### Bioconductor packages ###
bioc_pkgs <- c(
  "Signac",                       
  "DESeq2",                       
  "tximport",                    
  "ComplexHeatmap", "circlize",
  "EnsDb.Mmusculus.v79",
  "BSgenome.Mmusculus.UCSC.mm10",
  "GenomicRanges", "GenomeInfoDb", "rtracklayer",
  "BiocParallel",
  "scDblFinder",
  "ComplexUpset"
)

new_bioc <- bioc_pkgs[!sapply(bioc_pkgs, requireNamespace, quietly = TRUE)]
if (length(new_bioc) > 0) {
  BiocManager::install(new_bioc, update = FALSE, ask = FALSE)
}

cat("\n=== Installed ===\n")
installed.packages()[c(cran_pkgs, bioc_pkgs), c("Package", "Version")]

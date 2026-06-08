#!/usr/bin/env Rscript
# Split T-cell fragments by the final (Masopust) annotation for per-cell-type peak counts.
# Input tcells_final_annotated.rds; output split BEDs + DataList.txt.

suppressPackageStartupMessages({
  library(Signac)
  library(Seurat)
})

BASE_DIR     <- "/QRISdata/Q8448/Mouse_disease_data/Tcells"
IN_RDS       <- file.path(BASE_DIR, "tcells_final_annotated.rds")
FRAG_FILE    <- file.path(BASE_DIR, "atac_fragments.tsv.gz")
OUT_BED_DIR  <- file.path(BASE_DIR, "fragment_files_split_by_celltype_final")
GROUP_COL    <- "cell_type_final"

stopifnot(file.exists(IN_RDS), file.exists(FRAG_FILE))
dir.create(OUT_BED_DIR, showWarnings = FALSE, recursive = TRUE)

message("Loading ", IN_RDS)
tc <- readRDS(IN_RDS)
message("  cells = ", ncol(tc))

if (!GROUP_COL %in% colnames(tc@meta.data))
  stop("`", GROUP_COL, "` missing. Re-run Tcell_final_annotate.R first.")
message("Group column: ", GROUP_COL)
print(table(tc[[GROUP_COL]][[1]], useNA = "ifany"))

# Re-link the fragment file: cleaned objects often carry a stale absolute path
# baked in from the original Signac build.
DefaultAssay(tc) <- "ATAC"
new_frag <- CreateFragmentObject(path = FRAG_FILE, cells = colnames(tc))
Fragments(tc) <- NULL
Fragments(tc) <- new_frag
message("Fragment file relinked: ", FRAG_FILE)

# Drop cells with NA in the grouping column (SplitFragments errors otherwise)
keep <- !is.na(tc[[GROUP_COL]][[1]])
if (any(!keep)) {
  message(sprintf("Dropping %d cells with NA %s before split.",
                  sum(!keep), GROUP_COL))
  tc <- subset(tc, cells = colnames(tc)[keep])
}

message("\nSplitting fragments by ", GROUP_COL, " into:\n  ", OUT_BED_DIR)
SplitFragments(
  object   = tc,
  assay    = "ATAC",
  group.by = GROUP_COL,
  outdir   = OUT_BED_DIR,
  verbose  = TRUE
)

# Build DataList.txt for the SLURM array (one BED filename per line, no path)
bed_files <- list.files(OUT_BED_DIR, pattern = "\\.bed$", full.names = FALSE)
if (length(bed_files) == 0)
  stop("SplitFragments wrote no BED files. Check the run log above.")
writeLines(sort(bed_files), file.path(OUT_BED_DIR, "DataList.txt"))

message("\nWrote ", length(bed_files), " BED files:")
for (b in sort(bed_files)) message("  ", b)
message("\nDataList.txt: ", file.path(OUT_BED_DIR, "DataList.txt"))

#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --job-name=Lung_BED_HOMER
#SBATCH --time=12:00:00
#SBATCH --partition=general
#SBATCH --account=a_imb_ccbcd
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

set -euo pipefail

# ============================================================
# Step 1: Generate stable/NS BED files from DESeq2_all.tsv
# Step 2: Run HOMER for Lung (stable_bg and NS_bg)
# ============================================================

DAR_DIR="/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_pseudobulk_DESeq2/DAR_pseudobulk_Lung_DESeq2/Lung_DAR_tmp"
HOMER_HOME="/scratch/user/s4869245/homer"
GENOME="mm10"
CPU="${SLURM_CPUS_PER_TASK:-8}"
MIN_PEAKS=20

HOMER_STABLE="${DAR_DIR}/HOMER_stable_bg"
HOMER_NS="${DAR_DIR}/HOMER_NS_bg"
BED_DIR="${DAR_DIR}/DAR_BED_HOMER"

mkdir -p "${BED_DIR}" "${HOMER_STABLE}/logs" "${HOMER_NS}/logs"

echo "$(date)  Generating BED files..."

# ── Step 1: BED generation from _DESeq2_all.tsv ────────────
module load r/4.4.2

Rscript - <<'REOF'
library(dplyr)

dar_dir  <- Sys.getenv("DAR_DIR",  "/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_pseudobulk_DESeq2/DAR_pseudobulk_Lung_DESeq2/Lung_DAR_tmp")
bed_dir  <- file.path(dar_dir, "DAR_BED_HOMER")
contrast <- "Case_vs_Control"

# thresholds
padj_cut        <- 0.05
stable_padj_cut <- 0.9
stable_lfc_cut  <- 0.05

tsv_files <- list.files(
  file.path(dar_dir, "DAR_tables"),
  pattern = "_DESeq2_all\\.tsv$",
  full.names = TRUE
)
message("Found ", length(tsv_files), " DESeq2_all.tsv files")

peaks_to_bed <- function(peaks, out_file) {
  if (length(peaks) == 0) { message("  (0 peaks, skipping)"); return(invisible(NULL)) }
  df <- data.frame(
    Chr   = sub("(.*)-.*-.*",  "\\1", peaks),
    Start = sub(".*-(.*)-.*",  "\\1", peaks),
    End   = sub(".*-.*-(.*)",  "\\1", peaks),
    stringsAsFactors = FALSE
  )
  write.table(df, file = out_file, sep = "\t",
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  message("  Wrote ", nrow(df), " peaks -> ", basename(out_file))
}

for (f in tsv_files) {
  # parse cell type from filename: {ct}__{contrast}__{thresh}_DESeq2_all.tsv
  bname <- sub("_DESeq2_all\\.tsv$", "", basename(f))
  parts <- strsplit(bname, "__", fixed = TRUE)[[1]]
  if (length(parts) < 2) next
  ct       <- parts[1]
  contrast_f <- parts[2]
  if (contrast_f != contrast) next

  message("\n--- ", ct, " ---")
  res <- tryCatch(
    read.delim(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  if (is.null(res) || nrow(res) == 0) { message("  Empty, skip"); next }
  if (!all(c("peak","log2FoldChange","padj") %in% colnames(res))) {
    message("  Missing columns, skip"); next
  }
  res <- res[!is.na(res$padj), ]

  opening <- res[res$padj <  padj_cut & res$log2FoldChange >  0, ]
  closing  <- res[res$padj <  padj_cut & res$log2FoldChange <  0, ]
  stable   <- res[res$padj >  stable_padj_cut & abs(res$log2FoldChange) < stable_lfc_cut, ]
  ns       <- res[res$padj >  padj_cut, ]

  message("  opening=", nrow(opening), " closing=", nrow(closing),
          " stable=", nrow(stable), " NS=", nrow(ns))

  pfx <- file.path(bed_dir, paste0(ct, "__", contrast))
  peaks_to_bed(opening$peak, paste0(pfx, "__opening.bed"))
  peaks_to_bed(closing$peak, paste0(pfx, "__closing.bed"))
  peaks_to_bed(stable$peak,  paste0(pfx, "__stable.bed"))
  peaks_to_bed(ns$peak,      paste0(pfx, "__NS.bed"))
}
message("\nBED generation done.")
REOF

echo "$(date)  BED generation complete."

# ── Step 2: HOMER ──────────────────────────────────────────
export PATH="${HOMER_HOME}/bin:${PATH}"

CONTRAST="Case_vs_Control"

run_homer() {
  local fg="$1"
  local bg="$2"
  local out="$3"
  local label="$4"

  if [ ! -f "${fg}" ]; then
    echo "  SKIP (no foreground): ${label}"; return
  fi
  if [ ! -f "${bg}" ]; then
    echo "  SKIP (no background): ${label}"; return
  fi

  local fg_n=$(wc -l < "${fg}")
  local bg_n=$(wc -l < "${bg}")

  if [ "${fg_n}" -lt "${MIN_PEAKS}" ]; then
    echo "  SKIP (fg too few peaks=${fg_n}): ${label}"; return
  fi
  if [ "${bg_n}" -lt "${MIN_PEAKS}" ]; then
    echo "  SKIP (bg too few peaks=${bg_n}): ${label}"; return
  fi

  if [ -f "${out}/knownResults.txt" ]; then
    echo "  SKIP (already done): ${label}"; return
  fi

  mkdir -p "${out}"
  echo "  Running HOMER: ${label}  fg=${fg_n} bg=${bg_n}"
  "${HOMER_HOME}/bin/findMotifsGenome.pl" "${fg}" "${GENOME}" "${out}" \
    -bg "${bg}" -size given -p "${CPU}" -nomotif \
    > "${out}/homer.log" 2>&1 \
    && echo "  Done: ${label}" \
    || echo "  ERROR: ${label} (see ${out}/homer.log)"
}

echo "$(date)  Starting HOMER..."

for opening_bed in "${BED_DIR}"/*__${CONTRAST}__opening.bed; do
  [ -f "${opening_bed}" ] || continue
  base=$(basename "${opening_bed}" "__${CONTRAST}__opening.bed")
  ct="${base}"

  pfx="${BED_DIR}/${ct}__${CONTRAST}"
  opening="${pfx}__opening.bed"
  closing="${pfx}__closing.bed"
  stable="${pfx}__stable.bed"
  ns="${pfx}__NS.bed"

  thresh="005"
  s_pfx="${HOMER_STABLE}/${ct}__${CONTRAST}__${thresh}"
  n_pfx="${HOMER_NS}/${ct}__${CONTRAST}__${thresh}"

  echo ""
  echo "=== ${ct} ==="

  # stable background: 4 comparisons
  run_homer "${opening}" "${stable}" "${s_pfx}__opening_vs_stable" "${ct} opening_vs_stable"
  run_homer "${closing}"  "${stable}" "${s_pfx}__closing_vs_stable"  "${ct} closing_vs_stable"
  run_homer "${stable}"  "${opening}" "${s_pfx}__stable_vs_opening" "${ct} stable_vs_opening"
  run_homer "${stable}"  "${closing}"  "${s_pfx}__stable_vs_closing"  "${ct} stable_vs_closing"

  # NS background: 4 comparisons
  run_homer "${opening}" "${ns}"     "${n_pfx}__opening_vs_NS"     "${ct} opening_vs_NS"
  run_homer "${closing}"  "${ns}"     "${n_pfx}__closing_vs_NS"      "${ct} closing_vs_NS"
  run_homer "${ns}"      "${opening}" "${n_pfx}__NS_vs_opening"     "${ct} NS_vs_opening"
  run_homer "${ns}"      "${closing}"  "${n_pfx}__NS_vs_closing"      "${ct} NS_vs_closing"

done

echo ""
echo "$(date)  All done."

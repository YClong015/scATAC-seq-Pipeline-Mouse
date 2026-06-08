#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --job-name=Lung_HOMER
#SBATCH --time=6:00:00
#SBATCH --partition=general
#SBATCH --account=a_imb_ccbcd
#SBATCH --array=0-13
#SBATCH -o logs_slurm/%x_%A_%a.out
#SBATCH -e logs_slurm/%x_%A_%a.err

set -euo pipefail

## Lung full HOMER pipeline, one SLURM array task per cell type.
## Per task: generate opening/closing/stable/NS BEDs, then 10 HOMER runs
##   stable bg (4) + NS bg (4) + closing-vs-opening (2, identity TFs).

DAR_DIR="/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_pseudobulk_DESeq2/DAR_pseudobulk_Lung_DESeq2"
HOMER_HOME="/scratch/user/s4869245/homer"
GENOME="mm10"
CPU="${SLURM_CPUS_PER_TASK:-8}"
MIN_PEAKS=200
CONTRAST="Case_vs_Control"
CTS=(AT2 B Ciliated EC-vasc Eosinophils Fib Mac-alv Mac-inter Mesothelial Mo-Ly6c+ NK Pen SMCs T)

## HOMER's findMotifsGenome.pl calls child binaries by bare name, so its bin/
## must be on PATH or it silently fails and writes no knownResults.txt.
export PATH="${HOMER_HOME}/bin:${PATH}"

BED_DIR="${DAR_DIR}/DAR_BED_HOMER"
HOMER_STABLE="${DAR_DIR}/HOMER_stable_bg"
HOMER_NS="${DAR_DIR}/HOMER_NS_bg"
HOMER_CVO="${DAR_DIR}/HOMER_closing_vs_opening"
mkdir -p "${BED_DIR}" "${HOMER_STABLE}/logs" "${HOMER_NS}/logs" "${HOMER_CVO}/logs" logs_slurm

CT="${CTS[${SLURM_ARRAY_TASK_ID}]}"
echo "$(date)  Array task ${SLURM_ARRAY_TASK_ID}/${#CTS[@]}  CT=${CT}"

opening="${BED_DIR}/${CT}__${CONTRAST}__opening.bed"
closing="${BED_DIR}/${CT}__${CONTRAST}__closing.bed"
stable="${BED_DIR}/${CT}__${CONTRAST}__stable.bed"
ns="${BED_DIR}/${CT}__${CONTRAST}__NS.bed"

## Step 1: BED generation (only if missing for this CT)
if [ ! -f "${opening}" ] || [ ! -f "${closing}" ] || [ ! -f "${stable}" ] || [ ! -f "${ns}" ]; then
  echo "$(date)  Generating BEDs for ${CT}..."
  module load r/4.4.2
  export DAR_DIR CT CONTRAST
  Rscript - <<'REOF'
suppressPackageStartupMessages(library(dplyr))

dar_dir  <- Sys.getenv("DAR_DIR")
ct       <- Sys.getenv("CT")
contrast <- Sys.getenv("CONTRAST")
bed_dir  <- file.path(dar_dir, "DAR_BED_HOMER")

padj_cut        <- 0.05
stable_padj_cut <- 0.9
stable_lfc_cut  <- 0.05

f <- list.files(file.path(dar_dir, "DAR_tables"),
                pattern = paste0("^", ct, "__", contrast, "__.*_DESeq2_all\\.tsv$"),
                full.names = TRUE)
if (length(f) == 0) stop("No DESeq2 TSV found for ", ct, " / ", contrast)
res <- read.delim(f[1], stringsAsFactors = FALSE)
res <- res[!is.na(res$padj), ]

peaks_to_bed <- function(peaks, out) {
  if (length(peaks) == 0) { message("  (0 peaks, skip ", basename(out), ")"); return() }
  df <- data.frame(
    Chr   = sub("(.*)-.*-.*", "\\1", peaks),
    Start = sub(".*-(.*)-.*", "\\1", peaks),
    End   = sub(".*-.*-(.*)", "\\1", peaks)
  )
  write.table(df, out, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  message("  Wrote ", nrow(df), " peaks -> ", basename(out))
}

pfx <- file.path(bed_dir, paste0(ct, "__", contrast))
peaks_to_bed(res$peak[res$padj <  padj_cut & res$log2FoldChange >  0], paste0(pfx, "__opening.bed"))
peaks_to_bed(res$peak[res$padj <  padj_cut & res$log2FoldChange <  0], paste0(pfx, "__closing.bed"))
peaks_to_bed(res$peak[res$padj >  stable_padj_cut & abs(res$log2FoldChange) < stable_lfc_cut],
             paste0(pfx, "__stable.bed"))
peaks_to_bed(res$peak[res$padj >  padj_cut], paste0(pfx, "__NS.bed"))
REOF
else
  echo "$(date)  BEDs already exist for ${CT}, skipping generation."
fi

## HOMER runner: foreground vs background, skip if too few peaks or already done
run_homer() {
  local fg="$1" bg="$2" out="$3" label="$4"
  [ -f "${fg}" ] || { echo "  SKIP (no fg): ${label}"; return; }
  [ -f "${bg}" ] || { echo "  SKIP (no bg): ${label}"; return; }
  local fg_n bg_n
  fg_n=$(wc -l < "${fg}")
  bg_n=$(wc -l < "${bg}")
  if [ "${fg_n}" -lt "${MIN_PEAKS}" ]; then echo "  SKIP (fg n=${fg_n}<${MIN_PEAKS}): ${label}"; return; fi
  if [ "${bg_n}" -lt "${MIN_PEAKS}" ]; then echo "  SKIP (bg n=${bg_n}<${MIN_PEAKS}): ${label}"; return; fi
  if [ -f "${out}/knownResults.txt" ]; then echo "  SKIP (already done): ${label}"; return; fi
  mkdir -p "${out}"
  echo "  RUN HOMER: ${label}  fg=${fg_n} bg=${bg_n}"
  "${HOMER_HOME}/bin/findMotifsGenome.pl" "${fg}" "${GENOME}" "${out}" \
    -bg "${bg}" -size given -p "${CPU}" -nomotif -mask \
    > "${out}/homer.log" 2>&1 \
    && echo "  Done: ${label}" \
    || echo "  ERROR: ${label} (see ${out}/homer.log)"
}

s_pfx="${HOMER_STABLE}/${CT}__${CONTRAST}__005"
n_pfx="${HOMER_NS}/${CT}__${CONTRAST}__005"
c_pfx="${HOMER_CVO}/${CT}__${CONTRAST}"

echo ""
echo "=== ${CT} / ${CONTRAST}: 10 HOMER comparisons ==="

## stable background (4)
run_homer "${opening}" "${stable}"  "${s_pfx}__opening_vs_stable" "${CT} opening_vs_stable"
run_homer "${closing}" "${stable}"  "${s_pfx}__closing_vs_stable" "${CT} closing_vs_stable"
run_homer "${stable}"  "${opening}" "${s_pfx}__stable_vs_opening" "${CT} stable_vs_opening"
run_homer "${stable}"  "${closing}" "${s_pfx}__stable_vs_closing" "${CT} stable_vs_closing"

## NS background (4)
run_homer "${opening}" "${ns}" "${n_pfx}__opening_vs_NS" "${CT} opening_vs_NS"
run_homer "${closing}" "${ns}" "${n_pfx}__closing_vs_NS" "${CT} closing_vs_NS"
run_homer "${ns}" "${opening}" "${n_pfx}__NS_vs_opening" "${CT} NS_vs_opening"
run_homer "${ns}" "${closing}" "${n_pfx}__NS_vs_closing" "${CT} NS_vs_closing"

## closing-vs-opening reciprocal (2, identity TFs)
run_homer "${closing}" "${opening}" "${c_pfx}__closing_vs_opening" "${CT} closing_vs_opening"
run_homer "${opening}" "${closing}" "${c_pfx}__opening_vs_closing" "${CT} opening_vs_closing"

echo ""
echo "$(date)  Array task ${SLURM_ARRAY_TASK_ID} (${CT}) finished."

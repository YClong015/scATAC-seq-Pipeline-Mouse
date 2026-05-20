#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --job-name=HOMER_close_vs_open
#SBATCH --time=12:00:00
#SBATCH --partition=general
#SBATCH --account=a_imb_ccbcd
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

set -euo pipefail

module load bedtools/2.30.0-gcc-11.3.0
module load r/4.4.2

HOMER_HOME="/scratch/user/s4869245/homer"
GENOME="mm10"
CPU="${SLURM_CPUS_PER_TASK:-8}"
MIN_PEAKS=50
padj_cut=0.05

BASE="/QRISdata/Q8448/Mouse_disease_data/DAR"
TMP_DIR="${BASE}/DAR_closing_vs_opening/tmp_beds"
HOMER_OUT="${BASE}/DAR_closing_vs_opening/HOMER"
OUT_PLOT="${BASE}/DAR_closing_vs_opening/figures"

mkdir -p "${TMP_DIR}" "${HOMER_OUT}" "${OUT_PLOT}"
export PATH="${HOMER_HOME}/bin:${PATH}"

# ── Cell types ────────────────────────────────────────────────
# Format: TISSUE:CELLTYPE:CONTRAST:BASE_DIR:FILE_SUFFIX:FORMAT
# FORMAT: std  = col1=peak, col3=lfc, col7=padj  (Kidney/Lung)
#         aorta= col1=peak, col3=lfc, col7=padj, suffix=exp005
#         tcell= col7=peak, col2=lfc, col6=padj
CT_LIST=(
  "Kidney:DCT:Day42_vs_Sham:/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_pseudobulk_DESeq2/DAR_pseudobulk_Kidney_DESeq2:005_DESeq2_all.tsv:std"
  "Kidney:PC:Day42_vs_Sham:/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_pseudobulk_DESeq2/DAR_pseudobulk_Kidney_DESeq2:005_DESeq2_all.tsv:std"
  "Kidney:PT:Day42_vs_Sham:/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_pseudobulk_DESeq2/DAR_pseudobulk_Kidney_DESeq2:005_DESeq2_all.tsv:std"
  "Kidney:TAL:Day42_vs_Sham:/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_pseudobulk_DESeq2/DAR_pseudobulk_Kidney_DESeq2:005_DESeq2_all.tsv:std"
  "Kidney:Macrophages:Day42_vs_Sham:/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_pseudobulk_DESeq2/DAR_pseudobulk_Kidney_DESeq2:005_DESeq2_all.tsv:std"
  "Lung:AT2:Case_vs_Control:/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_pseudobulk_DESeq2/DAR_pseudobulk_Lung_DESeq2:005_DESeq2_all.tsv:std"
  "Lung:B:Case_vs_Control:/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_pseudobulk_DESeq2/DAR_pseudobulk_Lung_DESeq2:005_DESeq2_all.tsv:std"
  "Lung:Ciliated:Case_vs_Control:/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_pseudobulk_DESeq2/DAR_pseudobulk_Lung_DESeq2:005_DESeq2_all.tsv:std"
  "Lung:EC-vasc:Case_vs_Control:/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_pseudobulk_DESeq2/DAR_pseudobulk_Lung_DESeq2:005_DESeq2_all.tsv:std"
  "Lung:Eosinophils:Case_vs_Control:/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_pseudobulk_DESeq2/DAR_pseudobulk_Lung_DESeq2:005_DESeq2_all.tsv:std"
  "Lung:Fib:Case_vs_Control:/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_pseudobulk_DESeq2/DAR_pseudobulk_Lung_DESeq2:005_DESeq2_all.tsv:std"
  "Lung:Mac-alv:Case_vs_Control:/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_pseudobulk_DESeq2/DAR_pseudobulk_Lung_DESeq2:005_DESeq2_all.tsv:std"
  "Lung:Mac-inter:Case_vs_Control:/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_pseudobulk_DESeq2/DAR_pseudobulk_Lung_DESeq2:005_DESeq2_all.tsv:std"
  "Lung:NK:Case_vs_Control:/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_pseudobulk_DESeq2/DAR_pseudobulk_Lung_DESeq2:005_DESeq2_all.tsv:std"
  "Lung:Pen:Case_vs_Control:/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_pseudobulk_DESeq2/DAR_pseudobulk_Lung_DESeq2:005_DESeq2_all.tsv:std"
  "Lung:SMCs:Case_vs_Control:/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_pseudobulk_DESeq2/DAR_pseudobulk_Lung_DESeq2:005_DESeq2_all.tsv:std"
  "Lung:T:Case_vs_Control:/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_pseudobulk_DESeq2/DAR_pseudobulk_Lung_DESeq2:005_DESeq2_all.tsv:std"
  "Aorta:Macrophages:Challenge_vs_Control:/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_pseudobulk_DESeq2/DAR_pseudobulk_Aorta_DESeq2:exp005_DESeq2_all.tsv:tcell"
  "Aorta:Pericytes:Challenge_vs_Control:/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_pseudobulk_DESeq2/DAR_pseudobulk_Aorta_DESeq2:exp005_DESeq2_all.tsv:tcell"
  "Aorta:SMC:Challenge_vs_Control:/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_pseudobulk_DESeq2/DAR_pseudobulk_Aorta_DESeq2:exp005_DESeq2_all.tsv:tcell"
  "Tcell:Tcell:Young_chronic_vs_Young_control:/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_pseudobulk_DESeq2/DAR_pseudobulk_Tcells_DESeq2:005_DESeq2_all.tsv:tcell"
)

# ── Step 1: Generate BED files ────────────────────────────────
echo "$(date)  Step 1: Generating BED files..."

for entry in "${CT_LIST[@]}"; do
  TISSUE=$(echo "${entry}"   | cut -d: -f1)
  CT=$(echo "${entry}"       | cut -d: -f2)
  CONTRAST=$(echo "${entry}" | cut -d: -f3)
  BASE_DIR=$(echo "${entry}" | cut -d: -f4)
  SUFFIX=$(echo "${entry}"   | cut -d: -f5)
  FORMAT=$(echo "${entry}"   | cut -d: -f6)

  TSV="${BASE_DIR}/DAR_tables/${CT}__${CONTRAST}__${SUFFIX}"
  PREFIX="${TMP_DIR}/${TISSUE}_${CT}"
  LABEL="${TISSUE}_${CT}"

  if [ ! -f "${TSV}" ]; then
    echo "  SKIP (no TSV): ${LABEL} [${TSV}]"; continue
  fi

  if [ "${FORMAT}" = "tcell" ]; then
    # Tcell: col2=lfc, col6=padj, col7=peak(chr-start-end)
    awk -F'\t' -v cut="${padj_cut}" '
      NR>1 && $6~/^[0-9]/ && $6+0 < cut && $2+0 > 0 {
        split($7, a, "-");
        print a[1]"\t"a[2]"\t"a[3]
      }' "${TSV}" | sort -k1,1 -k2,2n > "${PREFIX}_opening.bed"

    awk -F'\t' -v cut="${padj_cut}" '
      NR>1 && $6~/^[0-9]/ && $6+0 < cut && $2+0 < 0 {
        split($7, a, "-");
        print a[1]"\t"a[2]"\t"a[3]
      }' "${TSV}" | sort -k1,1 -k2,2n > "${PREFIX}_closing.bed"
  else
    # std/aorta: col1=peak(chr-start-end), col3=lfc, col7=padj
    awk -F'\t' -v cut="${padj_cut}" '
      NR>1 && $7~/^[0-9]/ && $7+0 < cut && $3+0 > 0 {
        split($1, a, "-");
        print a[1]"\t"a[2]"\t"a[3]
      }' "${TSV}" | sort -k1,1 -k2,2n > "${PREFIX}_opening.bed"

    awk -F'\t' -v cut="${padj_cut}" '
      NR>1 && $7~/^[0-9]/ && $7+0 < cut && $3+0 < 0 {
        split($1, a, "-");
        print a[1]"\t"a[2]"\t"a[3]
      }' "${TSV}" | sort -k1,1 -k2,2n > "${PREFIX}_closing.bed"
  fi

  n_open=$(wc -l < "${PREFIX}_opening.bed")
  n_close=$(wc -l < "${PREFIX}_closing.bed")
  echo "  ${LABEL}: opening=${n_open}  closing=${n_close}"
done

# ── Step 2: Run HOMER (closing vs opening) ────────────────────
echo "$(date)  Step 2: Running HOMER..."

run_homer() {
  local fg="$1" bg="$2" out="$3" label="$4"

  [ -f "${fg}" ] || { echo "  SKIP (no fg): ${label}"; return; }
  [ -f "${bg}" ] || { echo "  SKIP (no bg): ${label}"; return; }

  local fg_n=$(wc -l < "${fg}")
  local bg_n=$(wc -l < "${bg}")

  [ "${fg_n}" -ge "${MIN_PEAKS}" ] || {
    echo "  SKIP (fg=${fg_n}<${MIN_PEAKS}): ${label}"; return; }
  [ "${bg_n}" -ge "${MIN_PEAKS}" ] || {
    echo "  SKIP (bg=${bg_n}<${MIN_PEAKS}): ${label}"; return; }

  if [ -f "${out}/knownResults.txt" ]; then
    echo "  SKIP (done): ${label}"; return
  fi

  mkdir -p "${out}"
  echo "  Running HOMER: ${label}  fg=${fg_n} bg=${bg_n}"
  "${HOMER_HOME}/bin/findMotifsGenome.pl" "${fg}" "${GENOME}" "${out}" \
    -bg "${bg}" -size given -p "${CPU}" -nomotif \
    > "${out}/homer.log" 2>&1 \
    && echo "  Done: ${label}" \
    || echo "  ERROR: ${label}"
}

for entry in "${CT_LIST[@]}"; do
  TISSUE=$(echo "${entry}"  | cut -d: -f1)
  CT=$(echo "${entry}"      | cut -d: -f2)
  PREFIX="${TMP_DIR}/${TISSUE}_${CT}"
  LABEL="${TISSUE}_${CT}"

  echo ""
  echo "=== ${LABEL} ==="

  run_homer \
    "${PREFIX}_closing.bed" \
    "${PREFIX}_opening.bed" \
    "${HOMER_OUT}/${LABEL}__closing_vs_opening" \
    "${LABEL} closing_vs_opening"
done

# ── Step 3: Plot ──────────────────────────────────────────────
echo "$(date)  Step 3: Plotting..."

Rscript - <<'REOF'
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

HOMER_OUT <- "/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_closing_vs_opening/HOMER"
OUT_PLOT  <- "/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_closing_vs_opening/figures"
TOP_N     <- 200   # take all, filter by significance below
CAP_P     <- 50

to_num      <- function(x) suppressWarnings(as.numeric(gsub(",|%","",as.character(x))))
short_motif <- function(x) sapply(strsplit(as.character(x),"/"),`[`,1)

known_files <- list.files(HOMER_OUT, pattern="^knownResults\\.txt$",
                          recursive=TRUE, full.names=TRUE)
message("Found ", length(known_files), " knownResults.txt files")

all_rows <- list()
for (fp in known_files) {
  dname <- basename(dirname(fp))
  parts <- strsplit(dname, "__", fixed=TRUE)[[1]]
  if (length(parts) < 2) next
  tissue_ct <- parts[1]

  tab <- tryCatch(
    read.delim(fp, header=TRUE, sep="\t",
               stringsAsFactors=FALSE, check.names=FALSE),
    error=function(e) NULL
  )
  if (is.null(tab) || nrow(tab)==0) next
  if (!all(c("Motif Name","P-value") %in% colnames(tab))) next

  tab$pval      <- to_num(tab[["P-value"]])
  tab$log10p    <- -log10(pmax(tab$pval, 1e-300))
  tab$motif     <- short_motif(tab[["Motif Name"]])
  tab$tissue_ct <- tissue_ct

  all_rows[[length(all_rows)+1]] <-
    tab[, c("motif","log10p","pval","tissue_ct")]
}

if (length(all_rows)==0) stop("No HOMER results found.")
df <- dplyr::bind_rows(all_rows)
message("Loaded: ", nrow(df), " rows, ", length(unique(df$tissue_ct)), " cell types")

# ── Dotplot: all cell types together ─────────────────────────
DOT_N <- 10  # top N per cell type for dotplot readability

dat_dedup <- df %>%
  group_by(tissue_ct, motif) %>%
  arrange(pval) %>%
  slice(1) %>%
  ungroup()

top_df <- dat_dedup %>%
  group_by(tissue_ct) %>%
  filter(pval < 0.05) %>%
  slice_max(order_by=log10p, n=DOT_N, with_ties=FALSE) %>%
  ungroup() %>%
  mutate(log10p_plot = pmin(log10p, CAP_P))

motif_ord <- top_df %>%
  group_by(motif) %>%
  summarise(m=max(log10p_plot), .groups="drop") %>%
  arrange(desc(m)) %>% pull(motif)
top_df$motif     <- factor(top_df$motif, levels=rev(motif_ord))
top_df$tissue_ct <- factor(top_df$tissue_ct, levels=sort(unique(top_df$tissue_ct)))

n_m <- length(unique(top_df$motif))
n_c <- length(unique(top_df$tissue_ct))
h   <- max(6, n_m*0.28 + 2)
w   <- max(8, n_c*1.3 + 3)

p_dot <- ggplot(top_df, aes(x=tissue_ct, y=motif,
                              size=log10p_plot, color=log10p_plot)) +
  geom_point() +
  scale_size_continuous(range=c(1,8), name="-log10(P)", limits=c(0,CAP_P)) +
  scale_color_gradient(low="grey88", high="#2166AC",
                       name="-log10(P)", limits=c(0,CAP_P)) +
  labs(title="TF motifs enriched in closing peaks\n(vs opening peaks — cell identity TFs)",
       x=NULL, y="TF Motif") +
  theme_bw(base_size=12) +
  theme(
    axis.text.x      = element_text(angle=45, hjust=1, size=10),
    axis.text.y      = element_text(size=8),
    plot.title       = element_text(face="bold", hjust=0.5),
    panel.grid.major = element_line(color="grey92")
  )

ggsave(file.path(OUT_PLOT, "CellIdentity_closing_vs_opening_dotplot.pdf"),
       p_dot, width=w, height=h, limitsize=FALSE)
ggsave(file.path(OUT_PLOT, "CellIdentity_closing_vs_opening_dotplot.png"),
       p_dot, width=w, height=h, dpi=300, limitsize=FALSE)
message("Saved dotplot")

# ── Barplot per cell type ─────────────────────────────────────
SHOW_N <- 30  # number of motifs to display in barplot
for (ct in sort(unique(df$tissue_ct))) {
  dat_ct <- dat_dedup %>%
    filter(tissue_ct==ct) %>%
    filter(pval < 0.05) %>%                              # significant only
    slice_max(order_by=log10p, n=SHOW_N, with_ties=FALSE) %>%
    mutate(log10p_plot=pmin(log10p, CAP_P))

  if (nrow(dat_ct)==0) next

  motif_ord_ct <- dat_ct %>%
    arrange(log10p_plot) %>% pull(motif)
  dat_ct$motif <- factor(dat_ct$motif, levels=motif_ord_ct)

  p_bar <- ggplot(dat_ct, aes(x=log10p_plot, y=motif)) +
    geom_col(fill="#2166AC") +
    geom_vline(xintercept=-log10(0.05), linetype="dashed",
               color="grey50", linewidth=0.4) +
    labs(
      title    = paste0(ct, " — closing vs opening (cell identity TFs)"),
      subtitle = "Foreground: closing peaks | Background: opening peaks",
      x="-log10(P)", y="TF Motif"
    ) +
    theme_bw(base_size=11) +
    theme(
      axis.text.y   = element_text(size=8),
      plot.title    = element_text(face="bold", hjust=0.5),
      plot.subtitle = element_text(hjust=0.5, size=9, color="grey40")
    )

  safe_ct <- gsub("[^A-Za-z0-9_]","_", ct)
  h_ct <- max(5, nrow(dat_ct)*0.28 + 2)
  ggsave(file.path(OUT_PLOT, paste0(safe_ct,"_closing_vs_opening_barplot.pdf")),
         p_bar, width=8, height=h_ct, limitsize=FALSE)
  ggsave(file.path(OUT_PLOT, paste0(safe_ct,"_closing_vs_opening_barplot.png")),
         p_bar, width=8, height=h_ct, dpi=300, limitsize=FALSE)
  message("Saved bar: ", ct)
}

message("All done. Output: ", OUT_PLOT)
REOF

echo "$(date)  All done."

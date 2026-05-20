#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --job-name=region_classification
#SBATCH --time=02:00:00
#SBATCH --partition=general
#SBATCH --account=a_imb_ccbcd
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

set -euo pipefail

module load bedtools/2.30.0-gcc-11.3.0
module load r/4.4.2

BASE="/QRISdata/Q8448/Mouse_disease_data/DAR"
TMP_DIR="${BASE}/DAR_science_comparison/tmp_beds"
OUT_DIR="${BASE}/DAR_science_comparison/results"
CLASS_DIR="${OUT_DIR}/classified_regions"
mkdir -p "${CLASS_DIR}"

# ════════════════════════════════════════════════════════════
# Step 1: Classify peaks into 6 categories per cell type
#   disease_specific_open  : disease open, NOT overlapping aging open
#   shared_open            : disease open AND aging open  (same direction)
#   aging_specific_open    : aging open,   NOT overlapping disease open
#   disease_specific_close : disease close, NOT overlapping aging close
#   shared_close           : disease close AND aging close (same direction)
#   aging_specific_close   : aging close,  NOT overlapping disease close
# ════════════════════════════════════════════════════════════
echo "$(date)  Step 1: Classifying regions..."

for aging_open_bed in "${TMP_DIR}"/*_aging_open.bed; do
    base=$(basename "${aging_open_bed}" _aging_open.bed)

    dis_open="${TMP_DIR}/${base}_disease_open.bed"
    dis_close="${TMP_DIR}/${base}_disease_close.bed"
    age_open="${TMP_DIR}/${base}_aging_open.bed"
    age_close="${TMP_DIR}/${base}_aging_close.bed"

    # Skip if no disease file (cell type not in our dataset)
    [ -f "${dis_open}" ] || continue

    n_do=$(wc -l < "${dis_open}")
    n_dc=$(wc -l < "${dis_close}" 2>/dev/null || echo 0)

    # ── Opening direction ──────────────────────────────────
    if [ "${n_do}" -gt 0 ]; then
        bedtools intersect -a "${dis_open}"  -b "${age_open}" -v \
            > "${CLASS_DIR}/${base}_disease_specific_open.bed"
        bedtools intersect -a "${dis_open}"  -b "${age_open}" -u \
            > "${CLASS_DIR}/${base}_shared_open.bed"
        bedtools intersect -a "${age_open}"  -b "${dis_open}" -v \
            > "${CLASS_DIR}/${base}_aging_specific_open.bed"
    else
        touch "${CLASS_DIR}/${base}_disease_specific_open.bed"
        touch "${CLASS_DIR}/${base}_shared_open.bed"
        touch "${CLASS_DIR}/${base}_aging_specific_open.bed"
    fi

    # ── Closing direction ──────────────────────────────────
    if [ "${n_dc}" -gt 0 ] && [ -s "${age_close}" ]; then
        bedtools intersect -a "${dis_close}" -b "${age_close}" -v \
            > "${CLASS_DIR}/${base}_disease_specific_close.bed"
        bedtools intersect -a "${dis_close}" -b "${age_close}" -u \
            > "${CLASS_DIR}/${base}_shared_close.bed"
        bedtools intersect -a "${age_close}" -b "${dis_close}" -v \
            > "${CLASS_DIR}/${base}_aging_specific_close.bed"
    else
        touch "${CLASS_DIR}/${base}_disease_specific_close.bed"
        touch "${CLASS_DIR}/${base}_shared_close.bed"
        touch "${CLASS_DIR}/${base}_aging_specific_close.bed"
    fi

    n_dso=$(wc -l < "${CLASS_DIR}/${base}_disease_specific_open.bed")
    n_sho=$(wc -l < "${CLASS_DIR}/${base}_shared_open.bed")
    n_aso=$(wc -l < "${CLASS_DIR}/${base}_aging_specific_open.bed")
    n_dsc=$(wc -l < "${CLASS_DIR}/${base}_disease_specific_close.bed")
    n_shc=$(wc -l < "${CLASS_DIR}/${base}_shared_close.bed")
    n_asc=$(wc -l < "${CLASS_DIR}/${base}_aging_specific_close.bed")
    echo "  ${base}: dis_specific_open=${n_dso}  shared_open=${n_sho}  aging_specific_open=${n_aso} | dis_specific_close=${n_dsc}  shared_close=${n_shc}  aging_specific_close=${n_asc}"
done

# ════════════════════════════════════════════════════════════
# Step 2: Compile statistics table (Python)
# ════════════════════════════════════════════════════════════
echo ""
echo "$(date)  Step 2: Compiling statistics table..."

PYTHON="/home/s4869245/.conda/envs/scanpy_env/bin/python"

${PYTHON} - <<'PYEOF'
import os
import pandas as pd

CLASS_DIR = "/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_science_comparison/results/classified_regions"
OUT_DIR   = "/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_science_comparison/results"

def count_bed(f):
    if not os.path.exists(f): return 0
    with open(f) as fh:
        return sum(1 for l in fh if l.strip())

categories = [
    "disease_specific_open",
    "shared_open",
    "aging_specific_open",
    "disease_specific_close",
    "shared_close",
    "aging_specific_close",
]

rows = []
bases = sorted(set(
    f.replace("_disease_specific_open.bed","")
     .replace("_shared_open.bed","")
     .replace("_disease_specific_close.bed","")
     .replace("_shared_close.bed","")
     .replace("_aging_specific_open.bed","")
     .replace("_aging_specific_close.bed","")
    for f in os.listdir(CLASS_DIR)
    if f.endswith(".bed")
))

for base in bases:
    tissue = base.split("_")[0]
    ct     = base[len(tissue)+1:]
    row    = dict(tissue=tissue, cell_type=ct)
    total_disease_open  = 0
    total_disease_close = 0
    for cat in categories:
        n = count_bed(f"{CLASS_DIR}/{base}_{cat}.bed")
        row[cat] = n
        if "open"  in cat and cat != "aging_specific_open":  total_disease_open  += n
        if "close" in cat and cat != "aging_specific_close": total_disease_close += n
    row["total_disease_open"]  = row["disease_specific_open"]  + row["shared_open"]
    row["total_disease_close"] = row["disease_specific_close"] + row["shared_close"]
    row["pct_open_shared"]  = (row["shared_open"]  / row["total_disease_open"]  * 100
                                if row["total_disease_open"]  > 0 else 0)
    row["pct_close_shared"] = (row["shared_close"] / row["total_disease_close"] * 100
                                if row["total_disease_close"] > 0 else 0)
    rows.append(row)

df = pd.DataFrame(rows)
out_csv = f"{OUT_DIR}/region_classification_summary.csv"
df.to_csv(out_csv, index=False)
print(f"Saved: {out_csv}")
print(df[["tissue","cell_type",
          "disease_specific_open","shared_open","aging_specific_open",
          "disease_specific_close","shared_close","aging_specific_close",
          "pct_open_shared","pct_close_shared"]].to_string(index=False))
PYEOF

# ════════════════════════════════════════════════════════════
# Step 3: Plots in R
#   A) Install packages if missing
#   B) Stacked bar (counts, opening + closing per cell type)
#   C) Proportion bar (% shared vs specific)
#   D) UpSet plot per tissue (ComplexHeatmap)
# ════════════════════════════════════════════════════════════
echo ""
echo "$(date)  Step 3: Installing R packages and plotting..."

Rscript - <<'REOF'

# ── A) Install missing packages ────────────────────────────
message("Checking / installing R packages...")

pkgs_cran <- c("ggplot2","dplyr","tidyr","patchwork","scales",
               "RColorBrewer","ggforce")
new_cran <- pkgs_cran[!sapply(pkgs_cran, requireNamespace, quietly = TRUE)]
if (length(new_cran) > 0) {
  message("  Installing from CRAN: ", paste(new_cran, collapse=", "))
  install.packages(new_cran, repos = "https://cloud.r-project.org", quiet = TRUE)
}

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org", quiet = TRUE)
if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  message("  Installing ComplexHeatmap from Bioconductor...")
  BiocManager::install("ComplexHeatmap", ask = FALSE, update = FALSE)
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(scales)
  library(ComplexHeatmap)
})

CLASS_DIR <- "/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_science_comparison/results/classified_regions"
OUT_DIR   <- "/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_science_comparison/results"
FIG_DIR   <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

df <- read.csv(file.path(OUT_DIR, "region_classification_summary.csv"),
               stringsAsFactors = FALSE)

# ── Keep only the selected cell types ──
KEEP <- list(
  Kidney = c("DCT", "PT", "TAL"),
  Lung   = c("AT2", "Fib", "Mac-alv")
)
df <- df %>%
  filter(mapply(function(tis, ct) ct %in% KEEP[[tis]], tissue, cell_type))
message("Cell types retained for plots: ", nrow(df))

# ── Fixed cell type order: Kidney first, then Lung, alpha within ──
ct_order <- df %>%
  arrange(tissue, cell_type) %>%
  mutate(label = paste0(tissue, "\n", cell_type)) %>%
  pull(label) %>% unique()
df$label <- factor(paste0(df$tissue, "\n", df$cell_type), levels = ct_order)

# Colour palette
cat_colors <- c(
  "Disease-specific open"  = "#D6604D",
  "Shared open"            = "#8B0000",
  "Aging-specific open"    = "#F4A582",
  "Disease-specific close" = "#4393C3",
  "Shared close"           = "#08306B",
  "Aging-specific close"   = "#9ECAE1"
)
cat_levels <- names(cat_colors)

# ── B) Stacked bar — absolute counts ──────────────────────
long_df <- df %>%
  select(label, tissue,
         `Disease-specific open`  = disease_specific_open,
         `Shared open`            = shared_open,
         `Aging-specific open`    = aging_specific_open,
         `Disease-specific close` = disease_specific_close,
         `Shared close`           = shared_close,
         `Aging-specific close`   = aging_specific_close) %>%
  pivot_longer(cols = -c(label, tissue),
               names_to = "category", values_to = "n") %>%
  mutate(category = factor(category, levels = cat_levels))

p_stack <- ggplot(long_df, aes(x = label, y = n, fill = category)) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = cat_colors) +
  scale_y_continuous(labels = comma) +
  facet_grid(. ~ tissue, scales = "free_x", space = "free_x") +
  labs(
    title    = "DAR classification: Disease-specific vs Shared with aging",
    subtitle = "Shared = peaks overlapping (≥1 bp) aging DARs in the same direction",
    x = NULL, y = "Number of peaks", fill = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 9),
    plot.title       = element_text(face = "bold", hjust = 0.5),
    plot.subtitle    = element_text(hjust = 0.5, size = 9, color = "grey40"),
    strip.background = element_rect(fill = "grey90"),
    strip.text       = element_text(face = "bold"),
    legend.position  = "bottom",
    legend.text      = element_text(size = 9)
  ) +
  guides(fill = guide_legend(nrow = 2))

ggsave(file.path(FIG_DIR, "region_stacked_counts.pdf"),  p_stack, width = 14, height = 6)
ggsave(file.path(FIG_DIR, "region_stacked_counts.png"),  p_stack, width = 14, height = 6, dpi = 300)
message("Saved: region_stacked_counts")

# ── C) Proportion bar — % of disease DARs that are shared ─
prop_df <- df %>%
  select(label, tissue,
         pct_open_shared, pct_close_shared) %>%
  pivot_longer(cols = c(pct_open_shared, pct_close_shared),
               names_to = "direction", values_to = "pct") %>%
  mutate(
    direction = ifelse(direction == "pct_open_shared", "Opening", "Closing"),
    direction = factor(direction, levels = c("Opening","Closing"))
  )

dir_colors2 <- c("Opening" = "#B2182B", "Closing" = "#2166AC")

p_prop <- ggplot(prop_df, aes(x = label, y = pct, fill = direction)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_hline(yintercept = c(25, 50, 75), linetype = "dashed",
             color = "grey70", linewidth = 0.3) +
  scale_fill_manual(values = dir_colors2) +
  scale_y_continuous(limits = c(0, 105), labels = function(x) paste0(x, "%")) +
  facet_grid(. ~ tissue, scales = "free_x", space = "free_x") +
  labs(
    title    = "% of disease DARs shared with aging DARs",
    subtitle = "Denominator = total significant disease DARs per direction per cell type",
    x = NULL, y = "% shared with aging", fill = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 9),
    plot.title       = element_text(face = "bold", hjust = 0.5),
    plot.subtitle    = element_text(hjust = 0.5, size = 9, color = "grey40"),
    strip.background = element_rect(fill = "grey90"),
    strip.text       = element_text(face = "bold"),
    legend.position  = "top"
  )

ggsave(file.path(FIG_DIR, "region_proportion_shared.pdf"), p_prop, width = 14, height = 6)
ggsave(file.path(FIG_DIR, "region_proportion_shared.png"), p_prop, width = 14, height = 6, dpi = 300)
message("Saved: region_proportion_shared")

# ── D) Shared DAR counts bar chart ────────────────────────
shared_df <- df %>%
  select(label, tissue, shared_open, shared_close) %>%
  pivot_longer(cols = c(shared_open, shared_close),
               names_to = "direction", values_to = "n") %>%
  mutate(
    direction = ifelse(direction == "shared_open", "Opening", "Closing"),
    direction = factor(direction, levels = c("Opening", "Closing"))
  )

p_shared <- ggplot(shared_df, aes(x = label, y = n, fill = direction)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_text(aes(label = n, group = direction),
            position = position_dodge(0.8),
            vjust = -0.3, size = 3, color = "grey30") +
  scale_fill_manual(values = dir_colors2) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.12))) +
  facet_grid(. ~ tissue, scales = "free_x", space = "free_x") +
  labs(
    title    = "Number of disease DARs shared with aging DARs",
    subtitle = "Shared = disease DAR overlapping aging DAR in the same direction (≥1 bp)",
    x = NULL, y = "Number of shared peaks", fill = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 9),
    plot.title       = element_text(face = "bold", hjust = 0.5),
    plot.subtitle    = element_text(hjust = 0.5, size = 9, color = "grey40"),
    strip.background = element_rect(fill = "grey90"),
    strip.text       = element_text(face = "bold"),
    legend.position  = "top"
  )

ggsave(file.path(FIG_DIR, "region_shared_counts.pdf"), p_shared, width = 14, height = 6)
ggsave(file.path(FIG_DIR, "region_shared_counts.png"), p_shared, width = 14, height = 6, dpi = 300)
message("Saved: region_shared_counts")

# ── E) UpSet plot per tissue (ComplexHeatmap) ──────────────
# For each cell type with enough peaks, build a binary matrix:
#   rows = peaks (union of disease & aging significant peaks)
#   cols = {disease_open, aging_open, disease_close, aging_close}
# Then use make_comb_mat() + UpSet()

read_bed <- function(f) {
  if (!file.exists(f)) return(character(0))
  lines <- readLines(f)
  lines <- lines[nzchar(lines)]
  if (length(lines) == 0) return(character(0))
  # Return "chr:start-end" as peak ID
  sapply(strsplit(lines, "\t"), function(x) paste0(x[1],":",x[2],"-",x[3]))
}

for (tis in unique(df$tissue)) {
  sub_df <- df %>% filter(tissue == tis)
  if (nrow(sub_df) == 0) next

  # Collect all peaks across cell types for this tissue
  all_peak_sets <- list()
  for (i in seq_len(nrow(sub_df))) {
    ct   <- sub_df$cell_type[i]
    base <- paste0(tis, "_", ct)

    do_peaks  <- read_bed(file.path(CLASS_DIR, paste0(base, "_disease_specific_open.bed")))
    sho_peaks <- read_bed(file.path(CLASS_DIR, paste0(base, "_shared_open.bed")))
    dc_peaks  <- read_bed(file.path(CLASS_DIR, paste0(base, "_disease_specific_close.bed")))
    shc_peaks <- read_bed(file.path(CLASS_DIR, paste0(base, "_shared_close.bed")))
    ao_peaks  <- read_bed(file.path(CLASS_DIR, paste0(base, "_aging_specific_open.bed")))
    ac_peaks  <- read_bed(file.path(CLASS_DIR, paste0(base, "_aging_specific_close.bed")))

    disease_open  <- c(do_peaks,  sho_peaks)
    disease_close <- c(dc_peaks,  shc_peaks)
    aging_open    <- c(sho_peaks, ao_peaks)
    aging_close   <- c(shc_peaks, ac_peaks)

    all_peak_sets[[paste0(ct," (dis-open)")]]   <- disease_open
    all_peak_sets[[paste0(ct," (aging-open)")]] <- aging_open
  }

  # Keep only sets with ≥5 peaks
  all_peak_sets <- all_peak_sets[sapply(all_peak_sets, length) >= 5]
  if (length(all_peak_sets) < 2) {
    message("  SKIP UpSet for ", tis, ": too few non-empty sets"); next
  }

  # Build binary matrix for ComplexHeatmap
  universe <- unique(unlist(all_peak_sets))
  if (length(universe) < 10) next

  mat <- sapply(all_peak_sets, function(s) as.integer(universe %in% s))
  rownames(mat) <- universe

  comb_mat <- make_comb_mat(mat)

  pdf(file.path(FIG_DIR, paste0("upset_", tis, ".pdf")), width = 14, height = 6)
  ht <- UpSet(
    comb_mat,
    set_order    = names(all_peak_sets),
    comb_order   = order(comb_size(comb_mat), decreasing = TRUE),
    top_annotation = upset_top_annotation(comb_mat, add_numbers = TRUE),
    left_annotation = upset_left_annotation(comb_mat, add_numbers = TRUE),
    column_title = paste0(tis, " — Disease vs Aging DAR overlap (opening)"),
    row_names_gp = gpar(fontsize = 8)
  )
  draw(ht)
  dev.off()
  message("Saved: upset_", tis, ".pdf")
}

message("\nAll classification plots done.")
message("Output: ", FIG_DIR)
REOF

echo "$(date)  All done."

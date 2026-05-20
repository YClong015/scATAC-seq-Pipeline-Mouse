#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --job-name=TF_motif_plots
#SBATCH --time=02:00:00
#SBATCH --partition=general
#SBATCH --account=a_imb_ccbcd
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

set -euo pipefail

module load r/4.4.2

OUT_DIR="/QRISdata/Q8448/Mouse_disease_data/DAR/TF_motif_plots"
mkdir -p "${OUT_DIR}"

Rscript - <<'REOF'

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

# ============================================================
# Config: same 4 tissues as heatmap
# ============================================================
cfg <- list(
  list(
    tissue    = "Lung",
    contrast  = "Case_vs_Control",
    base_dir  = "/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_pseudobulk_DESeq2/DAR_pseudobulk_Lung_DESeq2",
    cells     = c("AT2","B","Ciliated","EC-vasc","Fib","Mac-alv","Mac-inter","NK","SMCs","T")
  ),
  list(
    tissue    = "Kidney",
    contrast  = "Day42_vs_Sham",
    base_dir  = "/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_pseudobulk_DESeq2/DAR_pseudobulk_Kidney_DESeq2",
    cells     = c("DCT","Macrophages","PC","PT","TAL")
  ),
  list(
    tissue    = "Aorta",
    contrast  = "Challenge_vs_Control",
    base_dir  = "/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_pseudobulk_DESeq2/DAR_pseudobulk_Aorta_DESeq2",
    cells     = c("Macrophages","Pericytes","SMC","SMCs")
  ),
  list(
    tissue    = "Tcell",
    contrast  = "Young_chronic_vs_Young_control",
    base_dir  = "/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_pseudobulk_DESeq2/DAR_pseudobulk_Tcells_DESeq2",
    cells     = c("Tcell")
  )
)

# ── background mode: use stable ─────────────────────────────
MODE     <- "stable"   # change to "NS" if preferred
HOMER_SUBDIR <- ifelse(MODE == "stable", "HOMER_stable_bg", "HOMER_NS_bg")
BG_NAME  <- ifelse(MODE == "stable", "stable", "NS")
TOP_N    <- 40         # top motifs per cell type per direction

# ============================================================
# Helpers
# ============================================================
to_num <- function(x) suppressWarnings(as.numeric(gsub(",|%", "", as.character(x))))
short_motif <- function(x) sapply(strsplit(as.character(x), "/"), `[`, 1)
normalize_cell <- function(x) sub("_+$", "", x)

parse_dirname <- function(dname) {
  parts <- strsplit(dname, "__", fixed = TRUE)[[1]]
  if (length(parts) < 4) return(NULL)
  data.frame(
    cell_type  = normalize_cell(parts[1]),
    contrast   = parts[2],
    thresh     = parts[3],
    comparison = parts[4],
    stringsAsFactors = FALSE
  )
}

# ============================================================
# Read all HOMER results
# ============================================================
valid_opening <- paste0("opening_vs_", BG_NAME)
valid_closing <- paste0("closing_vs_", BG_NAME)

all_rows <- list()

for (x in cfg) {
  homer_dir <- file.path(x$base_dir, HOMER_SUBDIR)
  if (!dir.exists(homer_dir)) {
    message("SKIP (no dir): ", homer_dir); next
  }
  known_files <- list.files(homer_dir, pattern = "^knownResults\\.txt$",
                            recursive = TRUE, full.names = TRUE)
  message("Tissue=", x$tissue, "  files=", length(known_files))

  for (fp in known_files) {
    meta <- parse_dirname(basename(dirname(fp)))
    if (is.null(meta)) next
    if (meta$contrast  != x$contrast)              next
    if (!(meta$cell_type %in% x$cells))            next
    if (!(meta$comparison %in% c(valid_opening, valid_closing))) next

    tab <- tryCatch(
      read.delim(fp, header = TRUE, sep = "\t",
                 stringsAsFactors = FALSE, check.names = FALSE),
      error = function(e) NULL
    )
    if (is.null(tab) || nrow(tab) == 0)            next
    if (!all(c("Motif Name","P-value") %in% colnames(tab))) next

    tab$pval      <- to_num(tab[["P-value"]])
    tab$log10p    <- -log10(pmax(tab$pval, 1e-300))
    tab$motif     <- short_motif(tab[["Motif Name"]])
    tab$tissue    <- x$tissue
    tab$cell_type <- meta$cell_type
    tab$tissue_ct <- ifelse(x$tissue == meta$cell_type,
                            x$tissue,
                            paste0(x$tissue, "_", meta$cell_type))
    tab$direction <- ifelse(meta$comparison == valid_opening, "Opening", "Closing")

    all_rows[[length(all_rows)+1]] <-
      tab[, c("motif","log10p","pval","tissue","cell_type","tissue_ct","direction")]
  }
}

if (length(all_rows) == 0) stop("No HOMER results found.")
df <- bind_rows(all_rows)
message("Total rows loaded: ", nrow(df))

out_dir <- "/QRISdata/Q8448/Mouse_disease_data/DAR/TF_motif_plots"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. Cross-tissue Dot Plot
# ============================================================
message("Building cross-tissue dot plot...")

# Get top N motifs per tissue_ct per direction
top_df <- df %>%
  group_by(tissue_ct, direction) %>%
  slice_max(order_by = log10p, n = TOP_N, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(log10p_plot = pmin(log10p, 50))

# Count how many cell types each motif appears in (per direction)
motif_ct_count <- top_df %>%
  group_by(motif, direction) %>%
  summarise(n_celltypes = n_distinct(tissue_ct), .groups = "drop")
top_df <- top_df %>% left_join(motif_ct_count, by = c("motif","direction"))

# Split into specific (1-2 cell types) and shared (>=3 cell types)
specific_df <- top_df %>% filter(n_celltypes <= 2)
shared_df   <- top_df %>% filter(n_celltypes >= 3)

# Fixed column order: tissue → cell type
ct_order_all <- c()
for (x in cfg) {
  for (ct in x$cells) {
    label <- ifelse(x$tissue == ct, x$tissue, paste0(x$tissue, "_", ct))
    ct_order_all <- c(ct_order_all, label)
  }
}
ct_order_all <- unique(ct_order_all)

make_dot_plot <- function(dat, title_str, color_high) {
  if (nrow(dat) == 0) return(NULL)

  dat$tissue_ct <- factor(dat$tissue_ct,
    levels = ct_order_all[ct_order_all %in% unique(dat$tissue_ct)])

  motif_ord <- dat %>%
    group_by(motif) %>%
    summarise(max_p = max(log10p_plot)) %>%
    arrange(desc(max_p)) %>% pull(motif)
  dat$motif <- factor(dat$motif, levels = rev(motif_ord))

  n_motifs <- length(unique(dat$motif))
  n_ct     <- length(unique(dat$tissue_ct))
  h <- max(6,  n_motifs * 0.22 + 2)
  w <- max(10, n_ct * 0.7 + 4)

  p <- ggplot(dat, aes(x = tissue_ct, y = motif,
                       size = log10p_plot, color = log10p_plot)) +
    geom_point() +
    scale_size_continuous(range = c(0.5, 7), name = "-log10(P)") +
    scale_color_gradient(low = "grey88", high = color_high, name = "-log10(P)") +
    labs(title = title_str, x = NULL, y = "TF Motif") +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x      = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y      = element_text(size = 7),
      plot.title       = element_text(face = "bold", hjust = 0.5),
      panel.grid.major = element_line(color = "grey92")
    )
  list(plot = p, h = h, w = w)
}

save_dot <- function(dat, title_str, color_high, suffix) {
  if (nrow(dat) == 0) { message("No data: ", suffix); return() }
  obj <- make_dot_plot(dat, title_str, color_high)
  ggsave(file.path(out_dir, paste0(suffix, ".pdf")), obj$plot,
         width = obj$w, height = obj$h, limitsize = FALSE)
  ggsave(file.path(out_dir, paste0(suffix, ".png")), obj$plot,
         width = obj$w, height = obj$h, dpi = 300, limitsize = FALSE)
  message("Saved: ", suffix)
}

# ── Opening: specific (1-2 cell types) ──────────────────────
save_dot(
  filter(specific_df, direction == "Opening"),
  paste0("Cell type-SPECIFIC TF motifs — Opening (present in 1-2 cell types)"),
  "#B2182B", "Opening_specific_dotplot"
)

# ── Opening: shared (>=3 cell types) ────────────────────────
save_dot(
  filter(shared_df, direction == "Opening"),
  paste0("SHARED TF motifs — Opening (present in ≥3 cell types)"),
  "#B2182B", "Opening_shared_dotplot"
)

# ── Closing: specific ───────────────────────────────────────
save_dot(
  filter(specific_df, direction == "Closing"),
  paste0("Cell type-SPECIFIC TF motifs — Closing (present in 1-2 cell types)"),
  "#2166AC", "Closing_specific_dotplot"
)

# ── Closing: shared ─────────────────────────────────────────
save_dot(
  filter(shared_df, direction == "Closing"),
  paste0("SHARED TF motifs — Closing (present in ≥3 cell types)"),
  "#2166AC", "Closing_shared_dotplot"
)

# Combined signed dotplot
combined <- bind_rows(
  filter(plot_df, direction == "Opening") %>% mutate(log10p_signed =  log10p_plot),
  filter(plot_df, direction == "Closing") %>% mutate(log10p_signed = -log10p_plot)
)
motif_ord2 <- combined %>%
  group_by(motif) %>%
  summarise(m = max(abs(log10p_plot))) %>%
  arrange(desc(m)) %>% pull(motif)
combined$motif <- factor(combined$motif, levels = rev(motif_ord2))
combined$tissue_ct <- factor(combined$tissue_ct, levels = ct_order_all)

n_m2 <- length(unique(combined$motif))
n_c2 <- length(unique(combined$tissue_ct))
h2 <- max(7, n_m2 * 0.22 + 2)
w2 <- max(12, n_c2 * 0.7 + 5)

p_comb <- ggplot(combined,
                 aes(x = tissue_ct, y = motif,
                     size = abs(log10p_plot), color = log10p_signed)) +
  geom_point() +
  scale_size_continuous(range = c(0.5, 7), name = "-log10(P)") +
  scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                        midpoint = 0, name = "signed\n-log10(P)") +
  facet_wrap(~ direction, nrow = 1) +
  labs(title = paste0("All tissues — TF motifs (red=Opening, blue=Closing) [", MODE, " bg]"),
       x = NULL, y = "TF Motif") +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y      = element_text(size = 7),
    plot.title       = element_text(face = "bold", hjust = 0.5),
    strip.text       = element_text(face = "bold"),
    panel.grid.major = element_line(color = "grey92")
  )

ggsave(file.path(out_dir, "AllTissues_Combined_dotplot.pdf"),
       p_comb, width = w2, height = h2, limitsize = FALSE)
ggsave(file.path(out_dir, "AllTissues_Combined_dotplot.png"),
       p_comb, width = w2, height = h2, dpi = 300, limitsize = FALSE)
message("Saved: AllTissues_Combined_dotplot")

# ============================================================
# 2. Per-cell-type Bar Plots (top 40 motifs, opening & closing)
# ============================================================
message("Building per-cell-type bar plots...")

bar_df <- df %>%
  group_by(tissue_ct, direction) %>%
  slice_max(order_by = log10p, n = TOP_N, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(log10p_plot = pmin(log10p, 50))

for (ct_label in unique(bar_df$tissue_ct)) {
  dat_ct <- bar_df %>% filter(tissue_ct == ct_label)
  if (nrow(dat_ct) == 0) next

  # Order motifs: opening top, then closing top, no duplicates
  open_motifs  <- dat_ct %>% filter(direction == "Opening") %>%
    arrange(desc(log10p_plot)) %>% pull(motif)
  close_motifs <- dat_ct %>% filter(direction == "Closing") %>%
    arrange(desc(log10p_plot)) %>% pull(motif)
  motif_ord_ct <- unique(c(open_motifs, close_motifs))
  dat_ct$motif <- factor(dat_ct$motif, levels = rev(motif_ord_ct))

  dir_colors <- c(Opening = "#B2182B", Closing = "#2166AC")

  p_bar <- ggplot(dat_ct, aes(x = log10p_plot, y = motif, fill = direction)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = dir_colors) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed",
               color = "grey50", linewidth = 0.4) +
    labs(
      title = paste0(ct_label, " — Top ", TOP_N, " TF motifs (", MODE, " bg)"),
      x = "-log10(P)", y = "TF Motif", fill = "Direction"
    ) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.y  = element_text(size = 7),
      plot.title   = element_text(face = "bold", hjust = 0.5),
      legend.position = "top"
    )

  n_motifs_ct <- length(unique(dat_ct$motif))
  h_ct <- max(5, n_motifs_ct * 0.18 + 2)

  safe_label <- gsub("[^A-Za-z0-9_\\-]", "_", ct_label)
  ggsave(file.path(out_dir, paste0(safe_label, "_barplot.pdf")),
         p_bar, width = 8, height = h_ct, limitsize = FALSE)
  ggsave(file.path(out_dir, paste0(safe_label, "_barplot.png")),
         p_bar, width = 8, height = h_ct, dpi = 300, limitsize = FALSE)
  message("Saved bar plot: ", ct_label)
}

message("All done. Output: ", out_dir)
REOF

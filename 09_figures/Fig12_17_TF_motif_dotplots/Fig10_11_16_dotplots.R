#!/usr/bin/env Rscript
# ============================================================
# Fig 10: Shared TF motifs (≥3 cell type/tissue combinations)
# Fig 11: AllTissues comprehensive dotplot
# Fig 16: Tissue-specific TF motifs (unique to 1 tissue)
#
# All panels: ggplot2 dotplot (x = cell type, y = motif,
#             size = target%, colour = -log10(P))
# Faceted by tissue with tissue colour strip.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(scales)
})

# ---------------------------------------------------------------
# 1) Paths and tissue config
# ---------------------------------------------------------------
base_dar <- "/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_pseudobulk_DESeq2"

homer_dirs <- list(
  Lung   = file.path(base_dar, "DAR_pseudobulk_Lung_DESeq2",   "HOMER"),
  Aorta  = file.path(base_dar, "DAR_pseudobulk_Aorta_DESeq2",  "HOMER"),
  Kidney = file.path(base_dar, "DAR_pseudobulk_Kidney_DESeq2", "HOMER"),
  Tcell  = file.path(base_dar, "DAR_pseudobulk_Tcells_DESeq2", "HOMER")
)

tissue_levels <- c("Lung", "Aorta", "Kidney", "Tcell")

# Only retain these contrasts per tissue
contrast_keep <- list(
  Lung   = "Case_vs_Control",
  Aorta  = "Challenge_vs_Control",
  Kidney = c("Day14_vs_Sham", "Day42_vs_Sham"),
  Tcell  = "Young_chronic_vs_Young_control"
)

out_dir <- file.path(
  "/QRISdata/Q8448/Mouse_disease_data/DAR",
  "Combined_HOMER_Heatmap_Focused",
  "TF_motif_plots"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------
# 2) Helpers
# ---------------------------------------------------------------
to_num <- function(x) {
  suppressWarnings(as.numeric(gsub(",", "", gsub("%", "", as.character(x)))))
}

cap_log10 <- function(p, cap = 50) {
  pmin(-log10(pmax(p, 1e-300)), cap)
}

pick_col <- function(nm, pats) {
  hit <- nm[vapply(nm,
    function(x) any(grepl(paste(pats, collapse = "|"), x, ignore.case = TRUE)),
    logical(1))]
  if (length(hit) == 0) NA_character_ else hit[1]
}

short_motif <- function(x) sapply(strsplit(as.character(x), "/"), `[`, 1)

normalize_dir <- function(x) {
  x <- tolower(trimws(x))
  if (x %in% c("opening", "up"))   return("Opening")
  if (x %in% c("closing", "down")) return("Closing")
  NA_character_
}

parse_dirname <- function(dname) {
  parts <- strsplit(dname, "__", fixed = TRUE)[[1]]
  if (length(parts) < 3) return(NULL)
  dir_use <- normalize_dir(parts[3])
  if (is.na(dir_use)) return(NULL)
  data.frame(cell_type = parts[1], contrast = parts[2],
             direction = dir_use, stringsAsFactors = FALSE)
}

read_one <- function(fp, tissue) {
  meta <- parse_dirname(basename(dirname(fp)))
  if (is.null(meta)) return(NULL)

  tab <- tryCatch(
    read.delim(fp, header = TRUE, sep = "\t",
               stringsAsFactors = FALSE, check.names = TRUE),
    error = function(e) NULL
  )
  if (is.null(tab) || nrow(tab) < 1) return(NULL)

  nm         <- colnames(tab)
  motif_col  <- pick_col(nm, "^Motif\\.Name$")
  p_col      <- pick_col(nm, "^P\\.value$")
  q_col      <- pick_col(nm, c("q\\.value", "Benjamini"))
  tgt_pct_col <- pick_col(nm,
    c("^X\\.\\.of\\.Target\\.Sequences\\.with\\.Motif$"))

  if (is.na(motif_col) || is.na(p_col)) return(NULL)

  data.frame(
    tissue      = tissue,
    cell_type   = meta$cell_type,
    contrast    = meta$contrast,
    direction   = meta$direction,
    motif_raw   = tab[[motif_col]],
    p_value     = to_num(tab[[p_col]]),
    q_value     = if (!is.na(q_col))      to_num(tab[[q_col]])      else NA_real_,
    target_pct  = if (!is.na(tgt_pct_col)) to_num(tab[[tgt_pct_col]]) else NA_real_,
    rank        = seq_len(nrow(tab)),
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------
# 3) Load all tissues
# ---------------------------------------------------------------
all_rows <- list()
for (tis in names(homer_dirs)) {
  d <- homer_dirs[[tis]]
  if (!dir.exists(d)) { message("Skip missing: ", d); next }
  fps  <- list.files(d, "^knownResults\\.txt$", recursive = TRUE, full.names = TRUE)
  rows <- Filter(Negate(is.null), lapply(fps, read_one, tissue = tis))
  if (length(rows)) all_rows[[tis]] <- bind_rows(rows)
}
if (length(all_rows) == 0) stop("No HOMER results found.")

df <- bind_rows(all_rows) %>%
  mutate(
    q_value    = ifelse(is.na(q_value),    1, q_value),
    p_value    = ifelse(is.na(p_value),    1, p_value),
    target_pct = ifelse(is.na(target_pct), 0, target_pct),
    log10_p    = cap_log10(p_value),
    motif_label = short_motif(motif_raw),
    tissue      = factor(tissue, levels = tissue_levels)
  )

# Apply contrast filter per tissue
df <- df %>%
  rowwise() %>%
  filter(contrast %in% contrast_keep[[as.character(tissue)]]) %>%
  ungroup()

# Save merged table
write.csv(df, file.path(out_dir, "AllTissues_HOMER_merged.csv"), row.names = FALSE)

# ---------------------------------------------------------------
# 4) Select significant top motifs (unit = tissue × cell_type × contrast × direction)
# ---------------------------------------------------------------
q_cut  <- 0.05
top_n  <- 15   # top N per unit for initial candidate pool

df_sig <- df %>%
  filter(q_value < q_cut) %>%
  group_by(tissue, cell_type, contrast, direction, motif_raw, motif_label) %>%
  slice_min(rank, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  group_by(tissue, cell_type, contrast, direction) %>%
  slice_min(rank, n = top_n, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    unit = paste(tissue, cell_type, contrast, sep = "|"),
    # Column label for x-axis: CellType (no tissue prefix if same as tissue)
    col_label = if_else(
      as.character(tissue) == cell_type,
      as.character(tissue),
      cell_type
    )
  )

# ---------------------------------------------------------------
# 5) Common dotplot helpers
# ---------------------------------------------------------------
dot_theme <- function() {
  theme_bw(base_size = 11) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y     = element_text(size = 8),
    strip.text      = element_text(face = "bold", size = 9),
    strip.background = element_rect(fill = "grey92"),
    panel.grid      = element_blank(),
    plot.title      = element_text(face = "bold", size = 11, hjust = 0),
    legend.position = "right"
  )
}

make_dotplot <- function(plot_df, title_str, col_high = "#B30000",
                         x_var = "col_label") {
  ggplot(plot_df,
         aes(x = .data[[x_var]], y = motif_label,
             size = target_pct, colour = log10_p)) +
    geom_point(alpha = 0.85) +
    facet_grid(. ~ tissue, scales = "free_x", space = "free_x") +
    scale_size(range = c(1.5, 7), name = "Target %") +
    scale_colour_gradient(
      low = "grey85", high = col_high,
      limits = c(0, 50), oob = squish,
      name = "-log10(P)"
    ) +
    labs(title = title_str, x = NULL, y = "Motif") +
    dot_theme()
}

save_pair <- function(p, prefix, w = 14, h = 9) {
  ggsave(file.path(out_dir, paste0(prefix, ".pdf")), p, width = w, height = h)
  ggsave(file.path(out_dir, paste0(prefix, ".png")), p, width = w, height = h, dpi = 300)
}

# ---------------------------------------------------------------
# 6) FIG 10: Shared motifs (≥3 distinct tissue/cell_type units)
# ---------------------------------------------------------------
min_units <- 3

make_shared_panel <- function(dir_use, panel_letter) {
  sub <- df_sig %>% filter(direction == dir_use)
  if (nrow(sub) == 0) return(NULL)

  # How many distinct units does each motif appear in?
  motif_units <- sub %>%
    group_by(motif_raw, motif_label) %>%
    summarise(
      n_units = n_distinct(unit),
      n_tissues = n_distinct(as.character(tissue)),
      best_p  = max(log10_p),
      .groups = "drop"
    ) %>%
    filter(n_units >= min_units) %>%
    arrange(desc(n_tissues), desc(n_units), desc(best_p)) %>%
    slice_head(n = 30)

  if (nrow(motif_units) == 0) {
    message("No shared motifs (≥", min_units, " units) for ", dir_use)
    return(NULL)
  }

  shared_set <- motif_units$motif_raw
  label_order <- rev(motif_units$motif_label)   # best at top

  plot_df <- sub %>%
    filter(motif_raw %in% shared_set) %>%
    group_by(tissue, cell_type, contrast, col_label, motif_raw, motif_label) %>%
    summarise(log10_p = max(log10_p), target_pct = max(target_pct), .groups = "drop") %>%
    mutate(
      motif_label = factor(motif_label, levels = label_order),
      tissue      = factor(tissue, levels = tissue_levels)
    )

  col_high <- if (dir_use == "Opening") "#B30000" else "#2166AC"
  title <- paste0(panel_letter, "  Shared ", dir_use,
                  " motifs (≥", min_units, " cell types/tissues)")
  make_dotplot(plot_df, title, col_high)
}

p_shared_open  <- make_shared_panel("Opening", "A")
p_shared_close <- make_shared_panel("Closing", "B")

if (!is.null(p_shared_open))  save_pair(p_shared_open,  "Fig10_Opening_shared_dotplot",  w = 14, h = 9)
if (!is.null(p_shared_close)) save_pair(p_shared_close, "Fig10_Closing_shared_dotplot",  w = 14, h = 9)
message("Fig 10 saved.")

# ---------------------------------------------------------------
# 7) FIG 11: AllTissues comprehensive dotplot (top 25 per tissue)
# ---------------------------------------------------------------
make_all_panel <- function(dir_use, panel_letter) {
  sub <- df_sig %>% filter(direction == dir_use)
  if (nrow(sub) == 0) return(NULL)

  # Top 25 by best rank within each tissue, deduplicate motifs
  keep <- sub %>%
    group_by(tissue, motif_raw, motif_label) %>%
    summarise(best_rank = min(rank), best_p = max(log10_p), .groups = "drop") %>%
    group_by(tissue) %>%
    slice_min(best_rank, n = 25, with_ties = FALSE) %>%
    ungroup()

  all_motifs <- unique(keep$motif_raw)

  plot_df <- sub %>%
    filter(motif_raw %in% all_motifs) %>%
    group_by(tissue, cell_type, contrast, col_label, motif_raw, motif_label) %>%
    summarise(log10_p = max(log10_p), target_pct = max(target_pct), .groups = "drop") %>%
    mutate(tissue = factor(tissue, levels = tissue_levels))

  # Y-order: decreasing best log10_p overall
  motif_order <- plot_df %>%
    group_by(motif_label) %>%
    summarise(best = max(log10_p), .groups = "drop") %>%
    arrange(best) %>%   # rev() puts highest at top in ggplot discrete scale
    pull(motif_label)

  plot_df <- plot_df %>%
    mutate(motif_label = factor(motif_label, levels = motif_order))

  col_high <- if (dir_use == "Opening") "#B30000" else "#2166AC"
  title <- paste0(panel_letter, "  All tissues — ", dir_use,
                  " DAR motifs (top 25 per tissue)")
  make_dotplot(plot_df, title, col_high)
}

p_all_open  <- make_all_panel("Opening", "A")
p_all_close <- make_all_panel("Closing", "B")

if (!is.null(p_all_open))  save_pair(p_all_open,  "Fig11_AllTissues_Opening_dotplot",  w = 16, h = 11)
if (!is.null(p_all_close)) save_pair(p_all_close, "Fig11_AllTissues_Closing_dotplot",  w = 16, h = 11)
message("Fig 11 saved.")

# ---------------------------------------------------------------
# 8) FIG 16: Tissue-specific motifs (present in exactly 1 tissue)
# ---------------------------------------------------------------
make_specific_panel <- function(dir_use, panel_letter) {
  sub <- df_sig %>% filter(direction == dir_use)
  if (nrow(sub) == 0) return(NULL)

  # Count distinct tissues per motif
  motif_tissue_n <- sub %>%
    group_by(motif_raw) %>%
    summarise(n_tissues = n_distinct(as.character(tissue)), .groups = "drop")

  specific_set <- motif_tissue_n %>%
    filter(n_tissues == 1) %>%
    pull(motif_raw)

  if (length(specific_set) == 0) {
    message("No tissue-specific motifs for ", dir_use)
    return(NULL)
  }

  # Top 15 per tissue by best log10_p
  keep <- sub %>%
    filter(motif_raw %in% specific_set) %>%
    group_by(tissue, motif_raw, motif_label) %>%
    summarise(best_p = max(log10_p), .groups = "drop") %>%
    group_by(tissue) %>%
    slice_max(best_p, n = 15, with_ties = FALSE) %>%
    ungroup()

  keep_motifs <- unique(keep$motif_raw)

  plot_df <- sub %>%
    filter(motif_raw %in% keep_motifs) %>%
    group_by(tissue, cell_type, contrast, col_label, motif_raw, motif_label) %>%
    summarise(log10_p = max(log10_p), target_pct = max(target_pct), .groups = "drop") %>%
    mutate(tissue = factor(tissue, levels = tissue_levels))

  # Y-order: by tissue then decreasing best_p (tissue-grouped appearance)
  motif_order <- keep %>%
    arrange(tissue, desc(best_p)) %>%
    pull(motif_label) %>%
    unique()

  plot_df <- plot_df %>%
    mutate(motif_label = factor(motif_label, levels = rev(motif_order)))

  col_high <- if (dir_use == "Opening") "#B30000" else "#2166AC"
  title <- paste0(panel_letter, "  Tissue-specific ", dir_use,
                  " motifs (unique to 1 tissue)")
  make_dotplot(plot_df, title, col_high)
}

p_spec_open  <- make_specific_panel("Opening", "A")
p_spec_close <- make_specific_panel("Closing", "B")

if (!is.null(p_spec_open))  save_pair(p_spec_open,  "Fig16_Opening_specific_dotplot",  w = 14, h = 9)
if (!is.null(p_spec_close)) save_pair(p_spec_close, "Fig16_Closing_specific_dotplot",  w = 14, h = 9)
message("Fig 16 saved.")

message("All outputs → ", out_dir)

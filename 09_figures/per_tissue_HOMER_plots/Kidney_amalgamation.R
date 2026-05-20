#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

# ============================================================
# Kidney HOMER motif plots
#
# Opening:
#   TRUE amalgamation of top motifs across cell types
#   (union of each cell type's top10 motifs per contrast)
#
# Closing:
#   top 25 motifs for individual cell types
# ============================================================

# -----------------------------
# 1) Paths
# -----------------------------
out_dir <- paste0(
  "/QRISdata/Q8448/Mouse_disease_data/DAR/",
  "DAR_pseudobulk_DESeq2/DAR_pseudobulk_Kidney_DESeq2"
)

homer_dir <- file.path(out_dir, "HOMER")
plot_dir <- file.path(out_dir, "HOMER_Plots")

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 2) Helpers
# -----------------------------
parse_dirname <- function(dname) {
  parts <- strsplit(dname, "__", fixed = TRUE)[[1]]
  if (length(parts) < 3) return(NULL)
  
  data.frame(
    cell_type = gsub("_", " ", parts[1]),
    contrast = parts[2],
    direction = tolower(parts[3]),
    stringsAsFactors = FALSE
  )
}

to_num <- function(x) {
  x <- as.character(x)
  x <- gsub("%", "", x)
  x <- gsub(",", "", x)
  suppressWarnings(as.numeric(x))
}

cap_log10 <- function(p, cap = 300) {
  p <- pmax(p, 1e-300)
  v <- -log10(p)
  pmin(v, cap)
}

pick_col <- function(nm, patterns) {
  hit <- nm[
    vapply(
      nm,
      function(x) {
        any(grepl(
          paste(patterns, collapse = "|"),
          x,
          ignore.case = TRUE
        ))
      },
      logical(1)
    )
  ]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

read_known <- function(fp) {
  dname <- basename(dirname(fp))
  meta <- parse_dirname(dname)
  if (is.null(meta)) return(NULL)
  
  tab <- tryCatch(
    read.delim(
      fp,
      header = TRUE,
      sep = "\t",
      stringsAsFactors = FALSE,
      check.names = TRUE
    ),
    error = function(e) NULL
  )
  
  if (is.null(tab) || nrow(tab) < 1) return(NULL)
  
  nm <- colnames(tab)
  
  motif_col <- pick_col(nm, c("^Motif\\.Name$"))
  p_col <- pick_col(nm, c("^P\\.value$"))
  q_col <- pick_col(nm, c("q\\.value", "Benjamini"))
  
  tgt_n_col <- pick_col(
    nm,
    c("^X\\.of\\.Target\\.Sequences\\.with\\.Motif\\.of\\.")
  )
  tgt_pct_col <- pick_col(
    nm,
    c("^X\\.\\.of\\.Target\\.Sequences\\.with\\.Motif$")
  )
  bg_pct_col <- pick_col(
    nm,
    c("^X\\.\\.of\\.Background\\.Sequences\\.with\\.Motif$")
  )
  
  if (is.na(motif_col) || is.na(p_col)) return(NULL)
  
  out <- data.frame(
    motif_raw = tab[[motif_col]],
    p_value = to_num(tab[[p_col]]),
    q_value = if (!is.na(q_col)) to_num(tab[[q_col]]) else NA_real_,
    target_n = if (!is.na(tgt_n_col)) to_num(tab[[tgt_n_col]])
    else NA_real_,
    target_pct = if (!is.na(tgt_pct_col)) to_num(tab[[tgt_pct_col]])
    else NA_real_,
    bg_pct = if (!is.na(bg_pct_col)) to_num(tab[[bg_pct_col]])
    else NA_real_,
    rank = seq_len(nrow(tab)),
    stringsAsFactors = FALSE
  )
  
  cbind(meta, out)
}

short_motif <- function(x) {
  x <- as.character(x)
  sapply(strsplit(x, "/"), `[`, 1)
}

# -----------------------------
# 3) Read all knownResults.txt
# -----------------------------
known_files <- list.files(
  path = homer_dir,
  pattern = "^knownResults\\.txt$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(known_files) == 0) {
  stop("No knownResults.txt under: ", homer_dir)
}

lst <- lapply(known_files, read_known)
lst <- lst[!vapply(lst, is.null, logical(1))]

if (length(lst) == 0) {
  stop("All knownResults.txt unreadable or empty.")
}

df <- bind_rows(lst) %>%
  mutate(
    direction = ifelse(direction == "opening", "Opening", "Closing"),
    target_n = ifelse(is.na(target_n), 0, target_n),
    target_pct = ifelse(is.na(target_pct), 0, target_pct),
    bg_pct = ifelse(is.na(bg_pct), 0, bg_pct),
    p_value = ifelse(is.na(p_value), 1, p_value),
    q_value = ifelse(is.na(q_value), 1, q_value),
    log10_p = cap_log10(p_value, cap = 300),
    motif_label = short_motif(motif_raw)
  )

write.csv(
  df,
  file.path(plot_dir, "HOMER_knownResults_parsed_kidney.csv"),
  row.names = FALSE
)

# -----------------------------
# 4) Contrast order
# -----------------------------
contrast_levels <- c(
  "Day14_vs_Sham",
  "Day42_vs_Sham",
  "Day42_vs_Day14"
)

df$contrast <- factor(
  df$contrast,
  levels = intersect(
    contrast_levels,
    unique(as.character(df$contrast))
  )
)

# -----------------------------
# 5) Opening:
#    TRUE amalgamation
# -----------------------------
top_n_open <- 10
q_cut_open <- 0.05
max_union_per_contrast <- 25

open_df <- df %>%
  filter(direction == "Opening") %>%
  mutate(q_use = ifelse(is.na(q_value), 1, q_value))

# within each cell type + contrast + motif, keep best-ranked instance
open_df_u <- open_df %>%
  group_by(cell_type, contrast, motif_raw) %>%
  slice_min(rank, n = 1, with_ties = FALSE) %>%
  ungroup()

# each cell type contributes its own top10 significant motifs
open_top_raw <- open_df_u %>%
  filter(q_use < q_cut_open) %>%
  arrange(rank) %>%
  group_by(cell_type, contrast) %>%
  slice_head(n = top_n_open) %>%
  ungroup()

# TRUE amalgamation:
# take UNION across cell types within each contrast
# no n_ct >= 2 filtering
open_top <- open_top_raw %>%
  group_by(contrast, motif_raw) %>%
  summarise(
    motif_label = dplyr::first(motif_label),
    n_ct = dplyr::n_distinct(cell_type),
    mean_log10 = mean(log10_p, na.rm = TRUE),
    best_log10 = max(log10_p, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(contrast) %>%
  arrange(
    desc(n_ct),
    desc(mean_log10),
    desc(best_log10),
    .by_group = TRUE
  ) %>%
  slice_head(n = max_union_per_contrast) %>%
  ungroup()

# only show a point if that motif was truly in that cell type's own top10
open_plot <- open_top_raw %>%
  inner_join(
    open_top,
    by = c("contrast", "motif_raw", "motif_label")
  ) %>%
  group_by(cell_type, contrast, motif_raw, motif_label) %>%
  summarise(
    log10_p = max(log10_p, na.rm = TRUE),
    target_pct = max(target_pct, na.rm = TRUE),
    target_n = max(target_n, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    motif_facet = paste(
      motif_label,
      motif_raw,
      as.character(contrast),
      sep = "___"
    )
  )

open_levels <- open_top %>%
  mutate(
    motif_facet = paste(
      motif_label,
      motif_raw,
      as.character(contrast),
      sep = "___"
    )
  ) %>%
  arrange(
    contrast,
    desc(n_ct),
    desc(mean_log10),
    desc(best_log10)
  ) %>%
  pull(motif_facet)

open_plot$motif_facet <- factor(
  open_plot$motif_facet,
  levels = unique(open_levels)
)

p_open <- ggplot(
  open_plot,
  aes(
    x = cell_type,
    y = motif_facet,
    size = target_pct,
    colour = log10_p
  )
) +
  geom_point(alpha = 0.9) +
  facet_wrap(~ contrast, nrow = 1, scales = "free_y") +
  scale_y_discrete(
    labels = function(x) sub("___.*$", "", x)
  ) +
  scale_size(
    range = c(2, 10),
    name = "Target %"
  ) +
  scale_colour_gradient(
    low = "#FDD0D0",
    high = "#B30000",
    limits = c(0, 50),
    oob = scales::squish,
    name = "-log10(P)"
  ) +
  labs(
    title = paste0(
      "Kidney Opening DAR motifs ",
      "(true amalgamation of top motifs)"
    ),
    x = "Cell type",
    y = "Motif"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )

# save opening
ggsave(
  file.path(plot_dir, "Kidney_HOMER_opening_true_amalgamation.pdf"),
  p_open,
  width = 16,
  height = 8
)

ggsave(
  file.path(plot_dir, "Kidney_HOMER_opening_true_amalgamation.png"),
  p_open,
  width = 16,
  height = 8,
  dpi = 300
)

write.csv(
  open_top,
  file.path(plot_dir, "Kidney_opening_true_amalgamation_selected_motifs.csv"),
  row.names = FALSE
)

write.csv(
  open_plot,
  file.path(plot_dir, "Kidney_opening_true_amalgamation_plot_data.csv"),
  row.names = FALSE
)

# -----------------------------
# 6) Closing:
#    top 25 per individual cell type
# -----------------------------
top_n_cell <- 25
q_cut_cell <- 0.05

build_celltype_plot_data <- function(df_in, dir_use,
                                     top_n = 25,
                                     q_cut = 0.05) {
  sub_df <- df_in %>%
    filter(direction == dir_use) %>%
    mutate(q_use = ifelse(is.na(q_value), 1, q_value))
  
  sub_u <- sub_df %>%
    group_by(cell_type, contrast, motif_raw) %>%
    slice_min(rank, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  sub_top_raw <- sub_u %>%
    filter(q_use < q_cut) %>%
    arrange(rank) %>%
    group_by(cell_type, contrast) %>%
    slice_head(n = top_n) %>%
    ungroup()
  
  plot_df <- sub_top_raw %>%
    group_by(cell_type, contrast, motif_raw, motif_label) %>%
    summarise(
      log10_p = max(log10_p, na.rm = TRUE),
      target_pct = max(target_pct, na.rm = TRUE),
      target_n = max(target_n, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      motif_facet = paste(
        motif_label,
        motif_raw,
        cell_type,
        sep = "___"
      )
    )
  
  y_levels <- plot_df %>%
    group_by(cell_type, motif_facet) %>%
    summarise(
      best = max(log10_p, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(cell_type, desc(best)) %>%
    pull(motif_facet)
  
  plot_df$motif_facet <- factor(
    plot_df$motif_facet,
    levels = unique(y_levels)
  )
  
  plot_df
}

plot_celltype_motifs <- function(plot_df, dir_label) {
  ggplot(
    plot_df,
    aes(
      x = contrast,
      y = motif_facet,
      size = target_pct,
      colour = log10_p
    )
  ) +
    geom_point(alpha = 0.85) +
    facet_wrap(~ cell_type, ncol = 2, scales = "free_y") +
    scale_y_discrete(
      labels = function(x) sub("___.*$", "", x)
    ) +
    scale_size(
      range = c(2, 10),
      name = "Target %"
    ) +
    scale_colour_gradient(
      low = "#D0E1F2",
      high = "#084594",
      limits = c(0, 50),
      oob = scales::squish,
      name = "-log10(P)"
    ) +
    labs(
      title = paste0(
        dir_label,
        " DAR motifs (top 25 per cell type)"
      ),
      x = "Contrast",
      y = "Motif"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1
      ),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold")
    )
}

close_cell_plot <- build_celltype_plot_data(
  df_in = df,
  dir_use = "Closing",
  top_n = top_n_cell,
  q_cut = q_cut_cell
)

p_close_cell <- plot_celltype_motifs(
  close_cell_plot,
  dir_label = "Closing"
)

ggsave(
  file.path(plot_dir, "Kidney_HOMER_closing_celltype_top25.pdf"),
  p_close_cell,
  width = 14,
  height = 14
)

ggsave(
  file.path(plot_dir, "Kidney_HOMER_closing_celltype_top25.png"),
  p_close_cell,
  width = 14,
  height = 14,
  dpi = 300
)

write.csv(
  close_cell_plot,
  file.path(plot_dir, "Kidney_closing_celltype_top25_plot_data.csv"),
  row.names = FALSE
)

# -----------------------------
# 7) Show in RStudio
# -----------------------------
p_open
p_close_cell
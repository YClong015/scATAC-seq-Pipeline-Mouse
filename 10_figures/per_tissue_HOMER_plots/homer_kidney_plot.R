#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

# ============================================================
# Kidney HOMER motif bubble plots (adapted from motif_visualisation.R)
#
# Input directory structure:
#   HOMER_bg_tested_nomotif/
#     <CellType>__<Contrast>__opening/knownResults.txt
#     <CellType>__<Contrast>__closing/knownResults.txt
#
# Outputs:
#   HOMER_Plots/Kidney_HOMER_opening_bubble.{pdf,png}
#   HOMER_Plots/Kidney_HOMER_closing_bubble.{pdf,png}
#   HOMER_Plots/HOMER_knownResults_parsed_kidney.csv
# ============================================================

# -----------------------------
# 1) Paths (edit if needed)
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
        any(grepl(paste(patterns, collapse = "|"), x,
                  ignore.case = TRUE))
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
    log10_p = cap_log10(p_value, cap = 300),
    motif_label = short_motif(motif_raw)
  )

write.csv(
  df,
  file.path(plot_dir, "HOMER_knownResults_parsed_kidney.csv"),
  row.names = FALSE
)

# Contrast order (Kidney)
contrast_levels <- c(
  "Day14_vs_Sham",
  "Day42_vs_Sham",
  "Day42_vs_Day14"
)

df$contrast <- factor(
  df$contrast,
  levels = intersect(contrast_levels, unique(as.character(df$contrast)))
)
# -----------------------------
# 4) Plot A: Opening amalgamated top motifs per contrast
# -----------------------------
top_n_open <- 10
min_ct <- 2
q_cut <- 0.05

open_df <- df %>%
  filter(direction == "Opening") %>%
  mutate(q_use = ifelse(is.na(q_value), 1, q_value))

# For each cell_type + contrast, keep the motif with the best rank (lowest rank number)
open_df_u <- open_df %>%
  group_by(cell_type, contrast, motif_raw) %>%
  slice_min(rank, n = 1, with_ties = FALSE) %>%
  ungroup()

# Filter to top motifs per contrast, prioritizing those that appear in more cell types
open_top_raw <- open_df_u %>%
  filter(q_use < q_cut) %>%
  arrange(rank) %>%
  group_by(cell_type, contrast) %>%
  slice_head(n = top_n_open) %>%
  ungroup()

# contrasr and motif_new combination, count how many cell types have it, and get mean and best log10_p
open_top <- open_top_raw %>%
  group_by(contrast, motif_raw) %>%
  summarise(
    motif_label = dplyr::first(motif_label),
    n_ct = dplyr::n_distinct(cell_type),
    mean_log10 = mean(log10_p, na.rm = TRUE),
    best_log10 = max(log10_p, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_ct >= min_ct) %>%
  group_by(contrast) %>%
  arrange(
    desc(n_ct),
    desc(mean_log10),
    desc(best_log10),
    .by_group = TRUE
  ) %>%
  slice_head(n = top_n_open) %>%
  ungroup() %>%
  select(contrast, motif_raw, motif_label)

# Only keep motifs that are in the top list for each contrast, but keep all cell types that have it in their top list
open_keep <- open_top_raw %>%
  select(cell_type, contrast, motif_raw) %>%
  distinct()

open_plot <- open_df_u %>%
  inner_join(
    open_top,
    by = c("contrast", "motif_raw"),
    suffix = c(".df", ".top")
  ) %>%
  inner_join(
    open_keep,
    by = c("cell_type", "contrast", "motif_raw")
  ) %>%
  mutate(
    motif_label = dplyr::coalesce(motif_label.top, motif_label.df)
  ) %>%
  group_by(cell_type, contrast, motif_raw, motif_label) %>%
  summarise(
    log10_p = max(log10_p, na.rm = TRUE),
    target_pct = max(target_pct, na.rm = TRUE),
    target_n = max(target_n, na.rm = TRUE),
    .groups = "drop"
  )

# For plotting, create a combined motif-contrast facet variable, and order by best log10_p across cell types
open_plot <- open_plot %>%
  mutate(
    motif_facet = paste(
      motif_label,
      motif_raw,
      as.character(contrast),
      sep = "___"
    )
  )

open_levels <- open_plot %>%
  group_by(contrast, motif_facet) %>%
  summarise(
    best = max(log10_p, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(contrast, desc(best)) %>%
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
  geom_point(alpha = 0.85) +
  facet_wrap(~ contrast, nrow = 1, scales = "free_y") +
  scale_y_discrete(
    labels = function(x) sub("___.*$", "", x)
  ) +
  scale_size(
    range = c(2, 10),
    name = "Target %"
  ) +
  scale_colour_gradient(
    low = "grey85",
    high = "red",
    limits = c(0, 50),
    oob = scales::squish,
    name = "-log10(P)"
  ) +
  labs(
    title = "A  Opening DAR motifs (amalgamated across cell types)",
    x = "Cell type",
    y = "Motif"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
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

p_open


# -----------------------------
# Cell-type-specific top 25 motifs
# -----------------------------
top_n_cell <- 25
min_contrast <- 1
q_cut_cell <- 0.05

build_celltype_plot_data <- function(df_in, dir_use,
                                     top_n = 25,
                                     min_con = 1,
                                     q_cut = 0.05) {
  sub_df <- df_in %>%
    dplyr::filter(.data$direction == dir_use) %>%
    dplyr::mutate(
      q_use = ifelse(is.na(.data$q_value), 1, .data$q_value)
    )
  
  sub_u <- sub_df %>%
    dplyr::group_by(
      .data$cell_type,
      .data$contrast,
      .data$motif_raw
    ) %>%
    dplyr::slice_min(.data$rank, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()
  
  sub_top_raw <- sub_u %>%
    dplyr::filter(.data$q_use < q_cut) %>%
    dplyr::arrange(.data$rank) %>%
    dplyr::group_by(.data$cell_type, .data$contrast) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::ungroup()
  
  sub_top <- sub_top_raw %>%
    dplyr::group_by(.data$cell_type, .data$motif_raw) %>%
    dplyr::summarise(
      motif_label = dplyr::first(.data$motif_label),
      n_contrast = dplyr::n_distinct(.data$contrast),
      mean_log10 = mean(.data$log10_p, na.rm = TRUE),
      best_log10 = max(.data$log10_p, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::filter(.data$n_contrast >= min_con) %>%
    dplyr::group_by(.data$cell_type) %>%
    dplyr::arrange(
      dplyr::desc(.data$n_contrast),
      dplyr::desc(.data$mean_log10),
      dplyr::desc(.data$best_log10),
      .by_group = TRUE
    ) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::ungroup() %>%
    dplyr::select(
      .data$cell_type,
      .data$motif_raw,
      .data$motif_label
    )
  
  sub_keep <- sub_top_raw %>%
    dplyr::select(
      .data$cell_type,
      .data$contrast,
      .data$motif_raw
    ) %>%
    dplyr::distinct()
  
  plot_df <- sub_u %>%
    dplyr::inner_join(
      sub_top,
      by = c("cell_type", "motif_raw"),
      suffix = c(".df", ".top")
    ) %>%
    dplyr::inner_join(
      sub_keep,
      by = c("cell_type", "contrast", "motif_raw")
    ) %>%
    dplyr::mutate(
      motif_label = dplyr::coalesce(
        .data$motif_label.top,
        .data$motif_label.df
      )
    ) %>%
    dplyr::group_by(
      .data$cell_type,
      .data$contrast,
      .data$motif_raw,
      .data$motif_label
    ) %>%
    dplyr::summarise(
      log10_p = max(.data$log10_p, na.rm = TRUE),
      target_pct = max(.data$target_pct, na.rm = TRUE),
      target_n = max(.data$target_n, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      motif_facet = paste(
        .data$motif_label,
        .data$motif_raw,
        .data$cell_type,
        sep = "___"
      )
    )
  
  y_levels <- plot_df %>%
    dplyr::group_by(.data$cell_type, .data$motif_facet) %>%
    dplyr::summarise(
      best = max(.data$log10_p, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(
      .data$cell_type,
      dplyr::desc(.data$best)
    ) %>%
    dplyr::pull(.data$motif_facet)
  
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
      low = "grey85",
      high = ifelse(dir_label == "Opening", "red", "blue"),
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
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      ),
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

# Build plot data
open_cell_plot <- build_celltype_plot_data(
  df_in = df,
  dir_use = "Opening",
  top_n = top_n_cell,
  min_con = min_contrast,
  q_cut = q_cut_cell
)

close_cell_plot <- build_celltype_plot_data(
  df_in = df,
  dir_use = "Closing",
  top_n = top_n_cell,
  min_con = min_contrast,
  q_cut = q_cut_cell
)

# Draw
p_open_cell <- plot_celltype_motifs(
  open_cell_plot,
  dir_label = "Opening"
)

p_close_cell <- plot_celltype_motifs(
  close_cell_plot,
  dir_label = "Closing"
)

# Save
ggsave(
  file.path(plot_dir, "Kidney_HOMER_opening_celltype_top25.pdf"),
  p_open_cell,
  width = 14,
  height = 14
)

ggsave(
  file.path(plot_dir, "Kidney_HOMER_opening_celltype_top25.png"),
  p_open_cell,
  width = 14,
  height = 14,
  dpi = 300
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
# -----------------------------
# One plot per cell type
# -----------------------------
plot_one_celltype_motifs <- function(plot_df, cell_use, dir_label) {
  sub_df <- plot_df %>%
    dplyr::filter(.data$cell_type == cell_use)
  
  ggplot(
    sub_df,
    aes(
      x = contrast,
      y = motif_facet,
      size = target_pct,
      colour = log10_p
    )
  ) +
    geom_point(alpha = 0.85) +
    scale_y_discrete(
      labels = function(x) sub("___.*$", "", x)
    ) +
    scale_size(
      range = c(2, 10),
      name = "Target %"
    ) +
    scale_colour_gradient(
      low = "grey85",
      high = ifelse(dir_label == "Opening", "red", "blue"),
      limits = c(0, 50),
      oob = scales::squish,
      name = "-log10(P)"
    ) +
    labs(
      title = paste0(cell_use, " - ", dir_label, " DAR motifs"),
      x = "Contrast",
      y = "Motif"
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      ),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1
      ),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )
}

# Opening
cell_types_open <- unique(open_cell_plot$cell_type)

for (ct in cell_types_open) {
  p_ct <- plot_one_celltype_motifs(
    plot_df = open_cell_plot,
    cell_use = ct,
    dir_label = "Opening"
  )
  
  ggsave(
    file.path(
      plot_dir,
      paste0("Kidney_HOMER_opening_", gsub(" ", "_", ct), ".pdf")
    ),
    p_ct,
    width = 6,
    height = 8
  )
  
  ggsave(
    file.path(
      plot_dir,
      paste0("Kidney_HOMER_opening_", gsub(" ", "_", ct), ".png")
    ),
    p_ct,
    width = 6,
    height = 8,
    dpi = 300
  )
}

# Closing
cell_types_close <- unique(close_cell_plot$cell_type)

for (ct in cell_types_close) {
  p_ct <- plot_one_celltype_motifs(
    plot_df = close_cell_plot,
    cell_use = ct,
    dir_label = "Closing"
  )
  
  ggsave(
    file.path(
      plot_dir,
      paste0("Kidney_HOMER_closing_", gsub(" ", "_", ct), ".pdf")
    ),
    p_ct,
    width = 6,
    height = 8
  )
  
  ggsave(
    file.path(
      plot_dir,
      paste0("Kidney_HOMER_closing_", gsub(" ", "_", ct), ".png")
    ),
    p_ct,
    width = 6,
    height = 8,
    dpi = 300
  )
}

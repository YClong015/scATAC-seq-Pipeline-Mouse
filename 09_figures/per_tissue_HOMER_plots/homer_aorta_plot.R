#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

# ============================================================
# Aorta HOMER motif plots
#
# Opening:
#   top 10 motifs across cell types
#
# Closing:
#   top 25 motifs for individual cell types
# ============================================================

# -----------------------------
# 1) Paths
# -----------------------------
out_dir <- paste0(
  "/QRISdata/Q8448/Mouse_disease_data/DAR/",
  "DAR_pseudobulk_DESeq2/DAR_pseudobulk_Aorta_DESeq2"
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
  file.path(plot_dir, "HOMER_knownResults_parsed_aorta.csv"),
  row.names = FALSE
)

# Keep only expected contrast
df <- df %>%
  filter(contrast == "Challenge_vs_Control")

celltype_levels <- c("Macrophages", "Pericytes", "SMC")
df$cell_type <- factor(
  df$cell_type,
  levels = intersect(celltype_levels,
                     unique(as.character(df$cell_type)))
)

# -----------------------------
# 4) Opening:
#    top 10 across cell types
# -----------------------------
build_opening_plot_data <- function(df_in,
                                    top_n_each = 10,
                                    top_n_final = 10,
                                    q_cut = 0.05) {
  sub_df <- df_in %>%
    filter(direction == "Opening") %>%
    mutate(q_use = ifelse(is.na(q_value), 1, q_value))
  
  sub_u <- sub_df %>%
    group_by(cell_type, motif_raw) %>%
    slice_min(rank, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  sub_top_raw <- sub_u %>%
    filter(q_use < q_cut) %>%
    arrange(rank) %>%
    group_by(cell_type) %>%
    slice_head(n = top_n_each) %>%
    ungroup()
  
  sub_top <- sub_top_raw %>%
    group_by(motif_raw) %>%
    summarise(
      motif_label = dplyr::first(motif_label),
      n_ct = dplyr::n_distinct(cell_type),
      mean_log10 = mean(log10_p, na.rm = TRUE),
      best_log10 = max(log10_p, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(n_ct), desc(mean_log10), desc(best_log10)) %>%
    slice_head(n = top_n_final) %>%
    select(motif_raw, motif_label)
  
  sub_keep <- sub_top_raw %>%
    select(cell_type, motif_raw) %>%
    distinct()
  
  plot_df <- sub_u %>%
    inner_join(sub_top, by = "motif_raw",
               suffix = c(".df", ".top")) %>%
    inner_join(sub_keep, by = c("cell_type", "motif_raw")) %>%
    mutate(
      motif_label = dplyr::coalesce(motif_label.top, motif_label.df)
    ) %>%
    group_by(cell_type, motif_raw, motif_label) %>%
    summarise(
      log10_p = max(log10_p, na.rm = TRUE),
      target_pct = max(target_pct, na.rm = TRUE),
      target_n = max(target_n, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      motif_facet = paste(motif_label, motif_raw, sep = "___")
    )
  
  motif_levels <- plot_df %>%
    group_by(motif_facet) %>%
    summarise(
      n_ct = dplyr::n_distinct(cell_type),
      best = max(log10_p, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(n_ct), desc(best)) %>%
    pull(motif_facet)
  
  plot_df$motif_facet <- factor(
    plot_df$motif_facet,
    levels = unique(motif_levels)
  )
  
  plot_df
}

plot_opening <- function(plot_df) {
  ggplot(
    plot_df,
    aes(
      x = cell_type,
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
      high = "red",
      limits = c(0, 50),
      oob = scales::squish,
      name = "-log10(P)"
    ) +
    labs(
      title = "Aorta Opening DAR motifs (top 10 across cell types)",
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
      panel.grid.minor = element_blank()
    )
}

# -----------------------------
# 5) Closing:
#    top 25 for individual cell types
# -----------------------------
build_closing_plot_data <- function(df_in,
                                    top_n_each = 25,
                                    q_cut = 0.05) {
  sub_df <- df_in %>%
    filter(direction == "Closing") %>%
    mutate(q_use = ifelse(is.na(q_value), 1, q_value))
  
  sub_u <- sub_df %>%
    group_by(cell_type, motif_raw) %>%
    slice_min(rank, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  sub_top_raw <- sub_u %>%
    filter(q_use < q_cut) %>%
    arrange(rank) %>%
    group_by(cell_type) %>%
    slice_head(n = top_n_each) %>%
    ungroup()
  
  plot_df <- sub_top_raw %>%
    group_by(cell_type, motif_raw, motif_label) %>%
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
  
  motif_levels <- plot_df %>%
    group_by(cell_type, motif_facet) %>%
    summarise(
      best = max(log10_p, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(cell_type, desc(best)) %>%
    pull(motif_facet)
  
  plot_df$motif_facet <- factor(
    plot_df$motif_facet,
    levels = unique(motif_levels)
  )
  
  plot_df
}

plot_closing_facet <- function(plot_df) {
  ggplot(
    plot_df,
    aes(
      x = cell_type,
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
      high = "blue",
      limits = c(0, 50),
      oob = scales::squish,
      name = "-log10(P)"
    ) +
    labs(
      title = paste0(
        "Aorta Closing DAR motifs ",
        "(top 25 for individual cell types)"
      ),
      x = "Cell type",
      y = "Motif"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold")
    )
}

plot_one_closing_celltype <- function(plot_df, cell_use) {
  sub_df <- plot_df %>%
    filter(cell_type == cell_use)
  
  ggplot(
    sub_df,
    aes(
      x = cell_type,
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
      high = "blue",
      limits = c(0, 50),
      oob = scales::squish,
      name = "-log10(P)"
    ) +
    labs(
      title = paste0(
        "Aorta Closing DAR motifs - ",
        cell_use,
        " (top 25)"
      ),
      x = "Cell type",
      y = "Motif"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )
}

# -----------------------------
# 6) Build data and draw
# -----------------------------
open_plot <- build_opening_plot_data(
  df_in = df,
  top_n_each = 10,
  top_n_final = 10,
  q_cut = 0.05
)

close_plot <- build_closing_plot_data(
  df_in = df,
  top_n_each = 25,
  q_cut = 0.05
)

write.csv(
  open_plot,
  file.path(plot_dir, "Aorta_opening_plot_data_top10.csv"),
  row.names = FALSE
)

write.csv(
  close_plot,
  file.path(plot_dir, "Aorta_closing_plot_data_top25.csv"),
  row.names = FALSE
)

p_open <- plot_opening(open_plot)
p_close_facet <- plot_closing_facet(close_plot)

p_open
p_close_facet

# -----------------------------
# 7) Save opening
# -----------------------------
ggsave(
  file.path(plot_dir, "Aorta_HOMER_opening_top10.pdf"),
  p_open,
  width = 9,
  height = 6
)

ggsave(
  file.path(plot_dir, "Aorta_HOMER_opening_top10.png"),
  p_open,
  width = 9,
  height = 6,
  dpi = 300
)

# -----------------------------
# 8) Save closing facet
# -----------------------------
ggsave(
  file.path(plot_dir, "Aorta_HOMER_closing_top25_facet.pdf"),
  p_close_facet,
  width = 12,
  height = 12
)

ggsave(
  file.path(plot_dir, "Aorta_HOMER_closing_top25_facet.png"),
  p_close_facet,
  width = 12,
  height = 12,
  dpi = 300
)

# -----------------------------
# 9) Save closing per cell type
# -----------------------------
cell_types_close <- unique(as.character(close_plot$cell_type))

for (ct in cell_types_close) {
  p_ct <- plot_one_closing_celltype(close_plot, ct)
  
  ggsave(
    file.path(
      plot_dir,
      paste0("Aorta_HOMER_closing_", gsub(" ", "_", ct), "_top25.pdf")
    ),
    p_ct,
    width = 7,
    height = 8
  )
  
  ggsave(
    file.path(
      plot_dir,
      paste0("Aorta_HOMER_closing_", gsub(" ", "_", ct), "_top25.png")
    ),
    p_ct,
    width = 7,
    height = 8,
    dpi = 300
  )
}


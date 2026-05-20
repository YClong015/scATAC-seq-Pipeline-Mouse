#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

# ============================================================
# Tcell HOMER motif plots
#
# New logic:
# - Opening: top 10 motifs per contrast
# - Closing: top 25 motifs per contrast
# - Plot per contrast, not amalgamated across contrasts
# ============================================================

# -----------------------------
# 1) Paths
# -----------------------------
out_dir <- paste0(
  "/QRISdata/Q8448/Mouse_disease_data/DAR/",
  "DAR_pseudobulk_DESeq2/DAR_pseudobulk_Tcells_DESeq2"
)

homer_dir <- file.path(out_dir, "HOMER")
plot_dir <- file.path(out_dir, "HOMER_Plots")

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 2) Helpers
# -----------------------------
parse_dirname <- function(dname) {
  parts <- strsplit(dname, "__", fixed = TRUE)[[1]]
  
  if (length(parts) < 3) {
    return(NULL)
  }
  
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
  
  if (length(hit) == 0) {
    return(NA_character_)
  }
  
  hit[1]
}

short_motif <- function(x) {
  x <- as.character(x)
  sapply(strsplit(x, "/"), `[`, 1)
}

read_known <- function(fp) {
  dname <- basename(dirname(fp))
  meta <- parse_dirname(dname)
  
  if (is.null(meta)) {
    return(NULL)
  }
  
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
  
  if (is.null(tab) || nrow(tab) < 1) {
    return(NULL)
  }
  
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
  
  if (is.na(motif_col) || is.na(p_col)) {
    return(NULL)
  }
  
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

# -----------------------------
# 3) Read HOMER results
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
  file.path(plot_dir, "HOMER_knownResults_parsed_tcell.csv"),
  row.names = FALSE
)

# -----------------------------
# 4) Contrast order
# -----------------------------
contrast_levels <- c(
  "Aged_vs_Juvenile",
  "Aged_vs_Young_acute",
  "Aged_vs_Young_chronic",
  "Aged_vs_Young_control",
  "Juvenile_vs_Young_acute",
  "Juvenile_vs_Young_chronic",
  "Juvenile_vs_Young_control",
  "Young_acute_vs_Young_control",
  "Young_chronic_vs_Young_control"
)

df$contrast <- factor(
  df$contrast,
  levels = intersect(
    contrast_levels,
    unique(as.character(df$contrast))
  )
)

# -----------------------------
# 5) Build per-contrast top data
# -----------------------------
build_per_contrast_plot_data <- function(df_in, dir_use,
                                         top_n = 10,
                                         q_cut = 0.05) {
  sub_df <- df_in %>%
    filter(direction == dir_use) %>%
    mutate(q_use = ifelse(is.na(q_value), 1, q_value))
  
  sub_u <- sub_df %>%
    group_by(contrast, motif_raw) %>%
    slice_min(rank, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  plot_df <- sub_u %>%
    filter(q_use < q_cut) %>%
    arrange(rank) %>%
    group_by(contrast) %>%
    slice_head(n = top_n) %>%
    ungroup() %>%
    mutate(
      motif_facet = paste(
        motif_label,
        motif_raw,
        as.character(contrast),
        sep = "___"
      )
    )
  
  y_levels <- plot_df %>%
    group_by(contrast, motif_facet) %>%
    summarise(
      best = max(log10_p, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(contrast, best) %>%
    pull(motif_facet)
  
  plot_df$motif_facet <- factor(
    plot_df$motif_facet,
    levels = unique(y_levels)
  )
  
  plot_df
}

# -----------------------------
# 6) Plot functions
# -----------------------------
plot_facet_motifs <- function(plot_df, dir_label, top_n,
                              high_col) {
  ggplot(
    plot_df,
    aes(
      x = log10_p,
      y = motif_facet,
      size = target_pct,
      colour = log10_p
    )
  ) +
    geom_segment(
      aes(
        x = 0,
        xend = log10_p,
        y = motif_facet,
        yend = motif_facet
      ),
      colour = "grey80",
      linewidth = 0.4,
      inherit.aes = FALSE
    ) +
    geom_point(alpha = 0.9) +
    facet_wrap(~ contrast, ncol = 3, scales = "free_y") +
    scale_y_discrete(
      labels = function(x) sub("___.*$", "", x)
    ) +
    scale_size(
      range = c(2, 10),
      name = "Target %"
    ) +
    scale_colour_gradient(
      low = "grey85",
      high = high_col,
      limits = c(0, 50),
      oob = scales::squish,
      name = "-log10(P)"
    ) +
    labs(
      title = paste0(
        "Tcell ",
        dir_label,
        " DAR motifs (top ",
        top_n,
        " per contrast)"
      ),
      x = "-log10(P)",
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

plot_one_contrast <- function(plot_df, contrast_use,
                              dir_label, high_col) {
  sub_df <- plot_df %>%
    filter(as.character(contrast) == contrast_use)
  
  ggplot(
    sub_df,
    aes(
      x = log10_p,
      y = motif_facet,
      size = target_pct,
      colour = log10_p
    )
  ) +
    geom_segment(
      aes(
        x = 0,
        xend = log10_p,
        y = motif_facet,
        yend = motif_facet
      ),
      colour = "grey80",
      linewidth = 0.4,
      inherit.aes = FALSE
    ) +
    geom_point(alpha = 0.9) +
    scale_y_discrete(
      labels = function(x) sub("___.*$", "", x)
    ) +
    scale_size(
      range = c(2, 10),
      name = "Target %"
    ) +
    scale_colour_gradient(
      low = "grey85",
      high = high_col,
      limits = c(0, 50),
      oob = scales::squish,
      name = "-log10(P)"
    ) +
    labs(
      title = paste0(
        "Tcell ",
        dir_label,
        " DAR motifs - ",
        contrast_use
      ),
      x = "-log10(P)",
      y = "Motif"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      ),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )
}

# -----------------------------
# 7) Build data
# -----------------------------
open_plot <- build_per_contrast_plot_data(
  df_in = df,
  dir_use = "Opening",
  top_n = 10,
  q_cut = 0.05
)

close_plot <- build_per_contrast_plot_data(
  df_in = df,
  dir_use = "Closing",
  top_n = 25,
  q_cut = 0.05
)

write.csv(
  open_plot,
  file.path(plot_dir, "Tcell_opening_plot_data_top10_by_contrast.csv"),
  row.names = FALSE
)

write.csv(
  close_plot,
  file.path(plot_dir, "Tcell_closing_plot_data_top25_by_contrast.csv"),
  row.names = FALSE
)

# -----------------------------
# 8) Facet plots
# -----------------------------
p_open <- plot_facet_motifs(
  open_plot,
  dir_label = "Opening",
  top_n = 10,
  high_col = "red"
)

p_close <- plot_facet_motifs(
  close_plot,
  dir_label = "Closing",
  top_n = 25,
  high_col = "blue"
)

p_open
p_close

ggsave(
  file.path(
    plot_dir,
    "Tcell_HOMER_opening_top10_by_contrast_facet.pdf"
  ),
  p_open,
  width = 14,
  height = 10
)

ggsave(
  file.path(
    plot_dir,
    "Tcell_HOMER_opening_top10_by_contrast_facet.png"
  ),
  p_open,
  width = 14,
  height = 10,
  dpi = 300
)

ggsave(
  file.path(
    plot_dir,
    "Tcell_HOMER_closing_top25_by_contrast_facet.pdf"
  ),
  p_close,
  width = 16,
  height = 14
)

ggsave(
  file.path(
    plot_dir,
    "Tcell_HOMER_closing_top25_by_contrast_facet.png"
  ),
  p_close,
  width = 16,
  height = 14,
  dpi = 300
)

# -----------------------------
# 9) Individual contrast plots
# -----------------------------
contrast_vec_open <- unique(as.character(open_plot$contrast))
contrast_vec_close <- unique(as.character(close_plot$contrast))

for (cc in contrast_vec_open) {
  p_cc <- plot_one_contrast(
    plot_df = open_plot,
    contrast_use = cc,
    dir_label = "Opening",
    high_col = "red"
  )
  
  ggsave(
    file.path(
      plot_dir,
      paste0(
        "Tcell_HOMER_opening_",
        gsub(" ", "_", cc),
        "_top10.pdf"
      )
    ),
    p_cc,
    width = 7,
    height = 6
  )
  
  ggsave(
    file.path(
      plot_dir,
      paste0(
        "Tcell_HOMER_opening_",
        gsub(" ", "_", cc),
        "_top10.png"
      )
    ),
    p_cc,
    width = 7,
    height = 6,
    dpi = 300
  )
}

for (cc in contrast_vec_close) {
  p_cc <- plot_one_contrast(
    plot_df = close_plot,
    contrast_use = cc,
    dir_label = "Closing",
    high_col = "blue"
  )
  
  ggsave(
    file.path(
      plot_dir,
      paste0(
        "Tcell_HOMER_closing_",
        gsub(" ", "_", cc),
        "_top25.pdf"
      )
    ),
    p_cc,
    width = 8,
    height = 10
  )
  
  ggsave(
    file.path(
      plot_dir,
      paste0(
        "Tcell_HOMER_closing_",
        gsub(" ", "_", cc),
        "_top25.png"
      )
    ),
    p_cc,
    width = 8,
    height = 10,
    dpi = 300
  )
}


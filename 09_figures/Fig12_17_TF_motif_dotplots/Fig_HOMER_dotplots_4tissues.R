#!/usr/bin/env Rscript
# ============================================================
# HOMER Dot Plots — 4 tissues, opening + closing
# Based on teacher's template, adapted for project data
#
# Directory naming convention:
#   Kidney / Lung / Aorta : {CellType}__{contrast}__{direction}
#   Tcells                 : Tcell__{contrast}__{direction}
#
# X axis:
#   Kidney / Lung / Aorta : cell type  (fixed contrast per tissue)
#   Tcells                 : contrast   (fixed cell type = Tcell)
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(scales)
})

# ---------------------------------------------------------------
# 1) Paths
# ---------------------------------------------------------------
base_dar <- "/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_pseudobulk_DESeq2"

out_dir <- file.path(base_dar, "HOMER_dotplots_4tissues")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------
# 2) Tissue configurations
#    comparison.set : names  = X-axis display label
#                     values = directory stem (without __direction)
# ---------------------------------------------------------------
tissue_cfg <- list(

  # HOMER_NS_bg directories: {CellType}__{contrast}__005__{opening_vs_NS|closing_vs_NS}
  # directions to use: "opening_vs_NS" and "closing_vs_NS"
  Kidney = list(
    homer_base     = file.path(base_dar, "DAR_pseudobulk_Kidney_DESeq2", "HOMER_NS_bg"),
    comparison.set = c(
      "DCT"         = "DCT__Day42_vs_Sham__005",
      "Endothelial" = "Endothelial__Day42_vs_Sham__005",
      "IC"          = "IC__Day42_vs_Sham__005",
      "Macrophages" = "Macrophages__Day42_vs_Sham__005",
      "PC"          = "PC__Day42_vs_Sham__005",
      "PT"          = "PT__Day42_vs_Sham__005",
      "TAL"         = "TAL__Day42_vs_Sham__005"
    ),
    directions    = c("opening_vs_NS", "closing_vs_NS"),
    total_motifs  = 10,
    force_motifs  = "Bach2(bZIP)",   # partial match, always included
    x_label       = "Cell Type"
  ),

  Lung = list(
    homer_base      = file.path(base_dar, "DAR_pseudobulk_Lung_DESeq2", "HOMER_NS_bg"),
    comparison.set  = NULL,          # auto-detected at runtime
    contrast_filter = "Case_vs_Control",
    directions      = c("opening_vs_NS", "closing_vs_NS"),
    x_label         = "Cell Type"
  ),

  Aorta = list(
    homer_base     = file.path(base_dar, "DAR_pseudobulk_Aorta_DESeq2", "HOMER_NS_bg"),
    comparison.set = c(
      "Macrophages" = "Macrophages__Challenge_vs_Control__exp005",
      "Pericytes"   = "Pericytes__Challenge_vs_Control__exp005",
      "SMC"         = "SMC__Challenge_vs_Control__exp005"
    ),
    directions = c("opening_vs_NS", "closing_vs_NS"),
    x_label    = "Cell Type"
  ),

  # HOMER/ directories: Tcell__{contrast}__{opening|closing}
  Tcells = list(
    homer_base     = file.path(base_dar, "DAR_pseudobulk_Tcells_DESeq2", "HOMER"),
    comparison.set = c(
      "Young_chronic" = "Tcell__Young_chronic_vs_Young_control",
      "Aged"          = "Tcell__Aged_vs_Young_control",
      "Young_acute"   = "Tcell__Young_acute_vs_Young_control"
    ),
    directions = c("opening", "closing"),
    x_label    = "Contrast"
  )
)

# ---------------------------------------------------------------
# 3) Auto-detect comparison.set from directory names
#    Used for Lung (exact cell type names unknown)
# ---------------------------------------------------------------
auto_comparisons <- function(homer_base, direction, contrast_filter = NULL) {
  all_dirs <- list.dirs(homer_base, recursive = FALSE, full.names = FALSE)
  pattern  <- paste0("__", direction, "$")
  matched  <- all_dirs[grepl(pattern, all_dirs, perl = TRUE)]

  if (!is.null(contrast_filter))
    matched <- matched[grepl(contrast_filter, matched, fixed = TRUE)]

  if (length(matched) == 0) return(NULL)

  # Display name = cell type = everything before first __
  ct    <- sub("__.*$", "", matched)
  # Stem = full directory name minus the trailing __direction
  stems <- sub(paste0("__", direction, "$"), "", matched)
  setNames(stems, ct)
}
# NOTE: for HOMER_NS_bg dirs like {CellType}__{contrast}__005__{direction},
# auto_comparisons strips the trailing __{direction} leaving {CellType}__{contrast}__005

# ---------------------------------------------------------------
# 4) Read HOMER data and build ggData
# ---------------------------------------------------------------
make_ggdata <- function(homer_base, comparison.set, direction,
                        top_n_per_group = 10, total_motifs = 10,
                        force_motifs = NULL) {

  # Step 1: union top-N motifs across all comparisons
  # Also force-include any motifs matching force_motifs (partial match against full table)
  motifs.set <- c()
  for (comp in comparison.set) {
    fp <- file.path(homer_base, paste0(comp, "__", direction), "knownResults.txt")
    if (!file.exists(fp)) {
      message("    [SKIP] not found: ", basename(dirname(fp))); next
    }
    tab <- tryCatch(
      read.csv(fp, sep = "\t", header = TRUE, stringsAsFactors = FALSE),
      error = function(e) NULL
    )
    if (is.null(tab) || nrow(tab) == 0) next
    rownames(tab) <- make.unique(tab$Motif.Name)

    # Add top N
    top_rows <- rownames(tab)[seq_len(min(top_n_per_group, nrow(tab)))]
    motifs.set <- union(motifs.set, top_rows)

    # Force-include matching motifs from the full table (fixed string match)
    if (!is.null(force_motifs)) {
      forced_rows <- rownames(tab)[
        Reduce(`|`, lapply(force_motifs, function(p)
          grepl(p, rownames(tab), fixed = TRUE)))
      ]
      motifs.set <- union(motifs.set, forced_rows)
    }
  }

  if (length(motifs.set) == 0) return(NULL)

  # Step 2: read full table for each comparison
  ggData <- lapply(names(comparison.set), function(display_nm) {
    comp <- comparison.set[[display_nm]]
    fp   <- file.path(homer_base, paste0(comp, "__", direction), "knownResults.txt")
    if (!file.exists(fp)) return(NULL)

    tab <- tryCatch(
      read.csv(fp, sep = "\t", header = TRUE, stringsAsFactors = FALSE),
      error = function(e) NULL
    )
    if (is.null(tab) || nrow(tab) == 0) return(NULL)

    rownames(tab) <- make.unique(tab$Motif.Name)

    # Align to motifs.set (fill missing with NA)
    sub           <- tab[match(motifs.set, rownames(tab)), , drop = FALSE]
    rownames(sub) <- motifs.set

    pct_tgt <- as.numeric(gsub("%", "", sub$X..of.Target.Sequences.with.Motif))
    pct_bg  <- as.numeric(gsub("%", "", sub$X..of.Background.Sequences.with.Motif))
    fc      <- log2((pct_tgt + 0.001) / (pct_bg + 0.001) + 1)
    logp    <- -log10(pmax(sub$P.value, 1e-300))

    data.frame(
      Comparison = display_nm,
      Motif      = motifs.set,
      Target_pct = pct_tgt,
      FC         = fc,
      Num_hits   = sub[, 6],
      Log10_pval = logp,
      stringsAsFactors = FALSE
    )
  })

  ggData <- bind_rows(Filter(Negate(is.null), ggData))
  if (nrow(ggData) == 0) return(NULL)

  # Cap Inf
  finite_max <- max(ggData$Log10_pval[is.finite(ggData$Log10_pval)], na.rm = TRUE)
  ggData$Log10_pval[!is.finite(ggData$Log10_pval)] <- finite_max

  # Step 3: keep only global top total_motifs by max -log10(p) across groups
  # Force-include specific motifs (partial name match), then fill remaining slots
  ranked <- ggData %>%
    group_by(Motif) %>%
    summarise(best_p = max(Log10_pval, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(best_p))

  if (!is.null(force_motifs)) {
    forced <- ranked %>%
      filter(Reduce(`|`, lapply(force_motifs, function(p)
        grepl(p, Motif, fixed = TRUE)))) %>%
      pull(Motif)
    # top total_motifs from non-forced, then add forced on top (extra)
    remaining <- ranked %>%
      filter(!Motif %in% forced) %>%
      slice_head(n = total_motifs) %>%
      pull(Motif)
    top_motifs <- c(remaining, forced)   # forced appended at bottom of Y axis
  } else {
    top_motifs <- ranked %>% slice_head(n = total_motifs) %>% pull(Motif)
  }

  ggData     <- ggData %>% filter(Motif %in% top_motifs)
  motifs.set <- motifs.set[motifs.set %in% top_motifs]

  # Shorten motif names: keep only the part before the first "/"
  # Build a lookup from full name → short unique name, then apply to ggData
  shorten <- function(x) sapply(strsplit(as.character(x), "/"), `[`, 1)
  short_names           <- make.unique(shorten(motifs.set))
  names(short_names)    <- motifs.set          # full name → short name
  motifs.set            <- unname(short_names)
  ggData$Motif          <- short_names[as.character(ggData$Motif)]

  # Factor levels
  ggData$Comparison <- factor(ggData$Comparison, levels = names(comparison.set))
  ggData$Motif      <- factor(ggData$Motif, levels = rev(motifs.set))

  ggData
}

# ---------------------------------------------------------------
# 5) Plot function
# ---------------------------------------------------------------
make_dotplot <- function(ggData, title, x_label, colour_high = "red") {
  # Sort largest dots first so smaller dots render on top and stay visible
  ggData <- ggData %>% arrange(desc(Target_pct))

  ggplot(ggData, aes(x = Comparison,
                     y = reorder(Motif, Log10_pval, FUN = max),
                     size = Target_pct,
                     colour = Log10_pval)) +
    geom_point(alpha = 0.75) +
    ylab("") +
    xlab(x_label) +
    scale_size(
      range  = c(1, 10),
      name   = "Target %",
      limits = c(0, max(ggData$Target_pct, na.rm = TRUE))
    ) +
    scale_colour_gradientn(
      colours = c("#e3e1e1", adjustcolor(colour_high, alpha.f = 0.55), colour_high),
      values  = rescale(c(0, 10, 200), to = c(0, 1)),
      name    = "-log10(P)"
    ) +
    labs(title = title) +
    theme_classic(base_size = 15) +
    theme(
      plot.title   = element_text(face = "bold", hjust = 0.5, size = 14),
      axis.text.x  = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.text.y  = element_text(size = 9),
      legend.position = "right"
    )
}

# ---------------------------------------------------------------
# 6) Main loop — 4 tissues × 2 directions
# ---------------------------------------------------------------
all_plots <- list()   # collect all plots for interactive preview

for (tis in names(tissue_cfg)) {
  cfg <- tissue_cfg[[tis]]
  message("\n=== ", tis, " ===")

  for (direction in cfg$directions) {
    message("  Direction: ", direction)

    # Lung: auto-detect comparison.set
    comp_set <- cfg$comparison.set
    if (is.null(comp_set)) {
      comp_set <- auto_comparisons(cfg$homer_base, direction, cfg$contrast_filter)
      if (is.null(comp_set)) {
        message("    No directories found for ", direction, " — skipping"); next
      }
      message("    Auto-detected: ", paste(names(comp_set), collapse = ", "))
    }

    ggData <- make_ggdata(cfg$homer_base, comp_set, direction,
                          top_n_per_group = 10,
                          total_motifs    = if (!is.null(cfg$total_motifs)) cfg$total_motifs else 10,
                          force_motifs    = cfg$force_motifs)

    if (is.null(ggData) || nrow(ggData) == 0) {
      message("    No data — skipping"); next
    }

    # opening_vs_NS → red, closing_vs_NS → blue
    colour     <- if (grepl("opening", direction)) "red" else "blue"
    dir_label  <- if (grepl("opening", direction)) "opening" else "closing"
    title      <- paste0(tis, " — ", dir_label, " DARs (top 10 motifs per group)")
    file_label <- paste0(tis, "_", dir_label, "_dotplot")

    p <- make_dotplot(ggData, title = title,
                      x_label = cfg$x_label, colour_high = colour)

    all_plots[[file_label]] <- p   # store for preview: all_plots$Kidney_opening_dotplot

    prefix <- file.path(out_dir, file_label)
    ggsave(paste0(prefix, ".pdf"), p, width = 10, height = 7)
    ggsave(paste0(prefix, ".png"), p, width = 10, height = 7, dpi = 300)
    message("    Saved: ", basename(prefix))
  }
}

message("\nDone. Output -> ", out_dir)
message("Available plots: ", paste(names(all_plots), collapse = ", "))
message("Preview any plot with: all_plots[[\"Kidney_opening_dotplot\"]]")

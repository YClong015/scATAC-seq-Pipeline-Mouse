#!/usr/bin/env Rscript
# Fig 11 - Top opening-DAR TF motifs (NS background), one dotplot per tissue.
# x = cell type, y = top-N motifs, size = % targets, colour = -log10(P).

library(ggplot2)
library(dplyr)
library(scales)
library(patchwork)

OUT_DIR <- "/QRISdata/Q8448/Mouse_disease_data/DAR/Fig11_output"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
TOP_N       <- 10   # top motifs per cell type (union within a tissue)
MAX_MOTIFS  <- 10   # cap on motifs shown on the y-axis

## per-tissue config; Aorta dir names use "__exp005__", the others "__005__"
TISSUES <- list(
  Kidney = list(
    letter      = "A",
    title       = "Kidney",
    homer_root  = "/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_pseudobulk_Kidney_v5_DESeq2/HOMER_NS_bg",
    open_pat    = "__005__opening_vs_NS$",
    keep_contrast = c("Day42_vs_Sham"),
    ct_order    = c("PCT","PST","Injured_PT","TAL","DCT_CNT","PC_URO","LEUK","FIB")
  ),
  Lung = list(
    letter      = "B",
    title       = "Lung",
    homer_root  = "/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_pseudobulk_DESeq2/DAR_pseudobulk_Lung_DESeq2/HOMER_NS_bg",
    open_pat    = "__005__opening_vs_NS$",
    keep_contrast = c("Case_vs_Control"),
    ct_order    = NULL
  ),
  Aorta = list(
    letter      = "C",
    title       = "Aorta",
    homer_root  = "/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_pseudobulk_DESeq2/DAR_pseudobulk_Aorta_DESeq2/HOMER_NS_bg",
    open_pat    = "__005__opening_vs_NS$",
    keep_contrast = c("Challenge_vs_Control"),
    ct_order    = NULL
  ),
  Tcell = list(
    letter      = "D",
    title       = "T cells",
    homer_root  = "/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_pseudobulk_DESeq2/DAR_pseudobulk_Tcells_DESeq2/HOMER_NS_bg",
    open_pat    = "__005__opening_vs_NS$",
    # Each condition vs the Young_control baseline (one column per group).
    # Juvenile dropped (not relevant to age-related disease).
    keep_contrast = c("Young_acute_vs_Young_control",
                      "Young_chronic_vs_Young_control",
                      "Aged_vs_Young_control"),
    ct_order    = NULL,
    x_axis      = "Contrast",        # T cells: one pool, x-axis = group
    x_order     = c("Young acute", "Young chronic", "Aged")
  )
)

# Tissues without an explicit x_axis default to cell type on the x-axis.
for (nm in names(TISSUES))
  if (is.null(TISSUES[[nm]]$x_axis)) TISSUES[[nm]]$x_axis <- "CellType"

simplify.motif <- function(x) sub("/.*$", "", x)   # "Fra1(bZIP)/.../Homer" -> "Fra1(bZIP)"

## Discover (ct, contrast) units for one tissue
discover.units <- function(cfg) {
  if (!dir.exists(cfg$homer_root)) {
    message("  [SKIP] homer_root missing: ", cfg$homer_root); return(NULL)
  }
  bn   <- basename(list.dirs(cfg$homer_root, recursive = FALSE))
  hits <- bn[grepl(cfg$open_pat, bn)]
  if (length(hits) == 0) {
    message("  [SKIP] no dirs match ", cfg$open_pat, " in ", cfg$homer_root)
    return(NULL)
  }
  stripped <- sub(cfg$open_pat, "", hits)
  parts    <- strsplit(stripped, "__", fixed = TRUE)
  ct       <- vapply(parts, `[`, "", 1)
  contrast <- vapply(parts, function(p) if (length(p) >= 2) p[2] else NA_character_,
                     character(1))
  data.frame(dir_name = hits, ct = ct, contrast = contrast,
             stringsAsFactors = FALSE)
}

## Read one HOMER knownResults.txt
read.homer <- function(f, top_n = NULL) {
  if (!file.exists(f)) return(NULL)
  t <- read.csv(f, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  if (nrow(t) == 0) return(NULL)
  rownames(t) <- make.unique(t$Motif.Name)
  if (!is.null(top_n)) t <- t[seq_len(min(top_n, nrow(t))), ]
  t
}

## Build ggData for one tissue
build_ggData <- function(units, cfg) {
  # motif union (top-N per unit)
  motifs.set <- c()
  for (i in seq_len(nrow(units))) {
    f <- file.path(cfg$homer_root, units$dir_name[i], "knownResults.txt")
    t <- read.homer(f, top_n = TOP_N); if (is.null(t)) next
    motifs.set <- union(motifs.set, rownames(t))
  }
  if (length(motifs.set) == 0) return(NULL)

  gg <- c()
  for (i in seq_len(nrow(units))) {
    f <- file.path(cfg$homer_root, units$dir_name[i], "knownResults.txt")
    t <- read.homer(f, top_n = NULL); if (is.null(t)) next
    sub <- t[motifs.set, ]
    pct.t <- as.numeric(gsub("%", "", sub$X..of.Target.Sequences.with.Motif))
    gg <- rbind(gg, data.frame(
      CellType   = units$ct[i],
      Contrast   = units$contrast[i],
      Motif      = simplify.motif(sub$Motif.Name),
      Target_pct = pct.t,
      Log10_pval = -log10(sub$P.value),
      stringsAsFactors = FALSE
    ))
  }
  if (is.null(gg)) return(NULL)
  inf.idx <- which(gg$Log10_pval == Inf)
  if (length(inf.idx) > 0)
    gg$Log10_pval[inf.idx] <- max(gg$Log10_pval[-inf.idx], na.rm = TRUE)
  list(gg = gg, motifs = motifs.set)
}

## Plot one tissue panel
opening.pal <- c("#e3e1e1", "#f57f7f", "red")   # Ralph opening palette

make.panel <- function(cfg, panel_letter = NULL) {
  units <- discover.units(cfg); if (is.null(units)) return(NULL)
  if (!is.null(cfg$keep_contrast))
    units <- units[units$contrast %in% cfg$keep_contrast, , drop = FALSE]
  if (nrow(units) == 0) { message("  [SKIP] no units after contrast filter: ", cfg$title); return(NULL) }

  res <- build_ggData(units, cfg); if (is.null(res)) return(NULL)
  gg  <- res$gg

  # x-axis = cell type, or treatment group (before "_vs_") for the pooled T cells
  x_field <- cfg$x_axis
  gg$Xvar <- gg[[x_field]]
  if (x_field == "Contrast")
    gg$Xvar <- gsub("_", " ", sub("_vs_.*$", "", gg$Xvar))

  # x-axis order: explicit x_order wins, then ct_order (for cell types),
  # else alphabetical.
  if (!is.null(cfg$x_order)) {
    xlev <- intersect(cfg$x_order, unique(gg$Xvar))
  } else if (x_field == "CellType" && !is.null(cfg$ct_order)) {
    xlev <- intersect(cfg$ct_order, unique(gg$Xvar))
  } else {
    xlev <- sort(unique(gg$Xvar))
  }
  # Drop any cell type / group not in the chosen order, otherwise it would
  # become an NA factor level and render as a spurious "NA" column
  gg <- gg[gg$Xvar %in% xlev, , drop = FALSE]
  gg$Xvar <- factor(gg$Xvar, levels = xlev)

  # motif order: most significant at top. Cap the y-axis at MAX_MOTIFS
  # by keeping only the most significant motifs across this tissue.
  motif.rank <- gg %>% group_by(Motif) %>%
    summarise(m = max(Log10_pval, na.rm = TRUE), .groups = "drop") %>%
    arrange(m)                                  # ascending: top of plot = last
  if (nrow(motif.rank) > MAX_MOTIFS)
    motif.rank <- tail(motif.rank, MAX_MOTIFS)  # keep MAX_MOTIFS most significant
  gg <- gg[gg$Motif %in% motif.rank$Motif, ]
  gg$Motif <- factor(gg$Motif, levels = motif.rank$Motif)

  base.title <- paste0(cfg$title, " - opening DARs (top ", TOP_N,
                       " motifs per group)")
  title.str  <- if (is.null(panel_letter)) base.title
                else paste0(panel_letter, "   ", base.title)

  # X-axis caption flags whether columns are cell types or comparison
  # groups (T cells are a single pool compared across age/treatment).
  xlab.txt <- if (x_field == "Contrast") "Comparison group (vs Young control)"
              else "Cell type"

  ggplot(gg, aes(x = Xvar, y = Motif,
                 size = Target_pct, colour = Log10_pval)) +
    geom_point(alpha = 0.9) +
    scale_size(range = c(3, 18), name = "Target %",
               limits = c(0, max(gg$Target_pct, na.rm = TRUE))) +
    scale_colour_gradientn(colours = opening.pal,
                           values  = rescale(c(0, 10, 100), c(0, 1)),
                           name    = expression(-log[10] * "(" * italic(P) * ")")) +
    guides(colour = guide_colourbar(order = 1),
           size   = guide_legend(order = 2)) +
    labs(title = title.str, x = xlab.txt, y = NULL) +
    theme_classic(base_size = 20, base_family = "sans") +
    theme(
      plot.title          = element_text(face = "bold", size = 18, hjust = 0,
                                         margin = margin(b = 12),
                                         family = "sans", colour = "black"),
      plot.title.position = "plot",   # align title to the whole-plot left edge
      axis.text.x       = element_text(angle = 45, hjust = 1, vjust = 1,
                                       size = 22, colour = "black", family = "sans"),
      axis.text.y       = element_text(size = 21, colour = "black", family = "sans"),
      axis.title.x      = element_text(size = 20, face = "bold", colour = "black",
                                       family = "sans", margin = margin(t = 10)),
      panel.border      = element_rect(colour = "black", fill = NA, linewidth = 0.8),
      axis.line         = element_blank(),
      axis.ticks        = element_line(colour = "black", linewidth = 0.8),
      axis.ticks.length = unit(5, "pt"),
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      legend.title      = element_text(size = 18, family = "sans", colour = "black"),
      legend.text       = element_text(size = 16, family = "sans", colour = "black"),
      legend.key.size   = unit(0.8, "cm"),
      plot.margin       = margin(16, 20, 14, 18)
    )
}

## Build + save
solo.list  <- list()
combo.list <- list()
for (nm in names(TISSUES)) {
  cfg <- TISSUES[[nm]]
  message("=== ", cfg$title, " ===")
  p.solo  <- make.panel(cfg, panel_letter = NULL)
  p.combo <- make.panel(cfg, panel_letter = cfg$letter)
  if (is.null(p.solo)) next
  solo.list[[nm]]  <- p.solo
  combo.list[[nm]] <- p.combo
}
if (length(solo.list) == 0) stop("No tissue panels could be built. Check HOMER paths.")

# Per-tissue standalone PDFs (wider so bubbles spread out horizontally)
W_SOLO <- 14; H_SOLO <- 10
for (nm in names(solo.list)) {
  ggsave(file.path(OUT_DIR, sprintf("Fig11_opening_solo_%s.pdf", nm)),
         solo.list[[nm]], width = W_SOLO, height = H_SOLO, device = cairo_pdf)
  ggsave(file.path(OUT_DIR, sprintf("Fig11_opening_solo_%s.png", nm)),
         solo.list[[nm]], width = W_SOLO, height = H_SOLO, dpi = 400)
}

# Combined 2x2 with A/B/C/D tags (wider so bubbles + titles have room)
combined <- wrap_plots(combo.list, ncol = 2)
W_COMBO <- 30; H_COMBO <- 22
ggsave(file.path(OUT_DIR, "Fig11_opening_combined.pdf"),
       combined, width = W_COMBO, height = H_COMBO, device = cairo_pdf)
ggsave(file.path(OUT_DIR, "Fig11_opening_combined.png"),
       combined, width = W_COMBO, height = H_COMBO, dpi = 400)

message("\nSaved Fig 11 (", length(solo.list), " solo + 1 combined)  ->  ", OUT_DIR)
if (interactive()) print(combined)

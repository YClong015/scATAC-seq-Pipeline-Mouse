#!/usr/bin/env Rscript
# Fig 9 - DAR burden per tissue: diverging bars (opening up/red, closing down/blue).
# Counts from DESeq2_all.tsv (padj < 0.05).

library(dplyr)
library(ggplot2)
library(stringr)

HAVE_FREAD <- requireNamespace("data.table", quietly = TRUE)
if (!HAVE_FREAD) message("data.table not installed, falling back to read.delim")

OUT_DIR  <- "/QRISdata/Q8448/Mouse_disease_data/DAR/Fig9_DAR_burden"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PADJ_CUT <- 0.05

## per-tissue DAR_tables dir, contrast to keep, x-axis cell-type order, facet label
TISSUES <- list(
  Kidney = list(
    dar_tables  = "/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_pseudobulk_Kidney_v5_DESeq2/DAR_tables",
    contrast    = "Day42_vs_Sham",
    ct_order    = c("PCT","PST","Injured_PT","TAL","DCT_CNT","PC_URO","LEUK","FIB"),
    panel_label = "(a) Kidney, Day 42 vs Sham"
  ),
  Aorta = list(
    dar_tables  = "/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_pseudobulk_DESeq2/DAR_pseudobulk_Aorta_DESeq2/DAR_tables",
    contrast    = "Challenge_vs_Control",
    ct_order    = c("Macrophages","Pericytes","SMC"),
    panel_label = "(b) Aorta, Challenge vs Control"
  ),
  Lung = list(
    dar_tables  = "/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_pseudobulk_DESeq2/DAR_pseudobulk_Lung_DESeq2/DAR_tables",
    contrast    = "Case_vs_Control",
    ct_order    = c("AT2","B","Ciliated","EC-vasc","Eosinophils","Fib",
                    "Mac-alv","Mac-inter","Mo-Ly6c+","NK","Pen","SMCs","T"),
    panel_label = "(c) Lung, Case vs Control"
  ),
  Tcell = list(
    dar_tables  = "/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_pseudobulk_DESeq2/DAR_pseudobulk_Tcells_DESeq2/DAR_tables",
    contrast    = "Young_chronic_vs_Young_control",
    ct_order    = c("Tcell"),
    panel_label = "(d) T cells, Young chronic vs Young control"
  )
)

## Helper: count opening/closing DARs from one DESeq2_all.tsv
# Only reads padj + log2FoldChange via data.table::fread when available.
count.dars <- function(f, padj_cut = PADJ_CUT) {
  tab <- tryCatch({
    if (HAVE_FREAD) {
      data.table::fread(f, sep = "\t", header = TRUE,
                        select = c("padj", "log2FoldChange"),
                        showProgress = FALSE, data.table = FALSE)
    } else {
      read.delim(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    }
  }, error = function(e) NULL)
  if (is.null(tab) || nrow(tab) == 0) return(c(opening = NA, closing = NA))
  tab <- tab[!is.na(tab$padj), ]
  c(
    opening = sum(tab$padj < padj_cut & tab$log2FoldChange > 0),
    closing = sum(tab$padj < padj_cut & tab$log2FoldChange < 0)
  )
}

## Build long-format data frame
all.rows <- list()
for (tissue in names(TISSUES)) {
  cfg <- TISSUES[[tissue]]
  if (!dir.exists(cfg$dar_tables)) {
    message("[SKIP] no DAR_tables dir for ", tissue, ": ", cfg$dar_tables)
    next
  }
  tsvs <- list.files(cfg$dar_tables,
                     pattern = paste0("__", cfg$contrast, "__.*_DESeq2_all\\.tsv$"),
                     full.names = TRUE)
  if (length(tsvs) == 0) {
    message("[SKIP] no TSVs match contrast ", cfg$contrast, " for ", tissue)
    next
  }
  message("=== ", tissue, " (", length(tsvs), " TSVs) ===")
  for (i in seq_along(tsvs)) {
    f <- tsvs[i]
    bn <- sub("_DESeq2_all\\.tsv$", "", basename(f))
    parts <- strsplit(bn, "__", fixed = TRUE)[[1]]
    if (length(parts) < 2) next
    ct <- parts[1]
    cnt <- count.dars(f)
    message(sprintf("  [%d/%d] %-15s  opening=%6d  closing=%6d",
                    i, length(tsvs), ct, cnt["opening"], cnt["closing"]))
    all.rows[[length(all.rows) + 1]] <- data.frame(
      tissue      = tissue,
      cell_type   = ct,
      contrast    = cfg$contrast,
      opening     = cnt["opening"],
      closing     = cnt["closing"],
      panel_label = cfg$panel_label,
      stringsAsFactors = FALSE
    )
  }
}
wide <- do.call(rbind, all.rows)

# Apply CT order per tissue + drop any CT not in ct_order
wide.keep <- list()
for (tissue in names(TISSUES)) {
  cfg <- TISSUES[[tissue]]
  sub <- wide[wide$tissue == tissue, ]
  sub <- sub[sub$cell_type %in% cfg$ct_order, ]
  sub$cell_type <- factor(sub$cell_type, levels = cfg$ct_order)
  wide.keep[[tissue]] <- sub
}
wide <- do.call(rbind, wide.keep)

# Long format: opening positive, closing negative (for diverging bars)
long <- rbind(
  data.frame(wide[, c("tissue","cell_type","panel_label")],
             direction = "Opening",
             signed_n  =  wide$opening,
             n         =  wide$opening),
  data.frame(wide[, c("tissue","cell_type","panel_label")],
             direction = "Closing",
             signed_n  = -wide$closing,
             n         =  wide$closing)
)
long$direction   <- factor(long$direction, levels = c("Opening","Closing"))
long$panel_label <- factor(long$panel_label,
                           levels = sapply(TISSUES, `[[`, "panel_label"))

## Plot
PAL <- c("Opening" = "#D62728",   # red
         "Closing" = "#1F77B4")   # blue

fmt.n <- function(x) format(x, big.mark = ",", scientific = FALSE)

pl <- ggplot(long, aes(x = cell_type, y = signed_n, fill = direction)) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.3) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.45) +
  geom_text(aes(label = ifelse(n == 0, "", fmt.n(n)),
                vjust = ifelse(direction == "Opening", -0.45, 1.35)),
            size = 5.2, fontface = "bold", colour = "black", family = "sans") +
  scale_fill_manual(values = PAL, name = NULL,
                    labels = c(expression("Opening (log"[2]*"FC > 0)"),
                               expression("Closing (log"[2]*"FC < 0)"))) +
  scale_y_continuous(labels = function(x) fmt.n(abs(x)),
                     expand = expansion(mult = c(0.30, 0.30))) +
  facet_wrap(~ panel_label, scales = "free", nrow = 2) +
  coord_cartesian(clip = "off") +
  labs(x = NULL,
       y = expression("Number of DARs (adjusted " * italic(P) * " < 0.05)")) +
  theme_classic(base_size = 17, base_family = "sans") +
  theme(
    strip.text       = element_text(face = "bold", size = 18, hjust = 0,
                                    margin = margin(b = 12), family = "sans"),
    strip.background = element_blank(),
    panel.spacing.x  = unit(2.0, "lines"),
    panel.spacing.y  = unit(3.0, "lines"),
    axis.text.x      = element_text(angle = 45, hjust = 1, vjust = 1,
                                    size = 16, family = "sans", colour = "black"),
    axis.text.y      = element_text(family = "sans", colour = "black", size = 15),
    axis.title.y     = element_text(margin = margin(r = 12), family = "sans", size = 17),
    axis.ticks       = element_line(colour = "black", linewidth = 0.6),
    axis.ticks.length = unit(4, "pt"),
    axis.line        = element_line(colour = "black", linewidth = 0.6),
    legend.text      = element_text(family = "sans", size = 16),
    legend.title     = element_blank(),
    legend.position  = "bottom",
    legend.justification = "center",
    legend.box.margin = margin(t = 8, b = 0),
    legend.key.size  = unit(0.7, "cm"),
    legend.spacing.x = unit(0.6, "cm"),
    plot.margin      = margin(14, 18, 12, 14)
  )

## Save
W_CM <- 38
H_CM <- 28

ggsave(file.path(OUT_DIR, "Fig9_DAR_burden.pdf"),
       pl, width = W_CM, height = H_CM, units = "cm")
ggsave(file.path(OUT_DIR, "Fig9_DAR_burden.png"),
       pl, width = W_CM, height = H_CM, units = "cm", dpi = 300)

write.csv(wide, file.path(OUT_DIR, "Fig9_DAR_burden.csv"), row.names = FALSE)

message("\nDone. Output ->\n  ", OUT_DIR)

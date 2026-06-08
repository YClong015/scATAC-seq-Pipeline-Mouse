### Fig 4 - scATAC-seq quality-control metrics across the four disease tissues
### Reads per-cell metadata from cached CSVs if present, else from the RDS objects
library(Seurat)
library(Signac)
library(ggplot2)
library(dplyr)
library(patchwork)
library(scales)

out.dir  <- "/QRISdata/Q8448/Mouse_disease_data/QC_figures"
meta.dir <- file.path(out.dir, "metadata_cache")
dir.create(out.dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(meta.dir, showWarnings = FALSE, recursive = TRUE)

paths <- list(
  Kidney = "/QRISdata/Q8448/Mouse_disease_data/Kidney/kidney_merged_universal.rds",
  Lung   = "/QRISdata/Q8448/Mouse_disease_data/Lung/lung_universal_new_pruned.rds",
  Aorta  = "/QRISdata/Q8448/Mouse_disease_data/Aorta/Aorta_integrated_universal.rds",
  Tcells = "/QRISdata/Q8448/Mouse_disease_data/Tcells/tcells_universal.rds"
)

## load each tissue, extract per-cell metadata, cache to CSV
meta.list <- list()
for (tissue in names(paths)) {
  cache.file <- file.path(meta.dir, paste0(tissue, "_metadata.csv"))

  if (file.exists(cache.file)) {
    md <- read.csv(cache.file, row.names = 1, check.names = FALSE)
    md$Tissue <- tissue
    meta.list[[tissue]] <- md
    message(tissue, ": ", nrow(md), " cells (cache)")
    next
  }

  obj <- readRDS(paths[[tissue]])

  ## assay name differs by tissue
  assay.name <- if ("peaks_universal" %in% names(obj@assays)) {
    "peaks_universal"
  } else if ("ATAC" %in% names(obj@assays)) {
    "ATAC"
  } else if ("peaks" %in% names(obj@assays)) {
    "peaks"
  } else {
    names(obj@assays)[1]
  }
  DefaultAssay(obj) <- assay.name

  md <- obj@meta.data
  md$Tissue <- tissue

  ## standardise QC column names across tissues
  if ("peak_region_fragments" %in% colnames(md)) {
    md$nCount <- md$peak_region_fragments
  } else if ("nCount_ATAC" %in% colnames(md)) {
    md$nCount <- md$nCount_ATAC
  } else {
    cnt.col <- grep("^nCount", colnames(md), value = TRUE)[1]
    if (!is.na(cnt.col)) md$nCount <- md[[cnt.col]]
  }
  if (!"TSS.enrichment" %in% colnames(md)) {
    tss.col <- grep("TSS", colnames(md), ignore.case = TRUE, value = TRUE)[1]
    if (!is.na(tss.col)) md$TSS.enrichment <- md[[tss.col]]
  }
  if (!"nucleosome_signal" %in% colnames(md)) {
    ns.col <- grep("nucleosome", colnames(md), ignore.case = TRUE, value = TRUE)[1]
    if (!is.na(ns.col)) md$nucleosome_signal <- md[[ns.col]]
  }
  if (!"pct_reads_in_peaks" %in% colnames(md)) {
    frip.col <- grep("pct_reads|FRiP|frip", colnames(md), ignore.case = TRUE, value = TRUE)[1]
    if (!is.na(frip.col)) md$pct_reads_in_peaks <- md[[frip.col]]
  }

  write.csv(md, cache.file)
  message(tissue, ": ", nrow(md), " cells (", assay.name, ")")
  meta.list[[tissue]] <- md
  rm(obj); gc()
}

all.meta <- bind_rows(meta.list)
all.meta$Tissue <- factor(all.meta$Tissue, levels = c("Kidney", "Lung", "Aorta", "Tcells"))

tissue.colors <- c(Kidney = "#2166AC", Lung = "#B2182B", Aorta = "#4DAC26", Tcells = "#D6604D")

## condition column differs by tissue; collapse to one Group_label column
group.col.map <- list(
  Kidney = "condition",
  Lung   = "Group",
  Aorta  = "Group",
  Tcells = "deMultliplex2_final_mapped"
)
all.meta$Group_label <- NA_character_
for (tissue in names(group.col.map)) {
  col <- group.col.map[[tissue]]
  idx <- all.meta$Tissue == tissue
  if (col %in% colnames(all.meta)) all.meta$Group_label[idx] <- as.character(all.meta[[col]][idx])
}

## disease groups first, controls last (bottom of the ggplot2 stack)
group.order <- c(
  "Day42", "Case", "Challenge", "Aged",
  "Day14", "Young chronic", "Young acute", "Juvenile",
  "Sham", "Control", "Young control"
)
all.meta$Group_label <- factor(all.meta$Group_label, levels = group.order)

## blues = control/young, warm = disease/aged
group.colors <- c(
  "Sham"          = "#4393C3",
  "Control"       = "#4393C3",
  "Day14"         = "#FDAE61",
  "Day42"         = "#D73027",
  "Case"          = "#D73027",
  "Challenge"     = "#D73027",
  "Juvenile"      = "#74C476",
  "Young control" = "#2166AC",
  "Young acute"   = "#FEE090",
  "Young chronic" = "#F46D43",
  "Aged"          = "#A50026"
)

## violin + boxplot helper
violin.box <- function(data, yvar, ylabel, log.scale = FALSE, ylim = NULL, title = NULL) {
  p <- ggplot(data, aes(x = Tissue, y = .data[[yvar]], fill = Tissue)) +
    geom_violin(trim = TRUE, alpha = 0.85, linewidth = 0.25) +
    geom_boxplot(width = 0.10, fill = "white", outlier.shape = NA, linewidth = 0.35) +
    scale_fill_manual(values = tissue.colors) +
    labs(x = NULL, y = ylabel, title = title) +
    theme_classic(base_size = 9) +
    theme(
      plot.title   = element_text(face = "bold", hjust = 0.5, size = 9),
      axis.text.x  = element_text(angle = 35, hjust = 1, size = 8, color = "black"),
      axis.text.y  = element_text(size = 8, color = "black"),
      axis.title.y = element_text(size = 8),
      axis.line    = element_line(linewidth = 0.4),
      axis.ticks   = element_line(linewidth = 0.4),
      legend.position = "none",
      panel.grid   = element_blank()
    )
  if (log.scale) p <- p + scale_y_log10()
  if (!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)
  p
}

pub.theme <- theme_classic(base_size = 9) +
  theme(
    plot.title   = element_text(face = "bold", size = 9, hjust = 0),
    axis.text.x  = element_text(angle = 35, hjust = 1, size = 8, color = "black"),
    axis.text.y  = element_text(size = 8, color = "black"),
    axis.title.y = element_text(size = 8),
    axis.line    = element_line(linewidth = 0.4),
    axis.ticks   = element_line(linewidth = 0.4),
    legend.position = "none",
    panel.grid   = element_blank()
  )

## mice per group: Kidney/Lung from metadata, Aorta/Tcells from experimental design
sample.col.map <- list(
  Kidney = "orig.ident",
  Lung   = "SampleID",
  Aorta  = "sample_id",
  Tcells = "deMultliplex2_final_mapped"
)

# Tcells (Xiaoli Chen, pers. comm.): Juvenile 2 pooled, Young x3 each, Aged 4
tcell.sample.counts <- tibble(
  Tissue      = "Tcells",
  Group_label = c("Juvenile", "Young control", "Young acute", "Young chronic", "Aged"),
  n_samples   = c(2, 3, 3, 3, 4)
)
# Aorta (Zhang et al., ATVB 2023): n=3 mice per group, pooled into 1 sample
aorta.sample.counts <- tibble(
  Tissue      = "Aorta",
  Group_label = c("Control", "Challenge"),
  n_samples   = c(3, 3)
)

sample.counts <- bind_rows(
  bind_rows(lapply(c("Kidney", "Lung"), function(tissue) {
    col  <- sample.col.map[[tissue]]
    gcol <- group.col.map[[tissue]]
    md   <- meta.list[[tissue]]
    if (!col %in% colnames(md) || !gcol %in% colnames(md)) return(NULL)
    md %>%
      distinct(.data[[col]], .data[[gcol]]) %>%
      setNames(c("sample_id", "Group_label")) %>%
      mutate(Tissue = tissue)
  })) %>%
    filter(!is.na(Group_label)) %>%
    group_by(Tissue, Group_label) %>%
    summarise(n_samples = n(), .groups = "drop"),
  aorta.sample.counts,
  tcell.sample.counts
) %>%
  mutate(
    Tissue      = factor(Tissue, levels = c("Kidney", "Lung", "Aorta", "Tcells")),
    Group_label = factor(Group_label, levels = group.order)
  )

pA <- ggplot(sample.counts, aes(x = Tissue, y = n_samples, fill = Group_label)) +
  geom_col(width = 0.55, color = "white", linewidth = 0.3, position = "stack") +
  geom_text(aes(label = n_samples), position = position_stack(vjust = 0.5),
            size = 2.8, fontface = "bold", color = "white") +
  scale_fill_manual(values = group.colors, breaks = group.order, drop = TRUE) +
  scale_y_continuous(breaks = scales::pretty_breaks(), expand = expansion(mult = c(0, 0.08))) +
  labs(x = NULL, y = "Number of mice", fill = NULL, title = "A  Mice per experimental group") +
  pub.theme +
  theme(legend.position = "right", legend.text = element_text(size = 7),
        legend.key.size = unit(0.35, "cm"), legend.spacing.y = unit(0.15, "cm"))

## cell composition by group (stacked bar)
group.counts <- all.meta %>%
  filter(!is.na(Group_label)) %>%
  group_by(Tissue, Group_label) %>%
  summarise(n_cells = n(), .groups = "drop")

pA2 <- ggplot(group.counts, aes(x = Tissue, y = n_cells, fill = Group_label)) +
  geom_col(width = 0.55, color = "white", linewidth = 0.3, position = "stack") +
  geom_text(aes(label = scales::comma(n_cells)), position = position_stack(vjust = 0.5),
            size = 2.2, fontface = "bold", color = "white") +
  scale_fill_manual(values = group.colors, breaks = group.order, drop = TRUE) +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Number of cells", fill = NULL,
       title = "B  Cell composition by experimental group") +
  pub.theme +
  theme(legend.position = "right", legend.text = element_text(size = 7),
        legend.key.size = unit(0.35, "cm"), legend.spacing.y = unit(0.15, "cm"))

## three QC-metric violins
pB <- if ("nCount" %in% colnames(all.meta))
  violin.box(all.meta %>% filter(!is.na(nCount) & nCount > 0),
             "nCount", "Fragments in peaks (log10)", log.scale = TRUE, title = "A  Fragments per cell")
pC <- if ("TSS.enrichment" %in% colnames(all.meta))
  violin.box(all.meta %>% filter(!is.na(TSS.enrichment) & TSS.enrichment > 0),
             "TSS.enrichment", "TSS enrichment score", title = "B  TSS enrichment score")
pD <- if ("nucleosome_signal" %in% colnames(all.meta))
  violin.box(all.meta %>% filter(!is.na(nucleosome_signal) & nucleosome_signal > 0 & nucleosome_signal < 4),
             "nucleosome_signal", "Nucleosome signal", title = "C  Nucleosome signal")

## save panels
ggsave(file.path(out.dir, "QC_sample_counts.pdf"), pA, width = 4.5, height = 3.5, units = "in")
ggsave(file.path(out.dir, "QC_sample_counts.png"), pA, width = 4.5, height = 3.5, units = "in", dpi = 300)
ggsave(file.path(out.dir, "QC_group_composition.pdf"), pA2, width = 4.5, height = 3.5, units = "in")
ggsave(file.path(out.dir, "QC_group_composition.png"), pA2, width = 4.5, height = 3.5, units = "in", dpi = 300)

violin.panels <- Filter(Negate(is.null), list(pB, pC, pD))
if (length(violin.panels) > 0) {
  combined.violin <- wrap_plots(violin.panels, nrow = 1) +
    plot_annotation(
      title    = "scATAC-seq quality control metrics across tissues",
      subtitle = "Violin plots show per-cell distributions; box shows median and IQR",
      theme    = theme(plot.title    = element_text(face = "bold", hjust = 0.5, size = 13),
                       plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey40"))
    )
  ggsave(file.path(out.dir, "QC_metrics.pdf"), combined.violin,
         width = length(violin.panels) * 2.5, height = 3.5, units = "in")
  ggsave(file.path(out.dir, "QC_metrics.png"), combined.violin,
         width = length(violin.panels) * 2.5, height = 3.5, units = "in", dpi = 300)
}

## median QC summary table
summary.tbl <- all.meta %>%
  group_by(Tissue) %>%
  summarise(
    n_cells        = n(),
    median_nCount  = if ("nCount" %in% names(.)) median(nCount, na.rm = TRUE) else NA,
    median_TSS     = if ("TSS.enrichment" %in% names(.)) median(TSS.enrichment, na.rm = TRUE) else NA,
    median_nuc_sig = if ("nucleosome_signal" %in% names(.)) median(nucleosome_signal, na.rm = TRUE) else NA,
    .groups = "drop"
  )
print(summary.tbl)
write.csv(summary.tbl, file.path(out.dir, "QC_summary_table.csv"), row.names = FALSE)

## cell-type composition per sample (T cells excluded: cached labels are clusters,
## not the Masopust-2026 cell_type_final annotation used elsewhere)
composition.tissues <- c("Kidney", "Lung", "Aorta")
celltype.col.map <- list(Kidney = "cell_type", Lung = "cell_type", Aorta = "cell_type")

make.prop.panel <- function(tissue, show.y = TRUE) {
  md      <- meta.list[[tissue]]
  smp.col <- sample.col.map[[tissue]]
  ct.col  <- celltype.col.map[[tissue]]
  grp.col <- group.col.map[[tissue]]
  md$ct_label <- as.character(md[[ct.col]])
  if (!all(c(smp.col, grp.col) %in% colnames(md))) {
    message("  skip ", tissue, " (missing columns)"); return(NULL)
  }

  prop.df <- md %>%
    filter(!is.na(ct_label), !is.na(.data[[smp.col]])) %>%
    group_by(sample = .data[[smp.col]], condition = .data[[grp.col]], cell_type = ct_label) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(sample) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup() %>%
    mutate(condition = factor(condition, levels = group.order))

  ## order samples disease -> control
  smp.ord <- prop.df %>% distinct(sample, condition) %>% arrange(condition) %>% pull(sample) %>% unique()
  prop.df$sample <- factor(prop.df$sample, levels = smp.ord)

  all.ct  <- sort(unique(prop.df$cell_type))
  ct.cols <- setNames(colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(length(all.ct)), all.ct)

  ## colour x-axis labels by condition
  cond.key <- prop.df %>% distinct(sample, condition) %>%
    mutate(ax_col = group.colors[as.character(condition)])
  ax.cols <- cond.key$ax_col[match(levels(prop.df$sample), cond.key$sample)]

  ggplot(prop.df, aes(x = sample, y = prop, fill = cell_type)) +
    geom_col(width = 0.82, colour = NA) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       expand = expansion(mult = c(0, 0.02)), limits = c(0, 1)) +
    scale_fill_manual(values = ct.cols, name = NULL) +
    labs(x = NULL, y = if (show.y) "Cell type proportion" else NULL, title = tissue) +
    theme_classic(base_size = 9) +
    theme(
      plot.title   = element_text(face = "bold", size = 10, hjust = 0.5),
      axis.text.x  = element_text(angle = 45, hjust = 1, size = 7, colour = ax.cols),
      axis.text.y  = element_text(size = 7),
      axis.title.y = element_text(size = 8),
      axis.line    = element_line(linewidth = 0.4),
      axis.ticks   = element_line(linewidth = 0.4),
      legend.text  = element_text(size = 6),
      legend.key.size  = unit(0.26, "cm"),
      legend.spacing.y = unit(0.08, "cm"),
      panel.grid   = element_blank()
    )
}

prop.panels <- lapply(seq_along(composition.tissues),
                      function(i) make.prop.panel(composition.tissues[i], show.y = (i == 1)))
names(prop.panels) <- composition.tissues

pB.new <- wrap_plots(prop.panels, nrow = 1, widths = c(9, 6, 2)) +
  plot_annotation(title = "B  Cell type composition per sample",
                  theme = theme(plot.title = element_text(face = "bold", size = 10, hjust = 0)))
ggsave(file.path(out.dir, "QC_cell_composition.pdf"), pB.new, width = 16, height = 5)
ggsave(file.path(out.dir, "QC_cell_composition.png"), pB.new, width = 16, height = 5, dpi = 300)

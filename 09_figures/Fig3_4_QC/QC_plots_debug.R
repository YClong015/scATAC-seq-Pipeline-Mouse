suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(scales)
})

OUT_DIR  <- "/QRISdata/Q8448/Mouse_disease_data/QC_figures"
META_DIR <- file.path(OUT_DIR, "metadata_cache")  # cached CSVs to avoid reloading RDS
dir.create(OUT_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(META_DIR, showWarnings = FALSE, recursive = TRUE)

paths <- list(
  Kidney = "/QRISdata/Q8448/Mouse_disease_data/Kidney/kidney_merged_universal.rds",
  Lung   = "/QRISdata/Q8448/Mouse_disease_data/Lung/lung_universal_new_pruned.rds",
  Aorta  = "/QRISdata/Q8448/Mouse_disease_data/Aorta/Aorta_integrated_universal.rds",
  Tcells = "/QRISdata/Q8448/Mouse_disease_data/Tcells/tcells_universal.rds"
)

# ==============================================================================
# STEP 1: Load RDS files and extract metadata
#   - If a cached CSV exists for a tissue, load that instead (much faster)
#   - On first run: loads RDS, prints all column names, saves CSV cache
#   - On subsequent runs: reads CSV directly, skips RDS loading entirely
# ==============================================================================
meta_list <- list()

for (tissue in names(paths)) {
  cache_file <- file.path(META_DIR, paste0(tissue, "_metadata.csv"))

  # Load from cache if available (seconds vs. minutes)
  if (file.exists(cache_file)) {
    message("\n========== [CACHE] ", tissue, " ==========")
    md <- read.csv(cache_file, row.names = 1, check.names = FALSE)
    md$Tissue <- tissue
    meta_list[[tissue]] <- md
    message("  Loaded from cache: ", nrow(md), " cells")
    next
  }

  # First run: load full RDS
  message("\n========== [RDS] Loading: ", tissue, " ==========")
  obj <- readRDS(paths[[tissue]])

  # Auto-detect assay name
  assay_name <- if ("peaks_universal" %in% names(obj@assays)) {
    "peaks_universal"
  } else if ("ATAC" %in% names(obj@assays)) {
    "ATAC"
  } else if ("peaks" %in% names(obj@assays)) {
    "peaks"
  } else {
    names(obj@assays)[1]
  }
  DefaultAssay(obj) <- assay_name
  message("  Assay: ", assay_name, " | Cells: ", ncol(obj))

  md <- obj@meta.data
  md$Tissue <- tissue

  # Print all metadata columns for inspection
  message("  All metadata columns:")
  print(colnames(md))

  # Detect possible group / condition columns
  group_candidates <- grep(
    "condition|group|sample|disease|status|treatment|genotype|age|time|day|sham|case|ctrl|control",
    colnames(md), ignore.case = TRUE, value = TRUE
  )
  if (length(group_candidates) > 0) {
    message("  Possible group columns: ", paste(group_candidates, collapse = ", "))
    for (col in group_candidates) {
      message("    ", col, ": ", paste(unique(md[[col]]), collapse = " / "))
    }
  } else {
    message("  No obvious group column detected")
  }

  # Standardise QC columns across tissues
  # nCount: peak_region_fragments (Kidney/Lung/Aorta) or nCount_ATAC (T cells)
  if ("peak_region_fragments" %in% colnames(md)) {
    md$nCount <- md$peak_region_fragments
  } else if ("nCount_ATAC" %in% colnames(md)) {
    md$nCount <- md$nCount_ATAC
  } else {
    cnt_col <- grep("^nCount", colnames(md), value = TRUE)[1]
    if (!is.na(cnt_col)) md$nCount <- md[[cnt_col]]
  }

  if (!"TSS.enrichment" %in% colnames(md)) {
    tss_col <- grep("TSS", colnames(md), ignore.case = TRUE, value = TRUE)[1]
    if (!is.na(tss_col)) md$TSS.enrichment <- md[[tss_col]]
  }

  if (!"nucleosome_signal" %in% colnames(md)) {
    ns_col <- grep("nucleosome", colnames(md), ignore.case = TRUE, value = TRUE)[1]
    if (!is.na(ns_col)) md$nucleosome_signal <- md[[ns_col]]
  }

  if (!"pct_reads_in_peaks" %in% colnames(md)) {
    frip_col <- grep("pct_reads|FRiP|frip", colnames(md), ignore.case = TRUE, value = TRUE)[1]
    if (!is.na(frip_col)) md$pct_reads_in_peaks <- md[[frip_col]]
  }

  found_cols <- intersect(
    c("nCount", "TSS.enrichment", "nucleosome_signal", "pct_reads_in_peaks"),
    colnames(md)
  )
  message("  QC columns found: ", paste(found_cols, collapse = ", "))

  # Save cache so future runs skip RDS loading
  write.csv(md, cache_file)
  message("  Cache saved: ", cache_file)

  meta_list[[tissue]] <- md
  rm(obj); gc()
}

# ==============================================================================
# STEP 2: Combine metadata and set tissue factor order
#   After reviewing STEP 1 output, set GROUP_COL if you want to colour by
#   condition (e.g. GROUP_COL <- "condition"). Leave NULL to colour by tissue.
# ==============================================================================
all_meta <- bind_rows(meta_list)
all_meta$Tissue <- factor(all_meta$Tissue,
                          levels = c("Kidney", "Lung", "Aorta", "Tcells"))

tissue_colors <- c(
  Kidney = "#2166AC",
  Lung   = "#B2182B",
  Aorta  = "#4DAC26",
  Tcells = "#D6604D"
)

# Group column per tissue (confirmed from metadata inspection)
group_col_map <- list(
  Kidney = "condition",
  Lung   = "Group",
  Aorta  = "Group",
  Tcells = "deMultliplex2_final_mapped"
)

# Standardise group label into a new column "Group_label" (all Tcell groups kept)
all_meta$Group_label <- NA_character_
for (tissue in names(group_col_map)) {
  col <- group_col_map[[tissue]]
  idx <- all_meta$Tissue == tissue
  if (col %in% colnames(all_meta)) {
    all_meta$Group_label[idx] <- as.character(all_meta[[col]][idx])
  }
}

# Factor order: disease on top, control on bottom (bottom = first level in ggplot2 stack)
# Kidney: Day42 → Day14 → Sham (Sham at bottom)
# Lung:   Case → Control
# Aorta:  Challenge → Control
# Tcells: Aged → Young chronic → Young acute → Juvenile → Young control (control at bottom)
group_order <- c(
  "Day42", "Case", "Challenge", "Aged",
  "Day14", "Young chronic", "Young acute", "Juvenile",
  "Sham", "Control", "Young control"
)
all_meta$Group_label <- factor(all_meta$Group_label, levels = group_order)

# Publication-quality colorblind-friendly palette
# Blues = controls/young; warm = disease/aging progression
group_colors <- c(
  "Sham"          = "#4393C3",   # Kidney control
  "Control"       = "#4393C3",   # Lung / Aorta control
  "Day14"         = "#FDAE61",   # Kidney early disease
  "Day42"         = "#D73027",   # Kidney chronic disease
  "Case"          = "#D73027",   # Lung disease
  "Challenge"     = "#D73027",   # Aorta disease
  "Juvenile"      = "#74C476",   # Tcells: youngest
  "Young control" = "#2166AC",   # Tcells: young healthy
  "Young acute"   = "#FEE090",   # Tcells: young acute
  "Young chronic" = "#F46D43",   # Tcells: young chronic disease
  "Aged"          = "#A50026"    # Tcells: aged
)

# ==============================================================================
# STEP 3: Plot helper — violin + boxplot
# ==============================================================================
violin_box <- function(data, yvar, ylabel, log_scale = FALSE,
                       ylim = NULL, title = NULL) {
  p <- ggplot(data, aes(x = Tissue, y = .data[[yvar]], fill = Tissue)) +
    geom_violin(trim = TRUE, alpha = 0.85, linewidth = 0.25) +
    geom_boxplot(width = 0.10, fill = "white",
                 outlier.shape = NA, linewidth = 0.35) +
    scale_fill_manual(values = tissue_colors) +
    labs(x = NULL, y = ylabel, title = title) +
    theme_classic(base_size = 9) +
    theme(
      plot.title       = element_text(face = "bold", hjust = 0.5, size = 9),
      axis.text.x      = element_text(angle = 35, hjust = 1, size = 8, color = "black"),
      axis.text.y      = element_text(size = 8, color = "black"),
      axis.title.y     = element_text(size = 8),
      axis.line        = element_line(linewidth = 0.4),
      axis.ticks       = element_line(linewidth = 0.4),
      legend.position  = "none",
      panel.grid       = element_blank()
    )
  if (log_scale) p <- p + scale_y_log10()
  if (!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)
  p
}

# ==============================================================================
# STEP 4: Build panels
# ==============================================================================

# Publication theme base (used for all bar plots)
pub_theme <- theme_classic(base_size = 9) +
  theme(
    plot.title       = element_text(face = "bold", size = 9, hjust = 0),
    axis.text.x      = element_text(angle = 35, hjust = 1, size = 8, color = "black"),
    axis.text.y      = element_text(size = 8, color = "black"),
    axis.title.y     = element_text(size = 8),
    axis.line        = element_line(linewidth = 0.4),
    axis.ticks       = element_line(linewidth = 0.4),
    legend.position  = "none",
    panel.grid       = element_blank()
  )

# Panel A: number of biological samples per tissue per group
sample_col_map <- list(
  Kidney = "orig.ident",
  Lung   = "SampleID",
  Aorta  = "sample_id",
  Tcells = "deMultliplex2_final_mapped"
)

# Tcells: hard-coded from experimental design (Xiaoli Chen, pers. comm.)
# Juvenile: 2 pooled samples (Ab1=3 mice, Ab2=4 mice)
# Young control/acute/chronic: 3 individual mice each
# Aged: 4 individual mice
tcell_sample_counts <- tibble(
  Tissue      = "Tcells",
  Group_label = c("Juvenile", "Young control", "Young acute", "Young chronic", "Aged"),
  n_samples   = c(2, 3, 3, 3, 4)
)

# Aorta: hard-coded — 3 mice per group, pooled into 1 sequencing sample per group
# (Zhang et al., ATVB 2023: "samples harvested from each group (n=3) and pooled")
aorta_sample_counts <- tibble(
  Tissue      = "Aorta",
  Group_label = c("Control", "Challenge"),
  n_samples   = c(3, 3)
)

sample_counts <- bind_rows(
  # Kidney, Lung: derive from metadata (individual mice, 1 SRR per mouse)
  bind_rows(lapply(c("Kidney","Lung"), function(tissue) {
    col  <- sample_col_map[[tissue]]
    gcol <- group_col_map[[tissue]]
    md   <- meta_list[[tissue]]
    if (!col %in% colnames(md) || !gcol %in% colnames(md)) return(NULL)
    md %>%
      distinct(.data[[col]], .data[[gcol]]) %>%
      setNames(c("sample_id", "Group_label")) %>%
      mutate(Tissue = tissue)
  })) %>%
    filter(!is.na(Group_label)) %>%
    group_by(Tissue, Group_label) %>%
    summarise(n_samples = n(), .groups = "drop"),
  # Aorta: hard-coded (pooled, n=3 mice per group)
  aorta_sample_counts,
  # Tcells: hard-coded
  tcell_sample_counts
) %>%
  mutate(
    Tissue      = factor(Tissue, levels = c("Kidney","Lung","Aorta","Tcells")),
    Group_label = factor(Group_label, levels = group_order)
  )

pA <- ggplot(sample_counts,
             aes(x = Tissue, y = n_samples, fill = Group_label)) +
  geom_col(width = 0.55, color = "white", linewidth = 0.3,
           position = "stack") +
  geom_text(aes(label = n_samples),
            position = position_stack(vjust = 0.5),
            size = 2.8, fontface = "bold", color = "white") +
  scale_fill_manual(values = group_colors, breaks = group_order, drop = TRUE) +
  scale_y_continuous(breaks = scales::pretty_breaks(),
                     expand = expansion(mult = c(0, 0.08))) +
  labs(x = NULL, y = "Number of mice", fill = NULL,
       title = "A  Mice per experimental group") +
  pub_theme +
  theme(
    legend.position  = "right",
    legend.text      = element_text(size = 7),
    legend.key.size  = unit(0.35, "cm"),
    legend.spacing.y = unit(0.15, "cm")
  )

# Panel B: group composition stacked bar chart
group_counts <- all_meta %>%
  filter(!is.na(Group_label)) %>%
  group_by(Tissue, Group_label) %>%
  summarise(n_cells = n(), .groups = "drop")

pA2 <- ggplot(group_counts,
              aes(x = Tissue, y = n_cells, fill = Group_label)) +
  geom_col(width = 0.55, color = "white", linewidth = 0.3,
           position = "stack") +
  geom_text(aes(label = scales::comma(n_cells)),
            position = position_stack(vjust = 0.5),
            size = 2.2, fontface = "bold", color = "white") +
  scale_fill_manual(values = group_colors, breaks = group_order,
                    drop = TRUE) +
  scale_y_continuous(labels = scales::comma,
                     expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Number of cells", fill = NULL,
       title = "B  Cell composition by experimental group") +
  pub_theme +
  theme(
    legend.position  = "right",
    legend.text      = element_text(size = 7),
    legend.key.size  = unit(0.35, "cm"),
    legend.spacing.y = unit(0.15, "cm")
  )

# Panel B: fragments per cell (log10)
pB <- NULL
if ("nCount" %in% colnames(all_meta)) {
  pB <- violin_box(
    all_meta %>% filter(!is.na(nCount) & nCount > 0),
    yvar      = "nCount",
    ylabel    = "Fragments in peaks (log10)",
    log_scale = TRUE,
    title     = "A  Fragments per cell"
  )
}

# Panel C: TSS enrichment score
pC <- NULL
if ("TSS.enrichment" %in% colnames(all_meta)) {
  pC <- violin_box(
    all_meta %>% filter(!is.na(TSS.enrichment) & TSS.enrichment > 0),
    yvar   = "TSS.enrichment",
    ylabel = "TSS enrichment score",
    title  = "B  TSS enrichment score"
  )
}

# Panel D: nucleosome signal (filter to < 4 for consistent display across tissues)
pD <- NULL
if ("nucleosome_signal" %in% colnames(all_meta)) {
  pD <- violin_box(
    all_meta %>% filter(!is.na(nucleosome_signal) & nucleosome_signal > 0 & nucleosome_signal < 4),
    yvar   = "nucleosome_signal",
    ylabel = "Nucleosome signal",
    title  = "C  Nucleosome signal"
  )
}

# ==============================================================================
# STEP 5: Save Panel A separately, save Panels B+C+D together
# ==============================================================================

# --- Panel A: biological sample counts ---
ggsave(file.path(OUT_DIR, "QC_sample_counts.pdf"),
       pA, width = 4.5, height = 3.5, units = "in")
ggsave(file.path(OUT_DIR, "QC_sample_counts.png"),
       pA, width = 4.5, height = 3.5, units = "in", dpi = 300)
message("Saved: QC_sample_counts (.pdf + .png)")

# --- Panel B standalone: group composition (wider for legend) ---
ggsave(file.path(OUT_DIR, "QC_group_composition.pdf"),
       pA2, width = 4.5, height = 3.5, units = "in")
ggsave(file.path(OUT_DIR, "QC_group_composition.png"),
       pA2, width = 4.5, height = 3.5, units = "in", dpi = 300)
message("Saved: QC_group_composition (.pdf + .png)")

# --- Panels B + C + D combined ---
violin_panels <- Filter(Negate(is.null), list(pB, pC, pD))
n_violin <- length(violin_panels)

if (n_violin > 0) {
  combined_violin <- wrap_plots(violin_panels, nrow = 1) +
    plot_annotation(
      title    = "scATAC-seq quality control metrics across tissues",
      subtitle = "Violin plots show per-cell distributions; box shows median and IQR",
      theme    = theme(
        plot.title    = element_text(face = "bold", hjust = 0.5, size = 13),
        plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey40")
      )
    )
  ggsave(file.path(OUT_DIR, "QC_metrics.pdf"),
         combined_violin, width = n_violin * 2.5, height = 3.5, units = "in")
  ggsave(file.path(OUT_DIR, "QC_metrics.png"),
         combined_violin, width = n_violin * 2.5, height = 3.5, units = "in", dpi = 300)
  message("Saved: QC_metrics (.pdf + .png) -> ", OUT_DIR)
}

# Print and save summary table
summary_tbl <- all_meta %>%
  group_by(Tissue) %>%
  summarise(
    n_cells        = n(),
    median_nCount  = if ("nCount" %in% names(.)) median(nCount, na.rm = TRUE) else NA,
    median_TSS     = if ("TSS.enrichment" %in% names(.)) median(TSS.enrichment, na.rm = TRUE) else NA,
    median_nuc_sig = if ("nucleosome_signal" %in% names(.)) median(nucleosome_signal, na.rm = TRUE) else NA,
    .groups = "drop"
  )
print(summary_tbl)
write.csv(summary_tbl, file.path(OUT_DIR, "QC_summary_table.csv"), row.names = FALSE)
message("Saved: QC_summary_table.csv")

# ==============================================================================
# NEW Panel B: Cell type proportion per sample (atlas style)
# X axis  = each biological sample, coloured by condition
# Y axis  = proportion 0–100 %
# Fill    = cell type (each tissue gets its own palette)
# ==============================================================================

# Tcells: map seurat_clusters → cell type label (no cell_type column in cache)
tcell_ct_map <- c(
  "0"  = "Effector T cell",    "1"  = "Naive T cell",
  "2"  = "B cell",             "3"  = "Treg",
  "4"  = "CD8+ T cell",        "5"  = "Activated T cell",
  "6"  = "CD8+ effector",      "7"  = "Cycling T cell",
  "8"  = "Innate-like T cell", "9"  = "NK cell",
  "10" = "Memory T cell",      "11" = "Memory T cell",
  "12" = "Memory T cell",      "13" = "Memory T cell",
  "14" = "CD8+ effector"
)

celltype_col_map <- list(
  Kidney = "cell_type",
  Lung   = "cell_type",
  Aorta  = "cell_type",
  Tcells = "seurat_clusters"   # mapped via tcell_ct_map below
)

# ── build proportion panel for one tissue ──────────────────────────────────
make_prop_panel <- function(tissue, show_y = TRUE) {
  md      <- meta_list[[tissue]]
  smp_col <- sample_col_map[[tissue]]
  ct_col  <- celltype_col_map[[tissue]]
  grp_col <- group_col_map[[tissue]]

  # Resolve cell type label
  if (tissue == "Tcells") {
    md$ct_label <- unname(tcell_ct_map[as.character(md[[ct_col]])])
  } else {
    md$ct_label <- as.character(md[[ct_col]])
  }

  if (!all(c(smp_col, grp_col) %in% colnames(md))) {
    message("  [SKIP] missing columns for ", tissue); return(NULL)
  }

  # Compute per-sample proportions
  prop_df <- md %>%
    filter(!is.na(ct_label), !is.na(.data[[smp_col]])) %>%
    group_by(sample    = .data[[smp_col]],
             condition = .data[[grp_col]],
             cell_type = ct_label) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(sample) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup() %>%
    mutate(condition = factor(condition, levels = group_order))

  # Order samples: disease first → control last (matches stacked bar logic)
  smp_ord <- prop_df %>%
    distinct(sample, condition) %>%
    arrange(condition) %>%
    pull(sample) %>%
    unique()
  prop_df$sample <- factor(prop_df$sample, levels = smp_ord)

  # Per-tissue cell type colour palette (Paired → enough for ~15 types)
  all_ct  <- sort(unique(prop_df$cell_type))
  ct_cols <- setNames(
    colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(length(all_ct)),
    all_ct
  )

  # X-axis label colour = condition colour
  cond_key  <- prop_df %>% distinct(sample, condition) %>%
    mutate(ax_col = group_colors[as.character(condition)])
  ax_cols <- cond_key$ax_col[match(levels(prop_df$sample), cond_key$sample)]

  ggplot(prop_df, aes(x = sample, y = prop, fill = cell_type)) +
    geom_col(width = 0.82, colour = NA) +
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      expand = expansion(mult = c(0, 0.02)),
      limits = c(0, 1)
    ) +
    scale_fill_manual(values = ct_cols, name = NULL) +
    labs(x = NULL,
         y = if (show_y) "Cell type proportion" else NULL,
         title = tissue) +
    theme_classic(base_size = 9) +
    theme(
      plot.title       = element_text(face = "bold", size = 10, hjust = 0.5),
      axis.text.x      = element_text(angle = 45, hjust = 1, size = 7,
                                      colour = ax_cols),
      axis.text.y      = element_text(size = 7),
      axis.title.y     = element_text(size = 8),
      axis.line        = element_line(linewidth = 0.4),
      axis.ticks       = element_line(linewidth = 0.4),
      legend.text      = element_text(size = 6),
      legend.key.size  = unit(0.26, "cm"),
      legend.spacing.y = unit(0.08, "cm"),
      panel.grid       = element_blank()
    )
}

# ── build & preview each panel individually ────────────────────────────────
prop_panels <- lapply(seq_along(names(meta_list)), function(i) {
  make_prop_panel(names(meta_list)[i], show_y = (i == 1))
})
names(prop_panels) <- names(meta_list)

# Preview in RStudio Plots panel — run each line separately
print(prop_panels$Kidney)
print(prop_panels$Lung)
print(prop_panels$Aorta)
print(prop_panels$Tcells)

# ── combined figure (widths ∝ number of samples) ───────────────────────────
pB_new <- wrap_plots(prop_panels, nrow = 1, widths = c(9, 6, 2, 5)) +
  plot_annotation(
    title = "B  Cell type composition per sample",
    theme = theme(plot.title = element_text(face = "bold", size = 10, hjust = 0))
  )
print(pB_new)

ggsave(file.path(OUT_DIR, "QC_cell_composition.pdf"),
       pB_new, width = 16, height = 5)
ggsave(file.path(OUT_DIR, "QC_cell_composition.png"),
       pB_new, width = 16, height = 5, dpi = 300)
message("Saved: QC_cell_composition (.pdf + .png)")

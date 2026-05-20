#!/usr/bin/env Rscript
#
# Integration cluster re-annotation (simplified — no GeneActivity needed)
#
# HOW TO USE:
#   Step 1 — set PHASE <- 1, run script → review plots in DIAG_DIR
#   Step 2 — fill in cluster_map below using ONLY the 21 allowed cell types
#   Step 3 — set PHASE <- 2, run script → annotated UMAP saved to OUT_DIR
#
# Allowed cell types (21 total):
#   Lung_AT2, Lung_B, Lung_Ciliated, Lung_EC-vasc, Lung_Eosinophils,
#   Lung_Fib, Lung_Mac-alv, Lung_Mac-inter, Lung_NK, Lung_Pen, Lung_SMCs,
#   Lung_T, Aorta_Macrophages, Aorta_Pericytes, Aorta_SMC,
#   Kidney_DCT, Kidney_Macrophages, Kidney_PC, Kidney_PT, Kidney_TAL,
#   Tcell

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(ggrepel)
})

# ==============================================================================
# CONFIGURATION
# ==============================================================================
PHASE   <- 1   # <-- change to 2 after filling cluster_map

RDS     <- "/QRISdata/Q8448/Mouse_disease_data/Integrated/All_Tissues_Integrated_.rds"
OUT_DIR <- "/QRISdata/Q8448/Mouse_disease_data/QC_figures"
DIAG_DIR <- file.path(OUT_DIR, "annotation_diagnostics")
dir.create(DIAG_DIR, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# CLUSTER MAP — fill in AFTER reviewing Phase 1 plots
# Keys   = cluster number (string)
# Values = one of the 21 allowed types listed above
# ==============================================================================
cluster_map <- c(
  "0"  = "???",
  "1"  = "???",
  "2"  = "???",
  "3"  = "???",
  "4"  = "???",
  "5"  = "???",
  "6"  = "???",
  "7"  = "???",
  "8"  = "???",
  "9"  = "???",
  "10" = "???",
  "11" = "???",
  "12" = "???",
  "13" = "???",
  "14" = "???",
  "15" = "???",
  "16" = "???",
  "17" = "???",
  "18" = "???",
  "19" = "???",
  "20" = "???",
  "21" = "???",
  "22" = "???",
  "23" = "???",
  "24" = "???",
  "25" = "???",
  "26" = "???",
  "27" = "???",
  "28" = "???"
)

# ==============================================================================
# Load object
# ==============================================================================
message("Loading integrated object...")
obj <- readRDS(RDS)
message(sprintf("Clusters: %d  |  Cells: %d",
                length(unique(obj$seurat_clusters)), ncol(obj)))

# ==============================================================================
# PHASE 1: DIAGNOSTIC PLOTS (no GeneActivity — uses existing metadata only)
# ==============================================================================
if (PHASE == 1) {
  message("\n=== PHASE 1: Generating diagnostic plots ===")

  clusters <- as.character(obj$seurat_clusters)
  tissue   <- obj$Tissue
  ct       <- if ("cell_type" %in% colnames(obj@meta.data))
                as.character(obj$cell_type)
              else
                clusters

  # --------------------------------------------------------------------------
  # 1. Cluster × Tissue composition (fraction bar chart + counts CSV)
  # --------------------------------------------------------------------------
  comp_df <- as.data.frame(table(Cluster = clusters, Tissue = tissue))
  comp_df$Cluster <- factor(comp_df$Cluster,
                            levels = as.character(sort(as.numeric(unique(comp_df$Cluster)))))

  p_tissue <- ggplot(comp_df, aes(x = Cluster, y = Freq, fill = Tissue)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = c(Kidney="#E41A1C", Lung="#377EB8",
                                 Aorta="#4DAF4A", Tcells="#FF7F00")) +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(title = "Cluster composition by tissue",
         x = "Cluster", y = "Fraction") +
    theme_classic(base_size = 11) +
    theme(axis.text.x  = element_text(size = 8),
          plot.title   = element_text(face = "bold", hjust = 0.5))

  comp_wide <- comp_df %>%
    tidyr::pivot_wider(names_from = Tissue, values_from = Freq, values_fill = 0)
  write.csv(comp_wide,
            file.path(DIAG_DIR, "cluster_tissue_counts.csv"), row.names = FALSE)

  # --------------------------------------------------------------------------
  # 2. Cluster × existing cell_type (top cell type per cluster)
  # --------------------------------------------------------------------------
  ct_df <- as.data.frame(table(Cluster = clusters, CellType = ct))
  ct_df$Cluster <- factor(ct_df$Cluster,
                          levels = as.character(sort(as.numeric(unique(ct_df$Cluster)))))

  # Top cell type per cluster (for easy reading)
  top_ct <- ct_df %>%
    group_by(Cluster) %>%
    slice_max(order_by = Freq, n = 1, with_ties = FALSE) %>%
    rename(TopCellType = CellType, TopCount = Freq)

  ct_wide <- ct_df %>%
    tidyr::pivot_wider(names_from = CellType, values_from = Freq, values_fill = 0)
  write.csv(ct_wide,
            file.path(DIAG_DIR, "cluster_celltype_counts.csv"), row.names = FALSE)

  # Stacked bar: cell type fraction per cluster
  n_ct <- length(unique(ct))
  ct_pal <- c(
    "#E41A1C","#377EB8","#4DAF4A","#FF7F00","#984EA3",
    "#A65628","#F781BF","#8DD3C7","#BEBADA","#FB8072",
    "#80B1D3","#FDB462","#B3DE69","#FCCDE5","#BC80BD",
    "#CCEBC5","#FFED6F","#1F78B4","#33A02C","#FB9A99",
    "#6A3D9A","#B15928","#CAB2D6","#FFFF99","#B2DF8A",
    "#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854"
  )[seq_len(n_ct)]

  p_celltype <- ggplot(ct_df, aes(x = Cluster, y = Freq, fill = CellType)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = setNames(ct_pal, sort(unique(ct)))) +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(title = "Cluster composition by existing cell_type",
         x = "Cluster", y = "Fraction") +
    guides(fill = guide_legend(ncol = 2, keyheight = unit(0.35, "cm"))) +
    theme_classic(base_size = 10) +
    theme(axis.text.x = element_text(size = 7),
          plot.title  = element_text(face = "bold", hjust = 0.5),
          legend.text = element_text(size = 7))

  # --------------------------------------------------------------------------
  # 3. Summary table: cluster → tissue dominant + top cell type
  # --------------------------------------------------------------------------
  tiss_df <- as.data.frame(table(Cluster = clusters, Tissue = tissue))
  tiss_top <- tiss_df %>%
    group_by(Cluster) %>%
    slice_max(order_by = Freq, n = 1, with_ties = FALSE) %>%
    rename(DomTissue = Tissue, TissueCells = Freq)

  summary_tbl <- tiss_top %>%
    left_join(top_ct, by = "Cluster") %>%
    left_join(
      as.data.frame(table(Cluster = clusters)) %>% rename(Total = Freq),
      by = "Cluster"
    ) %>%
    mutate(TissuePct = round(100 * TissueCells / Total, 1)) %>%
    select(Cluster, Total, DomTissue, TissuePct, TopCellType, TopCount) %>%
    arrange(as.numeric(as.character(Cluster)))

  write.csv(summary_tbl,
            file.path(DIAG_DIR, "cluster_annotation_guide.csv"), row.names = FALSE)

  message("\n=== Cluster annotation guide ===")
  print(as.data.frame(summary_tbl), row.names = FALSE)

  # --------------------------------------------------------------------------
  # 4. UMAP panels: by cluster (labeled) / by tissue / by existing cell_type
  # --------------------------------------------------------------------------
  emb <- as.data.frame(Embeddings(obj, "umap.harmony"))
  colnames(emb) <- c("UMAP_1", "UMAP_2")
  emb$Cluster  <- clusters
  emb$Tissue   <- tissue
  emb$CellType <- ct
  set.seed(42); emb <- emb[sample(nrow(emb)), ]

  base_theme <- theme_classic(base_size = 10) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          axis.line = element_blank(), axis.title = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5, size = 11))

  cluster_centroids <- emb %>%
    group_by(Cluster) %>%
    summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = "drop")

  p_umap_cl <- ggplot(emb, aes(UMAP_1, UMAP_2, color = Cluster)) +
    geom_point(size = 0.2, alpha = 0.4, stroke = 0) +
    ggrepel::geom_label_repel(
      data = cluster_centroids, aes(label = Cluster),
      size = 3, fontface = "bold", color = "black",
      fill = alpha("white", 0.8), label.size = 0.2,
      max.overlaps = 60, show.legend = FALSE
    ) +
    ggtitle("Harmony clusters") +
    theme(legend.position = "none") +
    base_theme

  p_umap_tis <- ggplot(emb, aes(UMAP_1, UMAP_2, color = Tissue)) +
    geom_point(size = 0.2, alpha = 0.5, stroke = 0) +
    scale_color_manual(values = c(Kidney="#E41A1C", Lung="#377EB8",
                                  Aorta="#4DAF4A", Tcells="#FF7F00")) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    ggtitle("Tissue") + base_theme

  ct_centroids <- emb %>%
    group_by(CellType) %>%
    summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = "drop")

  p_umap_ct <- ggplot(emb, aes(UMAP_1, UMAP_2, color = CellType)) +
    geom_point(size = 0.2, alpha = 0.4, stroke = 0) +
    scale_color_manual(values = setNames(ct_pal, sort(unique(ct)))) +
    ggrepel::geom_label_repel(
      data = ct_centroids, aes(label = CellType),
      size = 2.3, fontface = "bold", color = "black",
      fill = alpha("white", 0.8), label.size = 0.15,
      max.overlaps = 50, show.legend = FALSE
    ) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1),
                                ncol = 2, keyheight = unit(0.35, "cm"))) +
    ggtitle("Existing cell_type") + base_theme +
    theme(legend.text = element_text(size = 7))

  combined_diag <- p_umap_cl + p_umap_tis + p_umap_ct +
    plot_layout(nrow = 1, widths = c(1, 1, 1.4))

  # --------------------------------------------------------------------------
  # Save
  # --------------------------------------------------------------------------
  ggsave(file.path(DIAG_DIR, "01_umap_clusters_tissue_celltype.pdf"),
         combined_diag, width = 18, height = 6, useDingbats = FALSE)
  ggsave(file.path(DIAG_DIR, "02_cluster_tissue_composition.pdf"),
         p_tissue, width = 12, height = 4, useDingbats = FALSE)
  ggsave(file.path(DIAG_DIR, "03_cluster_celltype_composition.pdf"),
         p_celltype, width = 14, height = 5, useDingbats = FALSE)

  message("\n=== Phase 1 complete ===")
  message("Files in: ", DIAG_DIR)
  message("  01_umap_clusters_tissue_celltype.pdf")
  message("  02_cluster_tissue_composition.pdf")
  message("  03_cluster_celltype_composition.pdf")
  message("  cluster_annotation_guide.csv  <-- start here")
  message("\nNext: fill cluster_map at top of script, set PHASE <- 2, re-run.")
  stop("Phase 1 done — set PHASE <- 2 to continue.", call. = FALSE)
}

# ==============================================================================
# PHASE 2: APPLY ANNOTATION AND GENERATE UMAP FIGURE
# ==============================================================================
if (PHASE == 2) {

  allowed <- c(
    "Lung_AT2", "Lung_B", "Lung_Ciliated", "Lung_EC-vasc",
    "Lung_Eosinophils", "Lung_Fib", "Lung_Mac-alv", "Lung_Mac-inter",
    "Lung_NK", "Lung_Pen", "Lung_SMCs", "Lung_T",
    "Aorta_Macrophages", "Aorta_Pericytes", "Aorta_SMC",
    "Kidney_DCT", "Kidney_Macrophages", "Kidney_PC", "Kidney_PT", "Kidney_TAL",
    "Tcell"
  )

  unfilled <- cluster_map[cluster_map == "???"]
  if (length(unfilled) > 0)
    stop("Unfilled clusters in cluster_map: ",
         paste(names(unfilled), collapse = ", "))

  invalid <- cluster_map[!cluster_map %in% allowed]
  if (length(invalid) > 0)
    stop("Invalid labels (not in allowed list): ",
         paste(unique(invalid), collapse = ", "))

  missing_cl <- setdiff(as.character(unique(obj$seurat_clusters)), names(cluster_map))
  if (length(missing_cl) > 0)
    stop("cluster_map missing clusters: ", paste(missing_cl, collapse = ", "))

  obj$cell_type_annotated <- cluster_map[as.character(obj$seurat_clusters)]
  message("Annotation applied:")
  print(table(obj$cell_type_annotated))

  ann_path <- file.path(dirname(RDS), "All_Tissues_Integrated_annotated.rds")
  message("Saving annotated object: ", ann_path)
  saveRDS(obj, ann_path)

  # --------------------------------------------------------------------------
  # Publication UMAP figure
  # --------------------------------------------------------------------------
  tissue_colors <- c(Kidney="#E41A1C", Lung="#377EB8",
                     Aorta="#4DAF4A", Tcells="#FF7F00")

  celltype_palette <- c(
    "#E41A1C","#377EB8","#4DAF4A","#FF7F00","#984EA3",
    "#A65628","#F781BF","#8DD3C7","#BEBADA","#FB8072",
    "#80B1D3","#FDB462","#B3DE69","#FCCDE5","#BC80BD",
    "#CCEBC5","#FFED6F","#1F78B4","#33A02C","#FB9A99",
    "#6A3D9A","#B15928","#CAB2D6","#FFFF99","#B2DF8A"
  )

  pub_theme <- theme_classic(base_size = 11) +
    theme(
      plot.title       = element_text(face = "bold", size = 12, hjust = 0.5,
                                      margin = margin(b = 4)),
      axis.text        = element_blank(),
      axis.ticks       = element_blank(),
      axis.line        = element_blank(),
      axis.title       = element_blank(),
      legend.text      = element_text(size = 9),
      legend.title     = element_blank(),
      legend.key.size  = unit(0.45, "cm"),
      legend.spacing.y = unit(0.15, "cm"),
      panel.grid       = element_blank(),
      plot.margin      = margin(6, 8, 6, 6)
    )

  lsi_emb <- as.data.frame(Embeddings(obj, "umap.lsi"))
  colnames(lsi_emb) <- c("UMAP_1", "UMAP_2")
  lsi_emb$Tissue <- obj$Tissue
  set.seed(42); lsi_emb <- lsi_emb[sample(nrow(lsi_emb)), ]

  har_emb <- as.data.frame(Embeddings(obj, "umap.harmony"))
  colnames(har_emb) <- c("UMAP_1", "UMAP_2")
  har_emb$Tissue    <- obj$Tissue
  har_emb$cell_type <- obj$cell_type_annotated
  set.seed(42); har_emb <- har_emb[sample(nrow(har_emb)), ]

  p_before <- ggplot(lsi_emb, aes(UMAP_1, UMAP_2, color = Tissue)) +
    geom_point(size = 0.3, alpha = 0.6, stroke = 0) +
    scale_color_manual(values = tissue_colors) +
    guides(color = guide_legend(override.aes = list(size = 4, alpha = 1),
                                ncol = 1)) +
    ggtitle("Before integration") + pub_theme

  p_after_tissue <- ggplot(har_emb, aes(UMAP_1, UMAP_2, color = Tissue)) +
    geom_point(size = 0.3, alpha = 0.6, stroke = 0) +
    scale_color_manual(values = tissue_colors) +
    guides(color = guide_legend(override.aes = list(size = 4, alpha = 1),
                                ncol = 1)) +
    ggtitle("After integration (tissue)") + pub_theme

  ct_levels <- sort(unique(har_emb$cell_type))
  ct_colors <- setNames(celltype_palette[seq_along(ct_levels)], ct_levels)

  centroids <- har_emb %>%
    group_by(cell_type) %>%
    summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = "drop")

  p_after_ct <- ggplot(har_emb, aes(UMAP_1, UMAP_2, color = cell_type)) +
    geom_point(size = 0.3, alpha = 0.6, stroke = 0) +
    scale_color_manual(values = ct_colors) +
    ggrepel::geom_label_repel(
      data = centroids, aes(label = cell_type),
      size = 2.8, fontface = "bold", color = "black",
      fill = alpha("white", 0.82), label.size = 0.2,
      label.padding = unit(0.18, "lines"), box.padding = 0.4,
      max.overlaps = 50, min.segment.length = 0.2,
      segment.size = 0.35, show.legend = FALSE
    ) +
    ggtitle("After integration (cell type)") +
    guides(color = guide_legend(override.aes = list(size = 4, alpha = 1),
                                ncol = 2, byrow = TRUE)) +
    pub_theme + theme(legend.position = "right")

  combined <- p_before + p_after_tissue + p_after_ct +
    plot_layout(nrow = 1, widths = c(1, 1, 1.5)) +
    plot_annotation(
      tag_levels = "A",
      theme = theme(plot.tag = element_text(size = 13, face = "bold"))
    )

  ggsave(file.path(OUT_DIR, "Integration_A_before.pdf"),
         p_before, width = 4.5, height = 4, useDingbats = FALSE)
  ggsave(file.path(OUT_DIR, "Integration_B_after_tissue.pdf"),
         p_after_tissue, width = 4.5, height = 4, useDingbats = FALSE)
  ggsave(file.path(OUT_DIR, "Integration_C_after_celltype.pdf"),
         p_after_ct, width = 6.5, height = 4, useDingbats = FALSE)
  ggsave(file.path(OUT_DIR, "Fig_Integration_UMAP.pdf"),
         combined, width = 14, height = 4.8, useDingbats = FALSE)
  ggsave(file.path(OUT_DIR, "Fig_Integration_UMAP.png"),
         combined, width = 14, height = 4.8, dpi = 300)

  message("Done. Saved: Fig_Integration_UMAP (.pdf + .png)")
}

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
})

OUT_DIR   <- "/QRISdata/Q8448/Mouse_disease_data/QC_figures"
CACHE_DIR <- file.path(OUT_DIR, "umap_cache")
dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)

# Updated RDS path: new object processed by Tcell_scATAC.R
# Cell type annotations are already stored in the 'cell_type' metadata column
# UMAP is based on Harmony-corrected LSI (not WNN, not raw LSI)
TCELLS_RDS <- "/QRISdata/Q8448/Mouse_disease_data/Tcells/tcells_processed.rds"

# ==============================================================================
# Load T cell object
# ==============================================================================
message("Loading T cell object...")
obj <- readRDS(TCELLS_RDS)

# Confirm cell_type column exists
if (!"cell_type" %in% colnames(obj@meta.data)) {
  stop("'cell_type' column not found. Run Tcell_scATAC.R first to generate annotations.")
}

message(sprintf("Total cells: %d | Cell types: %s",
                ncol(obj),
                paste(sort(unique(obj$cell_type)), collapse = ", ")))

# ==============================================================================
# Per-cluster protein marker heatmap (if TotalA columns exist)
# ==============================================================================
protein_cols <- intersect(
  c("CD3_TotalA", "CD4_TotalA", "CD8A_TotalA",
    "CD56_TotalA", "CD16_TotalA", "CD25_TotalA", "CD127_TotalA",
    "CD45RA_TotalA", "CD45RO_TotalA", "CD19_TotalA", "CD14_TotalA",
    "TIGIT_TotalA", "CD279_TotalA"),
  colnames(obj@meta.data)
)

if (length(protein_cols) > 0) {
  message("Generating protein marker heatmap per cell type...")
  md <- obj@meta.data
  md$cell_type <- as.character(md$cell_type)

  cluster_means <- md %>%
    group_by(cell_type) %>%
    summarise(across(all_of(protein_cols), mean, na.rm = TRUE),
              n_cells = n(), .groups = "drop")

  write.csv(cluster_means,
            file.path(OUT_DIR, "Tcells_celltype_protein_means.csv"),
            row.names = FALSE)

  mat <- cluster_means %>%
    select(cell_type, all_of(protein_cols)) %>%
    pivot_longer(-cell_type, names_to = "marker", values_to = "mean_expr") %>%
    group_by(marker) %>%
    mutate(scaled = rescale(mean_expr)) %>%
    ungroup()
  mat$marker <- gsub("_TotalA", "", mat$marker)

  p_heat <- ggplot(mat, aes(x = marker, y = cell_type, fill = scaled)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#D73027",
                         midpoint = 0.5, name = "Scaled\nexpression") +
    labs(x = NULL, y = "Cell type",
         title = "Protein marker expression per cell type") +
    theme_classic(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 7),
          plot.title  = element_text(face = "bold", hjust = 0.5))

  ggsave(file.path(OUT_DIR, "Tcells_celltype_marker_heatmap.pdf"),
         p_heat, width = 8, height = 5)
  ggsave(file.path(OUT_DIR, "Tcells_celltype_marker_heatmap.png"),
         p_heat, width = 8, height = 5, dpi = 300)
  message("Saved: Tcells_celltype_marker_heatmap")
} else {
  message("No TotalA protein columns found — skipping marker heatmap.")
}

# ==============================================================================
# Save UMAP cache
# ==============================================================================
emb    <- Embeddings(obj, reduction = "umap")
umap_df <- as.data.frame(emb)
colnames(umap_df) <- c("UMAP_1", "UMAP_2")
umap_df$label     <- as.character(obj@meta.data$cell_type)
umap_df$condition <- as.character(obj@meta.data$deMultliplex2_final_mapped)
umap_df$Tissue    <- "Tcells"

write.csv(umap_df, file.path(CACHE_DIR, "Tcells_umap.csv"))
message("Saved: Tcells_umap.csv")

rm(obj); gc()

# ==============================================================================
# UMAP plot — coloured by cell type
# ==============================================================================
p_umap <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = label)) +
  geom_point(size = 0.3, alpha = 0.6) +
  guides(color = guide_legend(override.aes = list(size = 3), ncol = 1)) +
  labs(title = "D  T cells — Harmony LSI UMAP",
       color = NULL, x = "UMAP 1", y = "UMAP 2") +
  theme_classic(base_size = 9) +
  theme(axis.text  = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 7),
        plot.title  = element_text(face = "bold", hjust = 0.5))

ggsave(file.path(OUT_DIR, "Fig5D_Tcells_UMAP.pdf"),  p_umap, width = 5, height = 4)
ggsave(file.path(OUT_DIR, "Fig5D_Tcells_UMAP.png"),  p_umap, width = 5, height = 4, dpi = 300)
message("Saved: Fig5D_Tcells_UMAP")

# ==============================================================================
# UMAP plot — coloured by HTO condition (for QC check)
# ==============================================================================
p_cond <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = condition)) +
  geom_point(size = 0.3, alpha = 0.6) +
  guides(color = guide_legend(override.aes = list(size = 3), ncol = 1)) +
  labs(title = "T cells — UMAP by condition (QC)",
       color = NULL, x = "UMAP 1", y = "UMAP 2") +
  theme_classic(base_size = 9) +
  theme(axis.text  = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 7),
        plot.title  = element_text(face = "bold", hjust = 0.5))

ggsave(file.path(OUT_DIR, "Fig5D_Tcells_UMAP_by_condition.pdf"), p_cond, width = 5, height = 4)
ggsave(file.path(OUT_DIR, "Fig5D_Tcells_UMAP_by_condition.png"), p_cond, width = 5, height = 4, dpi = 300)
message("Saved: Fig5D_Tcells_UMAP_by_condition (QC check)")

message("\nDone.")

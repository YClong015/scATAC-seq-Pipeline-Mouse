library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(ggrepel)

OUT_DIR  <- "/QRISdata/Q8448/Mouse_disease_data/QC_figures"
CACHE_DIR <- file.path(OUT_DIR, "umap_cache")
dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)

## Tissue configuration
tissue.configs <- list(
  Kidney = list(
    rds       = "/QRISdata/Q8448/Mouse_disease_data/Kidney/kidney_merged_universal.rds",
    reduction = "umap",
    label_col = "cell_type"
  ),
  Lung = list(
    rds       = "/QRISdata/Q8448/Mouse_disease_data/Lung/lung_universal_new_pruned.rds",
    reduction = "umap",
    label_col = "cell_type"
  ),
  Aorta = list(
    rds       = "/QRISdata/Q8448/Mouse_disease_data/Aorta/Aorta_integrated_universal.rds",
    reduction = "umap",
    label_col = "cell_type"
  ),
  Tcells = list(
    rds       = "/QRISdata/Q8448/Mouse_disease_data/Tcells/tcells_final_annotated.rds",
    reduction = "umap",
    label_col = "cell_type"
  )
)

## Phase 1: Load RDS and cache UMAP coordinates (skip if cache exists)
for (tissue in names(tissue.configs)) {
  cache.file <- file.path(CACHE_DIR, paste0(tissue, "_umap.csv"))
  if (file.exists(cache.file)) {
    message("Cache found, skipping RDS load: ", tissue)
    next
  }
  cfg <- tissue.configs[[tissue]]
  message("Loading RDS: ", tissue)
  obj <- readRDS(cfg$rds)

  # Extract UMAP coordinates
  emb <- Embeddings(obj, reduction = cfg$reduction)
  df  <- as.data.frame(emb)
  colnames(df) <- c("UMAP_1", "UMAP_2")

  # Extract label
  df$label  <- as.character(obj@meta.data[[cfg$label_col]])
  df$Tissue <- tissue

  write.csv(df, cache.file, row.names = TRUE)
  message("  Cache saved: ", cache.file)
  rm(obj); gc()
}

## Phase 2: Read caches and plot
umap.list <- list()
for (tissue in names(tissue.configs)) {
  cache.file <- file.path(CACHE_DIR, paste0(tissue, "_umap.csv"))
  df <- read.csv(cache.file, row.names = 1, check.names = FALSE)
  df$label  <- as.character(df$label)   # ensure discrete color mapping
  df$Tissue <- tissue
  umap.list[[tissue]] <- df
}

## Saturated qualitative color palette (up to 20 cell types)
umap.palette <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00", "#984EA3",
  "#A65628", "#F781BF", "#8DD3C7", "#BEBADA", "#FB8072",
  "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#BC80BD",
  "#CCEBC5", "#FFED6F", "#1F78B4", "#33A02C", "#FB9A99"
)

## Publication theme
pub.theme <- theme_classic(base_size = 9) +
  theme(
    plot.title      = element_text(face = "bold", size = 10, hjust = 0.5),
    axis.text       = element_blank(),
    axis.ticks      = element_blank(),
    axis.line       = element_line(linewidth = 0.4),
    axis.title      = element_text(size = 8),
    panel.grid      = element_blank(),
    legend.text     = element_text(size = 7),
    legend.title    = element_blank(),
    legend.key.size = unit(0.3, "cm"),
    legend.position = "right"
  )

panel.labels <- c("A", "B", "C", "D")

# Adaptive point size: smaller for dense tissues, larger for sparse
tissue.pt.size <- c(Kidney = 0.15, Lung = 0.5, Aorta = 0.5, Tcells = 0.4)

## Helper: compute centroids per tissue
make.centroids <- function(df) {
  df %>%
    group_by(label) %>%
    summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = "drop")
}

## Build panels - with and without centroid labels
plots.labeled   <- list()
plots.nolabel   <- list()

for (i in seq_along(names(tissue.configs))) {
  tissue <- names(tissue.configs)[i]
  df     <- umap.list[[tissue]]

  set.seed(1234)
  df <- df[sample(nrow(df)), ]

  n.labels <- length(unique(df$label))
  colors   <- umap.palette[seq_len(n.labels)]
  pt.size  <- tissue.pt.size[tissue]

  # Base plot (no labels)
  p.base <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = label)) +
    geom_point(size = pt.size, alpha = 0.85, stroke = 0) +
    scale_color_manual(values = colors) +
    labs(
      x     = "UMAP 1",
      y     = "UMAP 2",
      title = paste0(panel.labels[i], "  ", tissue)
    ) +
    guides(color = guide_legend(
      override.aes = list(size = 2.5, alpha = 1),
      ncol = 1
    )) +
    pub.theme

  plots.nolabel[[tissue]] <- p.base

  centroids <- make.centroids(df)
  p.labeled <- p.base +
    ggrepel::geom_label_repel(
      data         = centroids,
      aes(label    = label),
      size         = 2.5,
      fontface     = "bold",
      color        = "black",
      fill         = alpha("white", 0.7),
      label.size   = 0.2,
      box.padding  = 0.4,
      max.overlaps = 30,
      show.legend  = FALSE
    )

  plots.labeled[[tissue]] <- p.labeled
}

## Combine and save - labeled version
combined.labeled <- wrap_plots(plots.labeled, nrow = 2, ncol = 2)

ggsave(file.path(OUT_DIR, "Fig5_UMAP_annotation_labeled.pdf"),
       combined.labeled, width = 10, height = 8, units = "in")
ggsave(file.path(OUT_DIR, "Fig5_UMAP_annotation_labeled.png"),
       combined.labeled, width = 10, height = 8, units = "in", dpi = 300)
message("Saved: Fig5_UMAP_annotation_labeled (.pdf + .png)")

## Combine and save - no-label version (legend only)
combined.nolabel <- wrap_plots(plots.nolabel, nrow = 2, ncol = 2)

ggsave(file.path(OUT_DIR, "Fig5_UMAP_annotation_nolabel.pdf"),
       combined.nolabel, width = 10, height = 8, units = "in")
ggsave(file.path(OUT_DIR, "Fig5_UMAP_annotation_nolabel.png"),
       combined.nolabel, width = 10, height = 8, units = "in", dpi = 300)
message("Saved: Fig5_UMAP_annotation_nolabel (.pdf + .png)")

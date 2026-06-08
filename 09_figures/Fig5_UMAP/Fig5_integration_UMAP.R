#!/usr/bin/env Rscript
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(ggrepel)

RDS      <- "/QRISdata/Q8448/Mouse_disease_data/Integrated/All_Tissues_Integrated_.rds"
OUT_DIR  <- "/QRISdata/Q8448/Mouse_disease_data/QC_figures"
# Same cache directory used by Fig5_UMAP_annotation.R
CACHE_DIR <- file.path(OUT_DIR, "umap_cache")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

message("Loading integrated object...")
obj <- readRDS(RDS)

# cell_type from the Fig5_UMAP_annotation.R caches so labels match the per-tissue
# UMAPs; T cells are shown separately (Fig5_Tcells_annotation.R) and excluded here
tissues <- c("Kidney", "Aorta", "Lung")

label.lookup <- lapply(tissues, function(tis) {
  f <- file.path(CACHE_DIR, paste0(tis, "_umap.csv"))
  if (!file.exists(f)) stop("Missing annotation cache: ", f,
                            "\nRun Fig5_UMAP_annotation.R first.")
  df <- read.csv(f, row.names = 1, check.names = FALSE)
  data.frame(
    cell      = paste0(tis, "_", rownames(df)),
    cell_type = as.character(df$label),
    stringsAsFactors = FALSE
  )
})
label.df  <- do.call(rbind, label.lookup)
label.map <- setNames(label.df$cell_type, label.df$cell)

message("Label lookup built: ", nrow(label.df), " cells across ", length(tissues), " tissues")

## Color palettes
tissue.colors <- c(
  Kidney = "#E41A1C",
  Lung   = "#377EB8",
  Aorta  = "#4DAF4A"
)

celltype.palette <- c(
  "#E41A1C","#377EB8","#4DAF4A","#FF7F00","#984EA3",
  "#A65628","#F781BF","#8DD3C7","#BEBADA","#FB8072",
  "#80B1D3","#FDB462","#B3DE69","#FCCDE5","#BC80BD",
  "#CCEBC5","#FFED6F","#1F78B4","#33A02C","#FB9A99",
  "#6A3D9A","#B15928","#CAB2D6","#FFFF99","#B2DF8A",
  "#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854"
)

## Publication theme
pub.theme <- theme_classic(base_size = 11) +
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

# Extract UMAP coordinates; cell_type comes from the per-tissue lookup so Panel C
# matches Fig5_UMAP_annotation.R.
get.umap.df <- function(obj, reduction) {
  emb <- Embeddings(obj, reduction = reduction)
  df  <- as.data.frame(emb)
  colnames(df) <- c("UMAP_1", "UMAP_2")
  df$Tissue <- as.character(obj$Tissue)
  df$cell   <- rownames(df)

  # T cells are shown separately (Fig5_Tcells_annotation.R); exclude here
  df <- df[df$Tissue != "Tcells", ]

  # Map cached labels; fall back to seurat_clusters for any unmatched cells
  df$cell_type <- label.map[df$cell]
  n.missing <- sum(is.na(df$cell_type))
  if (n.missing > 0) {
    message("  ", n.missing, " cells without a cached label in ", reduction,
            " - falling back to seurat_clusters")
    missing <- is.na(df$cell_type)
    df$cell_type[missing] <- paste0("C", obj$seurat_clusters[df$cell[missing]])
  }

  set.seed(1234)
  df[sample(nrow(df)), ]
}

df.lsi     <- get.umap.df(obj, "umap.lsi")
df.harmony <- get.umap.df(obj, "umap.harmony")

## Panel A: Before Harmony - coloured by Tissue
p.before <- ggplot(df.lsi, aes(x = UMAP_1, y = UMAP_2, color = Tissue)) +
  geom_point(size = 0.3, alpha = 0.6, stroke = 0) +
  scale_color_manual(values = tissue.colors) +
  guides(color = guide_legend(
    override.aes = list(size = 4, alpha = 1),
    ncol = 1, byrow = TRUE)) +
  ggtitle("Before integration") +
  pub.theme

## Panel B: After Harmony - coloured by Tissue
p.after.tissue <- ggplot(df.harmony, aes(x = UMAP_1, y = UMAP_2, color = Tissue)) +
  geom_point(size = 0.3, alpha = 0.6, stroke = 0) +
  scale_color_manual(values = tissue.colors) +
  guides(color = guide_legend(
    override.aes = list(size = 4, alpha = 1),
    ncol = 1, byrow = TRUE)) +
  ggtitle("After integration (tissue)") +
  pub.theme

## Panel C: After Harmony - coloured by Cell Type (from Fig5 annotations)
cell.types <- sort(unique(df.harmony$cell_type))
n.ct       <- length(cell.types)
ct.colors  <- setNames(celltype.palette[seq_len(n.ct)], cell.types)

centroids <- df.harmony %>%
  group_by(cell_type) %>%
  summarise(UMAP_1 = median(UMAP_1),
            UMAP_2 = median(UMAP_2),
            .groups = "drop")

p.after.ct <- ggplot(df.harmony, aes(x = UMAP_1, y = UMAP_2, color = cell_type)) +
  geom_point(size = 0.3, alpha = 0.6, stroke = 0) +
  scale_color_manual(values = ct.colors) +
  ggrepel::geom_label_repel(
    data          = centroids,
    aes(label     = cell_type),
    size          = 2.8,
    fontface      = "bold",
    color         = "black",
    fill          = alpha("white", 0.82),
    label.size    = 0.2,
    label.padding = unit(0.18, "lines"),
    box.padding   = 0.4,
    max.overlaps  = 50,
    min.segment.length = 0.2,
    segment.size  = 0.35,
    show.legend   = FALSE
  ) +
  ggtitle("After integration (cell type)") +
  guides(color = guide_legend(
    override.aes = list(size = 4, alpha = 1),
    ncol = 2, byrow = TRUE)) +
  pub.theme +
  theme(legend.position = "right")

## Combine: A | B | C
combined <- p.before + p.after.tissue + p.after.ct +
  plot_layout(nrow = 1, widths = c(1, 1, 1.5)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(size = 13, face = "bold")
    )
  )

## Save
ggsave(file.path(OUT_DIR, "Integration_A_before.pdf"),
       p.before, width = 4.5, height = 4, useDingbats = FALSE)
ggsave(file.path(OUT_DIR, "Integration_B_after_tissue.pdf"),
       p.after.tissue, width = 4.5, height = 4, useDingbats = FALSE)
ggsave(file.path(OUT_DIR, "Integration_C_after_celltype.pdf"),
       p.after.ct, width = 6.5, height = 4, useDingbats = FALSE)

ggsave(file.path(OUT_DIR, "Fig_Integration_UMAP.pdf"),
       combined, width = 14, height = 4.8, useDingbats = FALSE)
ggsave(file.path(OUT_DIR, "Fig_Integration_UMAP.png"),
       combined, width = 14, height = 4.8, dpi = 300)

message("Saved: Fig_Integration_UMAP (.pdf + .png)")
message("Saved: Individual panels (A, B, C)")
message("Done.")

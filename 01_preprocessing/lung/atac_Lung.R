## Lung scATAC-seq workflow
# QC -> doublet removal -> per-sample LSI -> Harmony -> UMAP -> annotation -> split fragments

set.seed(1234)
library(Signac)
library(Seurat)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(patchwork)
library(scDblFinder)
library(harmony)
library(Matrix)
library(BiocGenerics)

## Genome annotation (mm10)
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

## Per-sample loading + QC + doublet filtering
samples <- data.frame(
  folder = c("CL100168054_L01", "CL100167942_L01",
             "CL100168054_L02", "CL100167942_L02",
             "CL100168078_L02", "CL100168078_L01"),
  sample_id = c("Control_F2", "Control_M1",
                "Case_F3", "Case_M2", "Case_F1", "Case_M3"),
  group = c("Control", "Control", "Case", "Case", "Case", "Case"),
  sex = c("Female", "Male", "Female", "Male", "Female", "Male")
)

base_dir <- "/QRISdata/Q8448/Mouse_disease_data/Lung/Lung_cellatac"
obj_list <- list()

for (i in 1:nrow(samples)) {
  cur_folder <- samples$folder[i]
  cur_id <- samples$sample_id[i]
  print(paste0("[", i, "/", nrow(samples), "] Processing: ", cur_id))

  matrix_dir <- file.path(base_dir, cur_folder, "outs", "filter_peak_matrix")
  meta_path <- file.path(base_dir, cur_folder, "outs", "singlecell.csv")
  frag_path <- file.path(base_dir, cur_folder, "outs", "fragments.tsv.gz")

  counts_mat <- readMM(file = file.path(matrix_dir, "matrix.mtx.gz"))
  barcodes <- read.table(file = file.path(matrix_dir, "barcodes.tsv.gz"), header = FALSE)
  peaks <- read.table(file = file.path(matrix_dir, "peaks.bed.gz"), header = FALSE)

  colnames(counts_mat) <- barcodes$V1
  peak_names <- paste0(peaks$V1, ":", peaks$V2, "-", peaks$V3)
  rownames(counts_mat) <- make.unique(peak_names)

  if (file.exists(meta_path)) {
    metadata <- read.csv(meta_path, header = TRUE, row.names = 1)
  } else {
    metadata <- NULL
  }

  chrom_assay <- CreateChromatinAssay(
    counts = counts_mat,
    sep = c(":", "-"),
    genome = 'mm10',
    fragments = frag_path,
    min.cells = 10,
    min.features = 200,
    annotation = annotations
  )

  cur_obj <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "peaks",
    meta.data = metadata,
    project = cur_id
  )

  cur_obj$SampleID <- cur_id
  cur_obj$Group <- samples$group[i]
  cur_obj$Sex <- samples$sex[i]
  cur_obj$Dataset <- "Lung"

  cur_obj <- NucleosomeSignal(object = cur_obj)
  cur_obj <- TSSEnrichment(object = cur_obj, fast = FALSE)
  cur_obj$peak_region_fragments <- cur_obj$nCount_peaks

  mtx <- GetAssayData(cur_obj, layer = "counts")
  dbl_results <- scDblFinder(mtx, verbose = FALSE)
  cur_obj$scDblFinder.class <- dbl_results$scDblFinder.class
  cur_obj$scDblFinder.score <- dbl_results$scDblFinder.score

  # Filter before merging (paper thresholds: frags > 1000, TSS > 4)
  cur_obj <- subset(
    x = cur_obj,
    subset = nCount_peaks > 1000 &
             nCount_peaks < 100000 &
             TSS.enrichment > 4 &
             nucleosome_signal < 4 &
             scDblFinder.class == "singlet"
  )
  print(paste("   ", cur_id, "cells after filter:", ncol(cur_obj)))

  obj_list[[i]] <- cur_obj
}

## Merge
lung <- merge(
  x = obj_list[[1]],
  y = obj_list[2:length(obj_list)],
  add.cell.ids = samples$sample_id
)
rm(obj_list); gc()
print(dim(lung))

VlnPlot(lung, features = c('nCount_peaks', 'TSS.enrichment', 'nucleosome_signal'),
        ncol = 3, pt.size = 0, group.by = "SampleID")

## LSI + Harmony batch correction
lung <- RunTFIDF(lung)
lung <- FindTopFeatures(lung, min.cutoff = 'q0')
lung <- RunSVD(lung)
DepthCor(lung)

lsi_embeddings <- Embeddings(lung, reduction = "lsi")
harmony_embeddings <- HarmonyMatrix(
  data_mat = lsi_embeddings,
  meta_data = lung@meta.data,
  vars_use = "SampleID",
  do_pca = FALSE,
  verbose = TRUE
)
lung[["harmony"]] <- CreateDimReducObject(
  embeddings = harmony_embeddings,
  key = "harmony_",
  assay = "peaks"
)

## Clustering + UMAP (resolution 0.8, matches paper)
lung <- RunUMAP(object = lung, reduction = 'harmony', dims = 2:30, seed.use = 1234)
lung <- FindNeighbors(object = lung, reduction = 'harmony', dims = 2:30)
lung <- FindClusters(object = lung, verbose = FALSE, algorithm = 3, resolution = 0.8, random.seed = 1234)
DimPlot(object = lung, label = TRUE, pt.size = 0.5) + ggtitle("UMAP after Harmony")

## Gene activities
gene.activites <- GeneActivity(lung)
lung[['RNA']] <- CreateAssayObject(counts = gene.activites)
lung <- NormalizeData(
  object = lung,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(lung$nCount_RNA)
)

## Marker dotplot (Figure S5)
marker_list_lung <- list(
  "AT1" = c("Ager"),
  "AT2" = c("Sftpb"),
  "Ciliated" = c("Foxj1"),
  "Fibroblast" = c("Gucy1a1", "Speg", "Cntn5"),
  "Endothelial" = c("Vwf", "Clec14a", "Dpep1"),
  "Alveolar Macrophage" = c("Mcemp1"),
  "Interstitial Macrophage" = c("C1qb"),
  "Monocytes" = c("Cd14"),
  "Dendritic Cells" = c("Ccl17", "Itgae"),
  "B Cells" = c("Cd79b"),
  "Plasma/B-Jchain" = c("Jchain"),
  "T Cells (General)" = c("Cd3g"),
  "CD4 T Cells" = c("Cd4"),
  "CD8 T Cells" = c("Cd8a"),
  "NK Cells" = c("Nkg7")
)

DefaultAssay(lung) <- 'RNA'
DotPlot(lung, features = marker_list_lung,
        cols = c("lightgrey", "red"), dot.scale = 8) +
  RotatedAxis() + ggtitle("Cell Type Markers (Figure S5)") +
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10))

FeaturePlot(lung, features = c("Sftpb", "Mcemp1", "Cd79b", "Cd3g", "Vwf", "Gucy1a1"),
            min.cutoff = "q9", ncol = 3)

## Cluster markers
all.markers <- FindAllMarkers(lung, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

library(dplyr)
top5_markers <- all.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)
print(as.data.frame(top5_markers))

my_object <- ScaleData(lung, features = top5_markers$gene)
DoHeatmap(my_object, features = top5_markers$gene) + NoLegend() +
  ggtitle("Top 5 Markers per Cluster") +
  theme(axis.text.y = element_text(size = 5))

## Manual cluster checks
genes_AT1       <- c("Ager", "Hopx", "Pdpn", "Aqp5")
genes_Ciliated  <- c("Foxj1", "Tubb4b")
genes_T_cell    <- c("Trac", "Cd3d", "Cd3e")
genes_ILC2      <- c("Il7r", "Gata3", "Rora", "Areg")
genes_Eosinophil<- c("Epx", "Prg2", "Siglecf")
genes_Mac       <- c("Lyz2", "Adgre1", "Marco", "Mrc1")
genes_Epi       <- c("Epcam", "Krt8", "Krt18", "Krt19")

# Epithelial (AT1 vs Ciliated)
print(FeaturePlot(lung, features = c(genes_AT1, genes_Ciliated),
                  cols = c("lightgrey", "red"), ncol = 3) + ggtitle("AT1 & Ciliated Markers"))
print(VlnPlot(lung, features = c(genes_AT1, genes_Ciliated), stack = TRUE, flip = TRUE) + NoLegend())

# Immune (Eosinophil vs ILC2 vs T)
print(FeaturePlot(lung, features = c(genes_Eosinophil, genes_ILC2, genes_T_cell),
                  cols = c("lightgrey", "red"), ncol = 4))
print(VlnPlot(lung, features = c(genes_Eosinophil, genes_ILC2, genes_T_cell), stack = TRUE, flip = TRUE) + NoLegend())

# Macrophage subtypes
print(FeaturePlot(lung, features = genes_Mac, cols = c("lightgrey", "red"), ncol = 2))
print(VlnPlot(lung, features = genes_Mac, stack = TRUE, flip = TRUE) + NoLegend())

# General epithelial
FeaturePlot(lung, features = genes_Epi, cols = c("lightgrey", "red"), ncol = 2)

## Validation: refine ambiguous clusters
# Cluster 4: T cell vs ILC2 (T = Cd3d/Cd3e/Trac high; ILC2 = Gata3/Il7r/Areg high, Cd3 negative)
genes_T     <- c("Trac", "Cd3d", "Cd3e")
genes_ILC2  <- c("Il7r", "Gata3", "Rora", "Areg")
print(VlnPlot(lung, features = c(genes_T, genes_ILC2), stack = TRUE, flip = TRUE, pt.size = 0) +
      ggtitle("Cluster 4 Check: T vs ILC2"))
print(FeaturePlot(lung, features = c("Cd3d", "Gata3", "Areg"), ncol = 3, cols = c("lightgrey", "red")))

# Cluster 5: Eosinophil vs Macrophage (Eos = Epx/Prg2/Ear1/Ear2/Alox15; exclude Adgre1/Lyz2/Csf1r)
genes_Eos <- c("Epx", "Prg2", "Ear1", "Ear2", "Alox15")
genes_Mac_Exclude <- c("Adgre1", "Lyz2", "Csf1r")
print(VlnPlot(lung, features = c(genes_Eos, genes_Mac_Exclude), stack = TRUE, flip = TRUE, pt.size = 0) +
      ggtitle("Cluster 5 Check: Eosinophil Identity"))
print(FeaturePlot(lung, features = c("Epx", "Adgre1", "Prg2"), ncol = 3, cols = c("lightgrey", "red")))

# Cluster 12: Endothelial (Pecam1/Kdr/Emcn/Cldn5)
genes_EC <- c("Pecam1", "Kdr", "Emcn", "Cldn5")
print(VlnPlot(lung, features = genes_EC, stack = TRUE, flip = TRUE, pt.size = 0) +
      ggtitle("Cluster 12 Check: Is it EC?"))
print(FeaturePlot(lung, features = genes_EC, ncol = 2, cols = c("lightgrey", "red")))

## Cluster annotation
new_cluster_ids <- c(
  "0"  = "Mac-alv",
  "1"  = "EC-vasc",
  "2"  = "Mo-Ly6c+",
  "3"  = "Mac-inter",
  "4"  = "B",
  "5"  = "Eosinophils",
  "6"  = "T",
  "7"  = "NK",
  "8"  = "Fib",
  "9"  = "Ciliated",
  "10" = "Fib",
  "11" = "Fib",
  "12" = "SMCs",
  "13" = "Pen",
  "14" = "AT2",
  "15" = "Fib",
  "16" = "Mesothelial",
  "17" = "Low Quality",
  "18" = "NK",
  "19" = "Low Quality"
)

lung <- RenameIdents(lung, new_cluster_ids)
lung$cell_type <- Idents(lung)
print(table(lung$cell_type))

# Drop Low Quality clusters
lung_clean <- subset(lung, idents = "Low Quality", invert = TRUE)
lung_clean$cell_type <- droplevels(lung_clean$cell_type)
Idents(lung_clean) <- lung_clean$cell_type
print(table(lung_clean$cell_type))

out_path <- "/QRISdata/Q8448/Mouse_disease_data/Lung/lung_integrated_clean_annotated.rds"
saveRDS(lung_clean, file = out_path)
message(paste0("Saved cleaned object to: ", out_path))

DimPlot(lung, label = TRUE, pt.size = 0.5, seed = 1234) + NoLegend()

## Split fragments by cell type
output_dir <- "/QRISdata/Q8448/Mouse_disease_data/Lung/Lung_fragments_file"
dir.create(output_dir, showWarnings = FALSE)

DefaultAssay(lung) <- 'peaks'
SplitFragments(
  object = lung,
  assay = 'peaks',
  group.by = 'cell_type',
  outdir = output_dir,
  verbose = TRUE
)

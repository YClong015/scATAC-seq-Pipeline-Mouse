## Aortic scATAC-seq Integration Workflow (Reproducible with Seed)
library(Signac)
library(Seurat)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(patchwork)
library(scDblFinder)
library(dplyr)
library(GenomeInfoDb)
library(GenomicRanges)
library(BiocParallel)

## Reproducibility (GLOBAL SEED)
SEED <- 1234L
set.seed(SEED)

# scDblFinder reproducibility: SerialParam(RNGseed=SEED) forces deterministic RNG
# under BiocParallel.
bp <- SerialParam(RNGseed = SEED)

## Annotation (mm10)
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "mm10"

## Sample Info
base_dir <- "/QRISdata/Q8448/Mouse_disease_data/Aorta/Aorta_cellranger_atac"

samples <- data.frame(
  sample_id = c(
    "SRR21686724_Aorta_atac",
    "SRR21686722_Aorta_atac"
  ),
  group = c("Control", "Challenge"),
  stringsAsFactors = FALSE
)

## Build Unified Peak Set
print("Building Unified Peak Set...")
peaks_list <- list()

for (i in 1:nrow(samples)) {
  sample_name <- samples$sample_id[i]
  peak_file <- file.path(base_dir, sample_name, "outs", "peaks.bed")
  
  peaks_df <- read.table(peak_file, col.names = c("chr", "start", "end"))
  peaks_list[[sample_name]] <- makeGRangesFromDataFrame(peaks_df)
}

combined.peaks <- reduce(x = c(peaks_list[[1]], peaks_list[[2]]))

peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths < 10000 & peakwidths > 20]

print(paste0("Unified Peak Set constructed: ", length(combined.peaks), " peaks."))

## Loop: Re-quantify -> QC -> Per-sample LSI
obj_list <- list()

for (i in 1:nrow(samples)) {
  
  sample_name <- samples$sample_id[i]
  group_name <- samples$group[i]
  print(paste0("Processing sample: ", sample_name, "..."))
  
  data_dir <- file.path(base_dir, sample_name, "outs")
  
  frag_path <- file.path(data_dir, "fragments.tsv.gz")
  
  meta <- read.csv(
    file.path(data_dir, "singlecell.csv"),
    header = TRUE,
    row.names = 1
  )
  
  cells_use <- rownames(meta)
  frag.object <- CreateFragmentObject(path = frag_path, cells = cells_use)
  
  print("  - Re-quantifying counts on unified peaks...")
  counts <- FeatureMatrix(
    fragments = frag.object,
    features = combined.peaks,
    cells = cells_use
  )
  
  chrom_assay <- CreateChromatinAssay(
    counts = counts,
    sep = c(":", "-"),
    genome = "mm10",
    fragments = frag.object,
    min.cells = 10,
    min.features = 200
  )
  
  cur_obj <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "peaks",
    meta.data = meta
  )
  
  cur_obj$Group <- group_name
  cur_obj$sample_id <- sample_name
  
  Annotation(cur_obj) <- annotations
  
  cur_obj <- NucleosomeSignal(cur_obj)
  cur_obj <- TSSEnrichment(cur_obj, fast = FALSE)
  
  if (!"peak_region_fragments" %in% colnames(cur_obj@meta.data)) {
    cur_obj$peak_region_fragments <- cur_obj$nCount_peaks
  }
  
  if ("passed_filters" %in% colnames(cur_obj@meta.data)) {
    cur_obj$pct_read_in_peaks <- cur_obj$peak_region_fragments /
      cur_obj$passed_filters * 100
  } else {
    cur_obj$pct_read_in_peaks <- NA
  }
  
  cur_obj$blacklist_ratio <- cur_obj$blacklist_region_fragments /
    cur_obj$peak_region_fragments
  
  ## Doublet Removal (Reproducible)
  counts_mat <- GetAssayData(
    cur_obj,
    assay = "peaks",
    layer = "counts"
  )
  
  # Ensure deterministic artificial doublet generation
  set.seed(SEED)
  dbl_results <- scDblFinder(counts_mat, BPPARAM = bp)
  
  cur_obj$scDblFinder.class <- dbl_results$scDblFinder.class
  
  print("  - Filtering...")
  cur_obj <- subset(
    cur_obj,
    subset = scDblFinder.class == "singlet" &
      peak_region_fragments > 3000 &
      peak_region_fragments < 100000 &
      pct_read_in_peaks > 15 &
      blacklist_ratio < 0.025 &
      nucleosome_signal < 10 &
      TSS.enrichment > 2
  )
  
  ## Per-sample LSI (Reproducible)
  # RunSVD uses irlba internally; set.seed helps keep it stable across runs.
  set.seed(SEED)
  cur_obj <- RunTFIDF(cur_obj)
  cur_obj <- FindTopFeatures(cur_obj, min.cutoff = "q0")
  cur_obj <- RunSVD(cur_obj)
  
  obj_list[[sample_name]] <- cur_obj
}

## Integration Anchors (rLSI)
common_features <- rownames(obj_list[[1]])
print(paste0("Common features for integration: ", length(common_features)))

# Not strictly necessary, but harmless to keep the whole pipeline deterministic
set.seed(SEED)
integration.anchors <- FindIntegrationAnchors(
  object.list = obj_list,
  anchor.features = common_features,
  reduction = "rlsi",
  dims = 2:30
)

print("Merging objects to create aortic object...")
aortic <- merge(
  x = obj_list[[1]],
  y = obj_list[-1]
)

DefaultAssay(aortic) <- "peaks"
set.seed(SEED)
aortic <- RunTFIDF(aortic)
aortic <- FindTopFeatures(aortic, min.cutoff = "q0")
aortic <- RunSVD(aortic)

print("Integrating Embeddings (Anchors -> integrated LSI)...")
aortic <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = aortic[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:30
)

## UMAP + Clustering (Reproducible)

aortic <- RunUMAP(
  aortic,
  reduction = "integrated_lsi",
  dims = 2:30,
  seed.use = SEED
)

aortic <- FindNeighbors(
  aortic,
  reduction = "integrated_lsi",
  dims = 2:30
)

aortic <- FindClusters(
  aortic,
  resolution = 0.6,
  random.seed = SEED,
  verbose = FALSE
)

## Gene Activities
DefaultAssay(aortic) <- "peaks"

gene.activities <- GeneActivity(
  aortic,
  extend.upstream = 2000,
  extend.downstream = 0
)

aortic[["RNA"]] <- CreateAssayObject(counts = gene.activities)

aortic <- NormalizeData(
  aortic,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = median(aortic$nCount_RNA)
)

marker_list_aortic <- list(
  "SMC (Contractile)" = c("Myh11", "Acta2", "Mir143hg", "Cnn1", "Tagln"),
  "Fibroblast" = c("Dcn", "Lum", "Pdgfra", "Serpinf1", "Col1a1"),
  "Pericyte" = c("Rgs5", "Cspg4", "Pdgfrb"),
  "Endothelial" = c("Pecam1", "Cdh5", "Egfl7", "Tie1"),
  "Macrophage" = c("Lyz2", "Adgre1", "C1qa", "C1qb", "Ptprc"),
  "T-cell" = c("Cd3d", "Cd3e", "Trac", "Cd247"),
  "Adaptive/Injury SMC" = c("Fn1", "Lox", "Eln", "Col3a1", "Spp1"),
  "Cycling" = c("Top2a", "Mki67", "Birc5")
)

DefaultAssay(aortic) <- "RNA"

p_markers <- DotPlot(
  object = aortic,
  features = marker_list_aortic,
  cols = c("lightgrey", "red"),
  dot.scale = 8
) +
  RotatedAxis() +
  ggtitle("Aortic Cell Type Markers (GeneActivity)") +
  theme(
    axis.text.x = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(size = 10)
  )

print(p_markers)

Idents(aortic) <- "seurat_clusters"

mk <- FindAllMarkers(
  aortic,
  only.pos = TRUE,
  min.pct = 0.2,
  logfc.threshold = 0.25
)

top5 <- mk %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5)

cluster_list <- split(top5, top5$cluster)
print(lapply(cluster_list, function(x) x[, c("gene", "avg_log2FC")]))

## Manual Annotation
print("Annotating Clusters...")

new_cluster_ids <- c(
  "0" = "SMC",
  "1" = "SMC",
  "2" = "Mac",
  "3" = "Fibroblast",
  "4" = "Fibroblast",
  "5" = "Fibroblast",
  "6" = "Mac",
  "7" = "Mac",
  "8" = "Pericyte",
  "9" = "Mac",
  "10" = "T-cell",
  "11" = "Endothelial"
)

aortic <- RenameIdents(aortic, new_cluster_ids)
aortic$cell_type <- Idents(aortic)

## Save annotated object (AFTER MARKING)
aortic@misc$pipeline_params <- list(
  seed = SEED,
  resolution = 0.6,
  umap_dims = 2:30,
  integration_dims = 2:30,
  geneactivity_upstream = 2000
)

save_dir <- paste0(
  "/QRISdata/Q8448/Mouse_disease_data/Aorta/",
  "objects"
)
dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)

obj_tag <- paste0(
  "aortic_integrated_res0.6_up2k_seed",
  SEED
)

saveRDS(
  aortic,
  file = file.path(save_dir, paste0(obj_tag, ".rds"))
)

print(paste0("Saved object to: ", file.path(save_dir, obj_tag)))

## Plots
p1 <- DimPlot(
  aortic,
  reduction = "umap",
  label = TRUE,
  pt.size = 0.5
) +
  ggtitle("Integrated Aortic Atlas (Unified Peaks + rLSI)")

p2 <- DimPlot(
  aortic,
  reduction = "umap",
  split.by = "Group",
  label = TRUE,
  pt.size = 0.5
)

print(p1)
print(p2)

## Split fragment files by cell type
output_dir <- paste0(
  "/QRISdata/Q8448/Mouse_disease_data/Aorta/",
  "fragment_files"
)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

DefaultAssay(aortic) <- "peaks"

print("Start splitting Fragment files by cell type...")
SplitFragments(
  object = aortic,
  assay = "peaks",
  group.by = "cell_type",
  outdir = output_dir,
  verbose = TRUE
)
print("Finished splitting Fragment files by cell type.")




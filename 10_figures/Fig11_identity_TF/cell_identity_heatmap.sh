#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --job-name=cell_identity_heatmap
#SBATCH --time=02:00:00
#SBATCH --partition=general
#SBATCH --account=a_imb_ccbcd
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

set -euo pipefail

module load r/4.4.2

Rscript - <<'REOF'

suppressPackageStartupMessages({ library(dplyr); library(tidyr) })

if (!requireNamespace("ComplexHeatmap", quietly=TRUE)) {
  if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager", repos="https://cloud.r-project.org")
  BiocManager::install("ComplexHeatmap", ask=FALSE)
}
if (!requireNamespace("circlize", quietly=TRUE))
  install.packages("circlize", repos="https://cloud.r-project.org")
suppressPackageStartupMessages({ library(ComplexHeatmap); library(circlize); library(grid) })

# ════════════════════════════════════════════════════════════════
# PATHS
# ════════════════════════════════════════════════════════════════
BASE_PB  <- "/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_pseudobulk_DESeq2"
BASE_DAR <- "/QRISdata/Q8448/Mouse_disease_data/DAR"
HOMER_DIRS <- list(
  "Aorta"  = file.path(BASE_PB, "DAR_pseudobulk_Aorta_DESeq2/HOMER_stable_bg"),
  "Lung"   = file.path(BASE_PB, "DAR_pseudobulk_Lung_DESeq2/HOMER_stable_bg"),
  "Kidney" = file.path(BASE_PB, "DAR_pseudobulk_Kidney_DESeq2/HOMER_stable_bg"),
  "Tcell"  = file.path(BASE_PB, "DAR_pseudobulk_Tcells_DESeq2/HOMER_stable_bg")
)
OUT_DIR <- file.path(BASE_DAR, "DAR_science_comparison/figures_Fig4G_style")
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

# ── Parameters ──────────────────────────────────────────────────
MIN_SIG_P       <- 0.05   # p-value for motif selection
MIN_LOG2FC      <- 0.5    # min log2FC for Cell Identity selection
CAP_LOG2FC      <- 3      # color scale cap (CM paper uses 2)
KIDNEY_CONTRAST <- "Day42_vs_Sham"

# ── AP-1 / CTCF patterns ────────────────────────────────────────
AP1_PATTERNS <- paste(c(
  "^Fra1\\(","^Fra2\\(","^Fos\\(","^FOS\\(","^c-Fos\\(",
  "^AP-1\\(","^Jun-AP1\\(","^Jun\\(","^JUN\\(",
  "^JunB\\(","^JunD\\(","^JUNB\\(","^JUND\\(",
  "^BATF\\(","^Fosl1\\(","^Fosl2\\(","^FOSL1\\(","^FOSL2\\(",
  "^Bach1\\(","^Bach2\\(","^BACH1\\(","^BACH2\\(",
  "^Atf3\\(","^ATF3\\("
), collapse="|")
CTCF_PATTERNS <- paste(c("^CTCF\\(","^BORIS\\(","^CTCF-Satellite"), collapse="|")

classify_motif <- function(m) {
  if (grepl(AP1_PATTERNS,  m, perl=TRUE)) return("1_AP-1")
  if (grepl(CTCF_PATTERNS, m, perl=TRUE)) return("2_CTCF")
  return("3_Cell Identity")
}

# ════════════════════════════════════════════════════════════════
# LOAD HOMER — only opening_vs_stable and closing_vs_stable
# color = log2(pct_target / pct_background), same as CM paper
# ════════════════════════════════════════════════════════════════
to_num      <- function(x) suppressWarnings(as.numeric(gsub(",|%","",as.character(x))))
short_motif <- function(x) trimws(sapply(strsplit(as.character(x),"/"),`[`,1))

all_rows <- list()
for (tissue in names(HOMER_DIRS)) {
  hdir <- HOMER_DIRS[[tissue]]
  if (!dir.exists(hdir)) { message("SKIP: ", hdir); next }
  fps <- list.files(hdir, pattern="^knownResults\\.txt$", recursive=TRUE, full.names=TRUE)
  message("Tissue: ", tissue, " — ", length(fps), " files")

  for (fp in fps) {
    dname     <- basename(dirname(fp))
    parts     <- strsplit(dname, "__", fixed=TRUE)[[1]]
    if (length(parts) < 4) next
    direction <- parts[length(parts)]
    if (!direction %in% c("opening_vs_stable","closing_vs_stable")) next
    if (tissue == "Kidney" && parts[2] != KIDNEY_CONTRAST) next

    tab <- tryCatch(read.delim(fp, header=TRUE, sep="\t",
                               stringsAsFactors=FALSE, check.names=FALSE),
                    error=function(e) NULL)
    if (is.null(tab) || nrow(tab) == 0) next
    if (!all(c("Motif Name","P-value") %in% colnames(tab))) next

    pct_t_col <- grep("% of Target",     colnames(tab), value=TRUE)[1]
    pct_b_col <- grep("% of Background", colnames(tab), value=TRUE)[1]
    if (is.na(pct_t_col) || is.na(pct_b_col)) next

    tab$motif     <- short_motif(tab[["Motif Name"]])
    tab$pval      <- to_num(tab[["P-value"]])
    pct_t         <- to_num(tab[[pct_t_col]])
    pct_b         <- to_num(tab[[pct_b_col]])
    # log2FC: same metric as CM paper (log2 of target% / background%)
    tab$log2fc    <- pmin(pmax(log2(pct_t / pmax(pct_b, 0.01)),
                               -CAP_LOG2FC), CAP_LOG2FC)
    tab$tissue    <- tissue
    tab$ct_key    <- paste0(tissue, "_", parts[1])
    tab$direction <- direction

    all_rows[[length(all_rows)+1]] <-
      tab[, c("motif","pval","log2fc","tissue","ct_key","direction")]
  }
}

df <- dplyr::bind_rows(all_rows) %>%
  group_by(motif, ct_key, direction) %>%
  arrange(pval) %>% slice(1) %>% ungroup()

message("\nTotal rows: ", nrow(df))
message("Cell types: ", paste(sort(unique(df$ct_key)), collapse=", "))

# ════════════════════════════════════════════════════════════════
# PLOT FUNCTION — per tissue
# Structure mirrors CM paper Figure 4G:
#   rows = TF motifs (AP-1 / CTCF / Cell Identity per cell type)
#   cols = Opening | Stable | Closing  (3 per cell type)
#   color = log2FC vs stable background
# ════════════════════════════════════════════════════════════════
plot_heatmap <- function(df_t, tissue_name) {
  message("\n=== ", tissue_name, " ===")
  cts      <- sort(unique(df_t$ct_key))
  ct_short <- gsub(paste0("^", tissue_name, "_"), "", cts)

  # ── Motif selection ──────────────────────────────────────────
  ap1_ctcf_motifs <- df_t %>%
    filter(grepl(AP1_PATTERNS, motif, perl=TRUE) |
           grepl(CTCF_PATTERNS, motif, perl=TRUE)) %>%
    pull(motif) %>% unique()

  identity_motifs <- df_t %>%
    filter(!grepl(AP1_PATTERNS, motif, perl=TRUE),
           !grepl(CTCF_PATTERNS, motif, perl=TRUE),
           pval < MIN_SIG_P, log2fc >= MIN_LOG2FC) %>%
    pull(motif) %>% unique()

  sig_motifs <- unique(c(ap1_ctcf_motifs, identity_motifs))
  message("  AP-1/CTCF: ", length(ap1_ctcf_motifs),
          "  Cell Identity: ", length(identity_motifs))
  if (length(sig_motifs) < 5) { message("  SKIP"); return(invisible(NULL)) }

  # ── Build matrix ─────────────────────────────────────────────
  # 3 columns per CT: Opening | Stable (=0, reference) | Closing
  df_open <- df_t %>%
    filter(direction=="opening_vs_stable", motif %in% sig_motifs) %>%
    mutate(col_id = paste0(ct_key, "__Opening")) %>%
    select(motif, col_id, log2fc)

  df_close <- df_t %>%
    filter(direction=="closing_vs_stable", motif %in% sig_motifs) %>%
    mutate(col_id = paste0(ct_key, "__Closing")) %>%
    select(motif, col_id, log2fc)

  # Stable = reference, all zero
  df_stable <- expand.grid(motif=sig_motifs, ct_key=cts,
                            stringsAsFactors=FALSE) %>%
    mutate(log2fc=0, col_id=paste0(ct_key, "__Stable")) %>%
    select(motif, col_id, log2fc)

  mat_wide <- bind_rows(df_open, df_stable, df_close) %>%
    pivot_wider(names_from=col_id, values_from=log2fc, values_fill=0)

  motifs <- mat_wide$motif
  mat    <- as.matrix(mat_wide[, -1])
  rownames(mat) <- motifs

  # Ordered columns: for each CT, [Opening, Stable, Closing]
  col_order <- unlist(lapply(cts, function(ct)
    c(paste0(ct,"__Opening"), paste0(ct,"__Stable"), paste0(ct,"__Closing"))
  ))
  mat <- mat[, col_order[col_order %in% colnames(mat)], drop=FALSE]

  col_type     <- gsub(".*__", "", colnames(mat))
  col_ct_short <- gsub(paste0("^",tissue_name,"_"),
                       "", gsub("__.*","", colnames(mat)))
  col_split    <- factor(col_ct_short, levels=unique(col_ct_short))

  # ── Column annotations ────────────────────────────────────────
  # Top bar: red=Opening, grey=Stable, blue=Closing (like CM paper)
  peak_colors <- c("Opening"="#D6604D", "Stable"="#AAAAAA", "Closing"="#4393C3")
  ha_top <- HeatmapAnnotation(
    "Peak" = col_type,
    col    = list("Peak" = peak_colors),
    annotation_name_gp = gpar(fontsize=8),
    simple_anno_size   = unit(4, "mm"),
    show_annotation_name = TRUE,
    show_legend = TRUE
  )

  # Column labels: show peak type name, colored (like CM paper axis labels)
  col_label_colors <- peak_colors[col_type]
  col_labels <- col_type  # "Opening", "Stable", "Closing"

  # ── Row grouping ─────────────────────────────────────────────
  # Assign each identity motif to the CT with highest closing log2FC
  closing_best <- df_t %>%
    filter(direction=="closing_vs_stable",
           motif %in% identity_motifs, log2fc > 0) %>%
    group_by(motif) %>%
    slice_max(log2fc, n=1, with_ties=FALSE) %>%
    ungroup() %>%
    mutate(best_ct = gsub(paste0("^",tissue_name,"_"), "", ct_key))

  row_group <- sapply(rownames(mat), function(m) {
    g <- classify_motif(m)
    if (g != "3_Cell Identity") return(g)
    ct <- closing_best$best_ct[closing_best$motif == m]
    if (length(ct) == 0 || is.na(ct[1])) return("3_zOther")
    paste0("3_", ct[1])
  })

  # Level order: AP-1 → CTCF → cell types in column order → Other
  id_levels   <- c(paste0("3_", ct_short), "3_zOther")
  id_present  <- id_levels[id_levels %in% unique(row_group)]
  all_levels  <- c("1_AP-1", "2_CTCF", id_present)
  row_class   <- factor(row_group, levels=all_levels)

  title_map  <- setNames(
    c("AP-1", "CTCF", gsub("^3_z?", "", id_present)),
    all_levels
  )
  n_id <- length(id_present)
  rt_cols <- c("#D6604D", "#4393C3",
               colorRampPalette(c("#333333","#888888"))(n_id))

  # ── Heatmap ───────────────────────────────────────────────────
  col_fn <- colorRamp2(c(-CAP_LOG2FC, 0, CAP_LOG2FC),
                       c("blue", "white", "red"))

  ht <- Heatmap(
    mat,
    name             = "log2(FC)\nvs stable",
    col              = col_fn,
    # Column labels colored like CM paper (red/grey/blue)
    column_labels    = col_labels,
    column_names_gp  = gpar(fontsize=7, col=col_label_colors, fontface="bold"),
    column_names_rot = 45,
    column_split     = col_split,
    column_title_gp  = gpar(fontsize=9, fontface="bold"),
    column_gap       = unit(2, "mm"),
    top_annotation   = ha_top,
    # Row split: AP-1 / CTCF / per-CT cell identity (fixed order)
    row_split              = row_class,
    row_title              = title_map[all_levels],
    row_title_gp           = gpar(fontsize=9, fontface="bold", col=rt_cols),
    cluster_row_slices     = FALSE,
    cluster_rows           = TRUE,
    cluster_columns        = FALSE,
    clustering_distance_rows = "pearson",
    clustering_method_rows   = "ward.D2",
    show_row_names   = TRUE,
    row_names_gp     = gpar(fontsize=5),
    border           = TRUE,
    use_raster       = TRUE,
    raster_quality   = 5,
    heatmap_legend_param = list(
      title_gp  = gpar(fontsize=8),
      labels_gp = gpar(fontsize=7),
      at        = c(-CAP_LOG2FC, -1, 0, 1, CAP_LOG2FC)
    )
  )

  n_row <- nrow(mat); n_col <- ncol(mat)
  h_in  <- max(8, n_row * 0.08 + 3)
  w_in  <- max(6, n_col * 0.28 + 3)

  title_str <- paste0(tissue_name, "  (", n_row, " motifs x ",
                      length(cts), " cell types)")
  for (ext in c("pdf","png")) {
    out_f <- file.path(OUT_DIR,
                       paste0(tissue_name, "_Fig4G_style.", ext))
    if (ext=="pdf") pdf(out_f, width=w_in, height=h_in)
    else            png(out_f, width=w_in*300, height=h_in*300, res=300)
    draw(ht, merge_legend=TRUE,
         column_title    = title_str,
         column_title_gp = gpar(fontsize=12, fontface="bold"))
    dev.off()
  }
  message("  Saved: ", tissue_name, "_Fig4G_style.pdf")
  message("  Matrix: ", n_row, " rows x ", n_col, " cols")
}

# ── Run ─────────────────────────────────────────────────────────
for (tissue in c("Aorta","Kidney","Lung","Tcell")) {
  df_t <- df %>% filter(tissue == !!tissue)
  if (nrow(df_t) == 0) next
  plot_heatmap(df_t, tissue)
}
message("\nDone. Output: ", OUT_DIR)
REOF

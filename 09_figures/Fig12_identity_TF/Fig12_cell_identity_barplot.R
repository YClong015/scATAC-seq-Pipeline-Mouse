library(ggplot2)
library(dplyr)

## HOMER sources: original = Lung/Aorta/Tcell, v5 = Kidney 11-type
HOMER_OUT_OLD <- "/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_closing_vs_opening/HOMER"
HOMER_OUT_V5  <- "/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_pseudobulk_Kidney_v5_DESeq2/HOMER_closing_vs_opening"

OUT_PLOT  <- "/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_closing_vs_opening/figures_filtered_v7_kidney11type"
dir.create(OUT_PLOT, showWarnings=FALSE, recursive=TRUE)

SHOW_N <- 10   # motifs shown per cell type
CAP_P  <- 50   # cap on -log10(P) for colour/size scales

## Read HOMER data
to.num      <- function(x) suppressWarnings(as.numeric(gsub(",|%","",as.character(x))))
# trimws() ensures leading/trailing whitespace in HOMER motif names is removed
short.motif <- function(x) trimws(sapply(strsplit(as.character(x),"/"),`[`,1))

# Source 1: original HOMER (Lung/Aorta/Tcell + old 5-type Kidney to drop)
known.files.old <- list.files(HOMER_OUT_OLD, pattern="^knownResults\\.txt$",
                              recursive=TRUE, full.names=TRUE)
# Source 2: v5 HOMER (11-type Kidney; we'll prepend "Kidney_" to cell type names)
known.files.v5  <- list.files(HOMER_OUT_V5,  pattern="^knownResults\\.txt$",
                              recursive=TRUE, full.names=TRUE)

# Tag each file with its source so we can name tissue_ct correctly
known.files <- c(
  setNames(known.files.old, rep("old", length(known.files.old))),
  setNames(known.files.v5,  rep("v5",  length(known.files.v5)))
)
message("Found ", length(known.files.old), " (old, original 4-tissue) + ",
        length(known.files.v5),  " (v5 Kidney 11-type) knownResults.txt files")

all.rows <- list()
for (i in seq_along(known.files)) {
  fp     <- known.files[[i]]
  source <- names(known.files)[i]
  dname  <- basename(dirname(fp))
  parts  <- strsplit(dname, "__", fixed=TRUE)[[1]]
  if (length(parts) < 2) next

  if (source == "v5") {
    # v5 folder names: <CellType>__Day42_vs_Sham__closing_vs_opening
    # -> prepend Kidney_ to match convention used elsewhere
    tissue_ct <- paste0("Kidney_", parts[1])
    # Drop Kidney_EC to stay consistent with the 8-CT order used in
    # Opening_dotplot.R and 4tissues_integrate_heatmap.R.
    if (tissue_ct == "Kidney_EC") next
  } else {
    # original folder names: <Tissue>_<CellType>__closing_vs_opening
    tissue_ct <- parts[1]
    # DROP old 5-type Kidney entries (v5 supersedes them)
    if (startsWith(tissue_ct, "Kidney_")) next
  }

  # Collapse self-redundant names, e.g. "Tcell_Tcell" -> "Tcell".
  tc.parts <- strsplit(tissue_ct, "_", fixed = TRUE)[[1]]
  if (length(tc.parts) == 2 && tc.parts[1] == tc.parts[2]) {
    tissue_ct <- tc.parts[1]
  }

  # The "Tcell" dataset is a CD8 T cell pool; display it as "T cells".
  if (tissue_ct == "Tcell") {
    tissue_ct <- "T cells"
  }

  tab <- tryCatch(
    read.delim(fp, header=TRUE, sep="\t",
               stringsAsFactors=FALSE, check.names=FALSE),
    error=function(e) NULL
  )
  if (is.null(tab) || nrow(tab)==0) next
  if (!all(c("Motif Name","P-value") %in% colnames(tab))) next

  tab$pval   <- to.num(tab[["P-value"]])
  tab$log10p <- -log10(pmax(tab$pval, 1e-300))
  tab$motif  <- short.motif(tab[["Motif Name"]])
  tab$tissue_ct <- tissue_ct

  pct.t.col <- grep("% of Target",     colnames(tab), value=TRUE)[1]
  pct.b.col <- grep("% of Background", colnames(tab), value=TRUE)[1]
  if (!is.na(pct.t.col) && !is.na(pct.b.col)) {
    pct.t <- to.num(tab[[pct.t.col]])
    pct.b <- to.num(tab[[pct.b.col]])
    tab$fold_enrich <- ifelse(pct.b > 0, pct.t / pct.b, pct.t / 0.01)
  } else {
    tab$fold_enrich <- NA_real_
  }

  all.rows[[length(all.rows)+1]] <-
    tab[, c("motif","log10p","pval","fold_enrich","tissue_ct")]
}

if (length(all.rows)==0) stop("No HOMER results found.")
df <- dplyr::bind_rows(all.rows)
message("Loaded: ", nrow(df), " rows, ",
        length(unique(df$tissue_ct)), " cell types")

## Keep the most significant occurrence of each motif per cell type, then keep
## significant motifs. No blacklist / whitelist / fold-enrichment curation.
dat.clean <- df %>%
  group_by(tissue_ct, motif) %>%
  slice_min(order_by = pval, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  filter(pval < 0.05)

message("After dedup + significance (pval < 0.05): ", nrow(dat.clean), " rows across ",
        length(unique(dat.clean$tissue_ct)), " cell types")

## Barplot per cell type
for (ct in sort(unique(dat.clean$tissue_ct))) {
  dat.ct <- dat.clean %>%
    filter(tissue_ct == ct) %>%
    slice_max(order_by=log10p, n=SHOW_N, with_ties=FALSE) %>%
    mutate(log10p_plot = pmin(log10p, CAP_P))

  if (nrow(dat.ct) == 0) { message("No motifs: ", ct); next }

  motif.ord <- dat.ct %>% arrange(log10p_plot) %>% pull(motif)
  dat.ct$motif <- factor(dat.ct$motif, levels=motif.ord)

  p <- ggplot(dat.ct, aes(x=log10p_plot, y=motif)) +
    geom_col(fill="#2166AC", width=0.72) +
    geom_vline(xintercept=-log10(0.05), linetype="dashed",
               color="grey25", linewidth=0.8) +
    labs(title=paste0(ct, ", cell identity TFs"),
         x = expression(-log[10] * "(" * italic(P) * ")"),
         y = "TF Motif") +
    theme_classic(base_size=20, base_family="sans") +
    theme(
      plot.title        = element_text(face="bold", size=26, hjust=0,
                                       margin=margin(b=10),
                                       family="sans", colour="black"),
      plot.subtitle     = element_text(hjust=0, size=14, colour="grey25",
                                       margin=margin(b=12), family="sans"),
      axis.text.x       = element_text(size=18, colour="black", family="sans"),
      axis.text.y       = element_text(size=16, colour="black", family="sans"),
      axis.title.x      = element_text(size=20, margin=margin(t=10),
                                       family="sans", colour="black"),
      axis.title.y      = element_text(size=20, margin=margin(r=10),
                                       family="sans", colour="black"),
      panel.border      = element_rect(colour="black", fill=NA, linewidth=0.8),
      axis.line         = element_blank(),
      axis.ticks        = element_line(colour="black", linewidth=0.8),
      axis.ticks.length = unit(5, "pt"),
      panel.grid.major.x = element_line(colour="grey90", linewidth=0.4),
      panel.grid.major.y = element_blank(),
      panel.grid.minor  = element_blank(),
      plot.margin       = margin(16, 22, 14, 18)
    )

  safe.ct <- gsub("[^A-Za-z0-9_]","_", ct)
  h <- max(8, nrow(dat.ct)*0.55 + 3)
  ggsave(file.path(OUT_PLOT, paste0(safe.ct,"_identity_v6_barplot.pdf")),
         p, width=12, height=h, limitsize=FALSE)
  ggsave(file.path(OUT_PLOT, paste0(safe.ct,"_identity_v6_barplot.png")),
         p, width=12, height=h, dpi=300, limitsize=FALSE)
  message("Saved: ", ct)
}

## Dotplot
DOT_N  <- 10
top.df <- dat.clean %>%
  group_by(tissue_ct) %>%
  slice_max(order_by=log10p, n=DOT_N, with_ties=FALSE) %>%
  ungroup() %>%
  mutate(log10p_plot = pmin(log10p, CAP_P))

if (nrow(top.df) > 0) {
  motif.ord <- top.df %>%
    group_by(motif) %>%
    summarise(m=max(log10p_plot), .groups="drop") %>%
    arrange(desc(m)) %>% pull(motif)
  top.df$motif     <- factor(top.df$motif, levels=rev(motif.ord))
  top.df$tissue_ct <- factor(top.df$tissue_ct,
                             levels=sort(unique(top.df$tissue_ct)))

  n.m <- length(unique(top.df$motif))
  n.c <- length(unique(top.df$tissue_ct))

  p.dot <- ggplot(top.df, aes(x=tissue_ct, y=motif,
                              size=log10p_plot, color=log10p_plot)) +
    geom_point() +
    scale_size_continuous(range=c(2,12), name=expression(-log[10]*"("*italic(P)*")"),
                          limits=c(0,CAP_P)) +
    scale_color_gradient(low="grey88", high="#2166AC",
                         name=expression(-log[10]*"("*italic(P)*")"),
                         limits=c(0,CAP_P)) +
    labs(
      title = "Cell identity TFs, closing vs opening",
      x = NULL, y = "TF Motif"
    ) +
    theme_classic(base_size=20, base_family="sans") +
    theme(
      plot.title         = element_text(face="bold", size=26, hjust=0,
                                        margin=margin(b=10),
                                        family="sans", colour="black"),
      plot.subtitle      = element_text(hjust=0, size=14, colour="grey25",
                                        margin=margin(b=12), family="sans"),
      axis.text.x        = element_text(angle=45, hjust=1, size=18,
                                        colour="black", family="sans"),
      axis.text.y        = element_text(size=16, colour="black", family="sans"),
      axis.title.y       = element_text(size=20, margin=margin(r=10),
                                        family="sans", colour="black"),
      panel.border       = element_rect(colour="black", fill=NA, linewidth=0.8),
      axis.line          = element_blank(),
      axis.ticks         = element_line(colour="black", linewidth=0.8),
      axis.ticks.length  = unit(5, "pt"),
      panel.grid.major   = element_line(colour="grey92", linewidth=0.4),
      panel.grid.minor   = element_blank(),
      legend.position    = "right",
      legend.text        = element_text(size=16, family="sans", colour="black"),
      legend.title       = element_text(size=18, family="sans", colour="black"),
      legend.key.size    = unit(0.8, "cm"),
      plot.margin        = margin(16, 22, 14, 18)
    )

  h <- min(40, max(10, n.m * 0.50 + 3))
  w <- min(36, max(12, n.c * 1.8 + 5))

  pdf.path <- file.path(OUT_PLOT, "CellIdentity_v6_dotplot.pdf")
  png.path <- file.path(OUT_PLOT, "CellIdentity_v6_dotplot.png")
  ggsave(pdf.path, p.dot, width = w, height = h, limitsize = FALSE)
  ggsave(png.path, p.dot, width = w, height = h, dpi = 300, limitsize = FALSE)

  message(sprintf("Saved combined dotplot (%.1f x %.1f in, %d motifs x %d CTs)",
                  w, h, n.m, n.c))
  message("  -> ", pdf.path)
  message("  -> ", png.path)

  if (interactive()) print(p.dot)
} else {
  message("Combined dotplot SKIPPED: top.df has 0 rows (no significant motifs).")
}

message("\nDone. Output: ", OUT_PLOT)

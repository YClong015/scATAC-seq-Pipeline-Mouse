library(ggplot2)
library(dplyr)

# ── Paths ───────────────────────────────────────────────────────
HOMER_OUT <- "/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_closing_vs_opening/HOMER"
OUT_PLOT  <- "/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_closing_vs_opening/figures_filtered_v6"
dir.create(OUT_PLOT, showWarnings=FALSE, recursive=TRUE)

SHOW_N    <- 20
MIN_SHOW  <- 30   # fallback target: try to reach this many motifs per CT
CAP_P     <- 50

# ── Fold enrichment threshold ───────────────────────────────────
MIN_FOLD <- 1.3   # uniform across all cell types

# ── Blacklist ───────────────────────────────────────────────────
BLACKLIST <- c(
  # General transcription machinery
  "TATA-Box", "^TBP", "^NFY", "^GFY", "^GFY-Staf",
  # Cell cycle
  "^E2F", "TFDP",
  # Ubiquitous / circadian bHLH
  "bHLHE", "USF", "Usf", "CLOCK", "BMAL",
  "^MAX\\(", "^Max\\(", "^MRE\\(",
  # Lysosomal / MiT-TFE bHLH (ubiquitous, not cell identity)
  "^TFE3\\(", "^TFEB\\(", "^MITF\\(", "^Mitf\\(",
  # Ubiquitous zinc fingers
  "^YY1", "^Maz\\(", "^SP[123]\\(", "^Sp[123]\\(",
  "^NRF1\\(", "^NRF2\\(", "^NRF\\(NRF\\)", "Ronin", "^CRE\\(bZIP\\)$",
  # ZNF proteins that appear as noise across unrelated cell types
  "^ZNF711\\(", "^Znf263\\(", "^ZNF519\\(", "^HINFP\\(",
  "^ZNF264\\(", "^BANP\\(", "^ZFX\\(",
  "^ZNF341\\(", "^ZNF415\\(", "^ZNF189\\(", "^ZNF692\\(",
  "^Zfp809\\(", "^ZFP3\\(", "^ZNF528\\(",
  # Chromatin architecture
  "^CTCF\\(", "^BORIS\\(",
  # Stem cell / reprogramming
  "^Nanog\\(", "^Esrrb\\(", "^ESRRB\\(",
  "OCT4-SOX2-TCF-NANOG", "^Zscan4",
  # Hormone response elements (ARE/GRE = motif sequences, not cell identity)
  # Covers both plain (ARE(NR)) and composite (ARE(NR),IR3 / GRE(NR),IR3) forms
  "^ARE\\(", "^GRE\\(", "^ARE$", "^GRE$",
  # Steroid receptors inappropriate in non-reproductive somatic tissues
  "^PGR\\(", "^AR\\(",
  # Erythroid-specific KLF (EKLF = KLF1, not expressed in non-erythroid cells)
  "^EKLF\\(",
  # Early B-cell factors appearing in non-B/non-immune cells
  "^EBF[12]?\\(", "^EBF\\(EBF\\)",
  # Ubiquitous ETS family (MAP kinase targets, broadly expressed)
  # Keep: ERG, FLI1 (endothelial), SPI1/PU.1 (macrophage), ETV2 (endothelial)
  "^GABPA\\(", "^Gabpa\\(",
  "^ELF1\\(", "^Elf1\\(",
  "^ELF5\\(", "^Elf5\\(",
  "^Elk1\\(", "^ELK1\\(",
  "^Elk4\\(", "^ELK4\\(",
  "^ETV1\\(", "^Etv1\\(",
  "^ETV4\\(", "^Etv4\\(",
  # Brain / neural TFs
  "^Lhx6\\b", "^NeuroD1\\(", "^Atoh1\\(", "^Atoh7\\(",
  "^Ascl[12]\\(", "^n-Myc\\(",
  # Skeletal muscle
  "^MyoD\\(", "^Myf5\\(", "^MYOD\\(",
  # HOX (inappropriate in adult somatic cells; covers 1- and 2-digit numbering,
  #       both lowercase mouse and uppercase human/HOMER naming)
  "^Hox[a-d][0-9]{1,2}\\(", "^HOX[A-D][0-9]{1,2}\\(",
  # Early embryonic / gonadal
  "^GSC\\(", "^Pax7\\(",
  # Cardiac / neural crest TFs appearing as noise
  "^Hand[12]\\(", "^HAND[12]\\(",
  "^GLI[123]\\(", "^Gli[123]\\(",
  "^Pit1", "^Otx[12]\\("
)

is_blacklisted <- function(x) grepl(paste(BLACKLIST, collapse="|"), x, perl=TRUE)

# ── Whitelists (used by CT_FILTER below) ────────────────────────

AT2_WHITELIST <- c(
  "NKX2", "TTF",          # NKX2-1/TTF1 — master AT2 TF
  "FOXA",                 # FOXA1/FOXA2 — pioneer TF for lung epithelium
  "SOX[29]",              # SOX2/SOX9 — lung progenitor
  "GATA[46]",             # GATA4/GATA6 — lung/cardiac lineage
  "CEBP", "C/EBP",        # C/EBPa/b — alveolar differentiation
  "KLF[245]",             # KLF4/5 — epithelial
  "SPDEF",                # SPDEF — goblet/club cell, lung epithelial
  "SMAD[234]?",           # SMAD2/3 — TGF-b, relevant to COPD AT2
  "HIF-1a", "HIF1A",      # hypoxia response, relevant to COPD
  "NFATC[124]?"           # NFAT — inflammatory response in AT2
)

MAC_WHITELIST <- c(
  # Core myeloid lineage
  "SPI1", "PU\\.1", "Spi1",         # PU.1/SPI1 — master macrophage TF
  "IRF[1248]",                       # IRF1/2/4/8 — macrophage identity & polarisation
  "CEBP", "C/EBP", "Cebp",          # C/EBP family — myeloid differentiation
  "RUNX[13]", "Runx[13]",            # RUNX1/3 — myeloid commitment
  # AP-1 superfamily
  "^Jun", "^JUN",                    # JUN/JUNB/JUND
  "^Fos", "^FOS",                    # FOS/FOSB
  "^Fra[12]\\(", "^FRA[12]\\(",      # Fra1/Fra2 (FOSL1/2)
  "^Fosl[12]\\(", "^FOSL[12]\\(",
  "^Atf[13]\\(", "^ATF[13]\\(",     # ATF1/3
  "BATF",                            # BATF
  # KLF
  "KLF[24]", "Klf[24]",
  # Nuclear receptors
  "Nur77", "NR4A", "Nr4a",
  "PPAR[gGdD]", "Pparg", "PPARG",   # PPARg — M2, foam cells
  "LXR", "NR1H[23]",                # LXR — cholesterol/lipid
  "RXR[AB]?", "Rxr[ab]?",
  # MAF family
  "MAFB", "MafB",
  "^c-Maf\\(", "^MAF\\(",
  # NF-kB / STAT / BACH / NFAT
  "RELA", "RelA", "NF-kB", "NFKB",
  "STAT[136]", "Stat[136]",
  "BACH[12]", "Bach[12]",
  "NFATC[1-4]?", "Nfatc[1-4]?",
  # Egr / ZEB / TEAD
  "Egr[12]\\(", "EGR[12]\\(",
  "ZEB[12]\\(", "Zeb[12]\\(",
  "TEAD[1234]?\\("
)

# B-cell-specific TFs exempted from the global blacklist
# (EBF1/2 are master B cell TFs; PAX5, BCL6, BACH2, OCT2 are B cell identity)
BCELL_EXEMPTIONS <- c(
  "^EBF",                            # EBF1/EBF2 — master B cell TF
  "PAX5", "Pax5",                    # B cell commitment
  "BCL6", "Bcl6",                    # germinal centre
  "IRF4", "IRF8",                    # plasma cell / B cell activation
  "BACH2", "Bach2",                  # B cell maturation
  "^OCT2\\(", "^POU2F2\\(",          # OCT2 — B cell
  "^Bob1\\(", "^OBF1\\("            # OCT co-activators
)

# ── Per-cell-type filter configuration ──────────────────────────
# mode options:
#   "default"     global blacklist + fold enrichment  (most cell types)
#   "whitelist"   only motifs matching `patterns` pass; no blacklist used
#   "passthrough" p-value filter only; no blacklist, no fold enrichment
#   "custom"      global blacklist + fold, but motifs matching `exemptions`
#                 are allowed through even if blacklisted
#
# Cell types not listed here automatically use "default".

CT_FILTER <- list(
  # ── T cells: p-value only (low peak counts; no spurious noise) ─
  "Tcell_Tcell" = list(mode = "passthrough"),
  "Lung_T"      = list(mode = "passthrough"),
  
  # ── B cells: standard filters + exempt B-cell-specific TFs ────
  "Lung_B" = list(mode = "custom", exemptions = BCELL_EXEMPTIONS),
  
  # ── AT2: whitelist (ubiquitous ETS noise dominates otherwise) ──
  "Lung_AT2" = list(mode = "whitelist", patterns = AT2_WHITELIST),
  
  # ── Macrophages: whitelist (epithelial/stromal NR noise) ───────
  "Kidney_Macrophages" = list(mode = "whitelist", patterns = MAC_WHITELIST),
  "Aorta_Macrophages"  = list(mode = "whitelist", patterns = MAC_WHITELIST),
  "Lung_Mac-alv"       = list(mode = "whitelist", patterns = MAC_WHITELIST),
  "Lung_Mac_alv"       = list(mode = "whitelist", patterns = MAC_WHITELIST)
)

# Vectorised filter: returns logical vector, one value per row
passes_ct_filter <- function(motif_vec, ct_vec, fold_vec) {
  mapply(function(m, ct, fe) {
    cfg     <- if (!is.null(CT_FILTER[[ct]])) CT_FILTER[[ct]] else list(mode = "default")
    fold_ok <- is.na(fe) | fe >= MIN_FOLD
    switch(cfg$mode,
           "passthrough" = TRUE,
           "whitelist"   = grepl(paste(cfg$patterns,    collapse="|"), m,
                                 perl=TRUE, ignore.case=TRUE),
           "default"     = !is_blacklisted(m) & fold_ok,
           "custom"      = {
             exempt <- grepl(paste(cfg$exemptions, collapse="|"), m, perl=TRUE)
             (!is_blacklisted(m) | exempt) & fold_ok
           }
    )
  }, motif_vec, ct_vec, fold_vec, SIMPLIFY=TRUE, USE.NAMES=FALSE)
}

# ── Read HOMER data ──────────────────────────────────────────────
to_num      <- function(x) suppressWarnings(as.numeric(gsub(",|%","",as.character(x))))
# trimws() ensures leading/trailing whitespace in HOMER motif names does not
# break the ^-anchored blacklist patterns
short_motif <- function(x) trimws(sapply(strsplit(as.character(x),"/"),`[`,1))

known_files <- list.files(HOMER_OUT, pattern="^knownResults\\.txt$",
                          recursive=TRUE, full.names=TRUE)
message("Found ", length(known_files), " knownResults.txt files")

all_rows <- list()
for (fp in known_files) {
  dname <- basename(dirname(fp))
  parts <- strsplit(dname, "__", fixed=TRUE)[[1]]
  if (length(parts) < 2) next
  tissue_ct <- parts[1]
  
  tab <- tryCatch(
    read.delim(fp, header=TRUE, sep="\t",
               stringsAsFactors=FALSE, check.names=FALSE),
    error=function(e) NULL
  )
  if (is.null(tab) || nrow(tab)==0) next
  if (!all(c("Motif Name","P-value") %in% colnames(tab))) next
  
  tab$pval   <- to_num(tab[["P-value"]])
  tab$log10p <- -log10(pmax(tab$pval, 1e-300))
  tab$motif  <- short_motif(tab[["Motif Name"]])
  tab$tissue_ct <- tissue_ct
  
  pct_t_col <- grep("% of Target",     colnames(tab), value=TRUE)[1]
  pct_b_col <- grep("% of Background", colnames(tab), value=TRUE)[1]
  if (!is.na(pct_t_col) && !is.na(pct_b_col)) {
    pct_t <- to_num(tab[[pct_t_col]])
    pct_b <- to_num(tab[[pct_b_col]])
    tab$fold_enrich <- ifelse(pct_b > 0, pct_t / pct_b, pct_t / 0.01)
  } else {
    tab$fold_enrich <- NA_real_
  }
  
  all_rows[[length(all_rows)+1]] <-
    tab[, c("motif","log10p","pval","fold_enrich","tissue_ct")]
}

if (length(all_rows)==0) stop("No HOMER results found.")
df <- dplyr::bind_rows(all_rows)
message("Loaded: ", nrow(df), " rows, ",
        length(unique(df$tissue_ct)), " cell types")

# ── Per-CT fallback filter ───────────────────────────────────────
# Level 1 (strict)       : pval < 0.05 + full CT_FILTER (blacklist/whitelist + fold)
# Level 2 (relaxed fold) : pval < 0.05 + CT_FILTER with fold check bypassed
# Level 3 (wl→default)   : pval < 0.05 + blacklist-only (for whitelist-mode CTs only)
# Level 4 (relaxed p)    : pval < 0.10 + blacklist-only, no fold
# Each CT uses the first level that yields >= MIN_SHOW motifs.

filter_with_fallback <- function(df_ct, ct_name) {
  base <- df_ct %>%
    group_by(motif) %>% arrange(pval) %>% slice(1) %>% ungroup()
  
  # Level 1: strict
  r <- base %>%
    filter(pval < 0.05) %>%
    mutate(pass = passes_ct_filter(motif, ct_name, fold_enrich)) %>%
    filter(pass) %>% select(-pass) %>%
    mutate(filter_level = "strict")
  if (nrow(r) >= MIN_SHOW) return(r)
  
  # Level 2: bypass fold enrichment check (pass NA so fold_ok = TRUE)
  r <- base %>%
    filter(pval < 0.05) %>%
    mutate(pass = passes_ct_filter(motif, ct_name, NA_real_)) %>%
    filter(pass) %>% select(-pass) %>%
    mutate(filter_level = "relaxed_fold")
  if (nrow(r) >= MIN_SHOW) return(r)
  
  # Level 3: for whitelist-mode CTs, fall back to blacklist-only
  cfg <- CT_FILTER[[ct_name]]
  if (!is.null(cfg) && cfg$mode == "whitelist") {
    r <- base %>%
      filter(pval < 0.05) %>%
      filter(!is_blacklisted(motif)) %>%
      mutate(filter_level = "whitelist_fallback")
    if (nrow(r) >= MIN_SHOW) return(r)
  }
  
  # Level 4: relax p-value to 0.10, blacklist-only (custom-mode CTs keep their exemptions)
  cfg <- CT_FILTER[[ct_name]]
  if (!is.null(cfg) && cfg$mode == "custom") {
    r <- base %>%
      filter(pval < 0.10) %>%
      mutate(exempt = grepl(paste(cfg$exemptions, collapse="|"), motif, perl=TRUE)) %>%
      filter(!is_blacklisted(motif) | exempt) %>%
      select(-exempt) %>%
      mutate(filter_level = "relaxed_p")
  } else {
    r <- base %>%
      filter(pval < 0.10) %>%
      filter(!is_blacklisted(motif)) %>%
      mutate(filter_level = "relaxed_p")
  }
  r
}

# Apply per-CT fallback and combine
dat_clean <- dplyr::bind_rows(
  lapply(split(df, df$tissue_ct), function(df_ct) {
    ct <- unique(df_ct$tissue_ct)
    filter_with_fallback(df_ct, ct)
  })
)

message("After all filters: ", nrow(dat_clean), " rows across ",
        length(unique(dat_clean$tissue_ct)), " cell types")
dat_clean %>%
  count(tissue_ct, filter_level) %>%
  arrange(n) %>%
  { message("Motifs per cell type:\n",
            paste(sprintf("  %-30s %3d  [%s]", .$tissue_ct, .$n, .$filter_level),
                  collapse="\n")) }

# ── Barplot per cell type ────────────────────────────────────────
for (ct in sort(unique(dat_clean$tissue_ct))) {
  dat_ct <- dat_clean %>%
    filter(tissue_ct == ct) %>%
    slice_max(order_by=log10p, n=SHOW_N, with_ties=FALSE) %>%
    mutate(log10p_plot = pmin(log10p, CAP_P))

  if (nrow(dat_ct) == 0) { message("No motifs: ", ct); next }

  flevel <- unique(dat_ct$filter_level)[1]
  filter_note <- switch(flevel,
                        "strict"            = paste0("Fold enrichment ≥", MIN_FOLD, "× | Generic TFs removed"),
                        "relaxed_fold"      = "Fold filter relaxed | Generic TFs removed",
                        "whitelist_fallback"= paste0("Fold enrichment ≥", MIN_FOLD, "× | Generic TFs removed"),
                        "relaxed_p"         = "p-value threshold relaxed to 0.10 | Generic TFs removed",
                        flevel
  )
  subtitle_str <- paste0(
    "Foreground: closing peaks | Background: opening peaks\n",
    filter_note
  )

  motif_ord <- dat_ct %>% arrange(log10p_plot) %>% pull(motif)
  dat_ct$motif <- factor(dat_ct$motif, levels=motif_ord)

  # Separator line between rank-10 and rank-11 from top
  sep_y <- nrow(dat_ct) - 10 + 0.5

  p <- ggplot(dat_ct, aes(x=log10p_plot, y=motif)) +
    geom_col(fill="#2166AC") +
    geom_hline(yintercept=sep_y, color="#D73027",
               linetype="dashed", linewidth=0.5) +
    annotate("text", x=Inf, y=sep_y + 0.4, label="Top 10",
             hjust=1.1, size=2.8, color="#D73027", fontface="bold") +
    geom_vline(xintercept=-log10(0.05), linetype="dashed",
               color="grey50", linewidth=0.4) +
    labs(title=paste0(ct, " — cell identity TFs"),
         subtitle=subtitle_str, x="-log10(P)", y="TF Motif") +
    theme_bw(base_size=11) +
    theme(
      axis.text.y   = element_text(size=9),
      plot.title    = element_text(face="bold", hjust=0.5),
      plot.subtitle = element_text(hjust=0.5, size=9, color="grey40")
    )
  
  safe_ct <- gsub("[^A-Za-z0-9_]","_", ct)
  h <- max(5, nrow(dat_ct)*0.32 + 2)
  ggsave(file.path(OUT_PLOT, paste0(safe_ct,"_identity_v6_barplot.pdf")),
         p, width=8, height=h, limitsize=FALSE)
  ggsave(file.path(OUT_PLOT, paste0(safe_ct,"_identity_v6_barplot.png")),
         p, width=8, height=h, dpi=300, limitsize=FALSE)
  message("Saved: ", ct)
}

# ── Dotplot ──────────────────────────────────────────────────────
DOT_N  <- 10
top_df <- dat_clean %>%
  group_by(tissue_ct) %>%
  slice_max(order_by=log10p, n=DOT_N, with_ties=FALSE) %>%
  ungroup() %>%
  mutate(log10p_plot = pmin(log10p, CAP_P))

if (nrow(top_df) > 0) {
  motif_ord <- top_df %>%
    group_by(motif) %>%
    summarise(m=max(log10p_plot), .groups="drop") %>%
    arrange(desc(m)) %>% pull(motif)
  top_df$motif     <- factor(top_df$motif, levels=rev(motif_ord))
  top_df$tissue_ct <- factor(top_df$tissue_ct,
                             levels=sort(unique(top_df$tissue_ct)))
  
  n_m <- length(unique(top_df$motif))
  n_c <- length(unique(top_df$tissue_ct))
  
  p_dot <- ggplot(top_df, aes(x=tissue_ct, y=motif,
                              size=log10p_plot, color=log10p_plot)) +
    geom_point() +
    scale_size_continuous(range=c(1,8), name="-log10(P)", limits=c(0,CAP_P)) +
    scale_color_gradient(low="grey88", high="#2166AC",
                         name="-log10(P)", limits=c(0,CAP_P)) +
    labs(
      title    = "Cell identity TFs — closing vs opening",
      subtitle = paste0("Blacklist (generic + context-inappropriate TFs removed)",
                        " | Fold enrichment ≥", MIN_FOLD, "×"),
      x=NULL, y="TF Motif"
    ) +
    theme_bw(base_size=12) +
    theme(
      axis.text.x      = element_text(angle=45, hjust=1, size=10),
      axis.text.y      = element_text(size=8),
      plot.title       = element_text(face="bold", hjust=0.5),
      plot.subtitle    = element_text(hjust=0.5, size=9, color="grey40"),
      panel.grid.major = element_line(color="grey92")
    )
  
  h <- max(6, n_m*0.28 + 2)
  w <- max(8, n_c*1.3 + 3)
  ggsave(file.path(OUT_PLOT, "CellIdentity_v6_dotplot.pdf"),
         p_dot, width=w, height=h, limitsize=FALSE)
  ggsave(file.path(OUT_PLOT, "CellIdentity_v6_dotplot.png"),
         p_dot, width=w, height=h, dpi=300, limitsize=FALSE)
  message("Saved dotplot")
}

message("\nDone. Output: ", OUT_PLOT)


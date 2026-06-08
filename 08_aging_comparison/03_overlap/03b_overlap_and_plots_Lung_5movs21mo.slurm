#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --job-name=overlap_Lung_5movs21mo
#SBATCH --time=04:00:00
#SBATCH --partition=general
#SBATCH --account=a_imb_ccbcd
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

set -euo pipefail

module load bedtools/2.30.0-gcc-11.3.0
module load r/4.4.2

TISSUE="Lung"
BASE="/QRISdata/Q8448/Mouse_disease_data/DAR"
AGING_DIR="${BASE}/DAR_science_comparison/aging_DARs_5movs21mo"
TMP_DIR="${BASE}/DAR_science_comparison/tmp_beds_Lung_5movs21mo"
OUT_DIR="${BASE}/DAR_science_comparison/results_Lung_5movs21mo"
mkdir -p "${TMP_DIR}" "${OUT_DIR}/figures"

CONTRAST="Case_vs_Control"
DAR_TABLES="${BASE}/DAR_pseudobulk_DESeq2/DAR_pseudobulk_Lung_DESeq2/DAR_tables"

declare -A CT_MAP=()

PADJ_CUT=0.05

### Step 1: Convert peaks to BED files ###
echo "$(date)  Step 1: Generating BED files..."

for aging_tsv in "${AGING_DIR}"/*_Aged_vs_Adult_DAR.tsv; do
    base=$(basename "${aging_tsv}" _Aged_vs_Adult_DAR.tsv)
    tissue="${base%%_*}"
    [[ "${tissue}" != "${TISSUE}" ]] && continue
    ct="${base#*_}"

    # -- Aging BED (Science paper peaks: chr1:start-end) ------
    # Columns: baseMean(1) log2FC(2) lfcSE(3) stat(4) pvalue(5) padj(6) peak(7) cell_type(8) tissue(9)
    awk -F'\t' -v cut="${PADJ_CUT}" '
      NR>1 && $6~/^[0-9]/ && $6+0 < cut && $2+0 > 0 {
        gsub(/:/, "\t", $7); gsub(/-/, "\t", $7); print $7
      }' "${aging_tsv}" | sort -k1,1 -k2,2n \
      > "${TMP_DIR}/${base}_aging_open.bed"

    awk -F'\t' -v cut="${PADJ_CUT}" '
      NR>1 && $6~/^[0-9]/ && $6+0 < cut && $2+0 < 0 {
        gsub(/:/, "\t", $7); gsub(/-/, "\t", $7); print $7
      }' "${aging_tsv}" | sort -k1,1 -k2,2n \
      > "${TMP_DIR}/${base}_aging_close.bed"

    n_ao=$(wc -l < "${TMP_DIR}/${base}_aging_open.bed")
    n_ac=$(wc -l < "${TMP_DIR}/${base}_aging_close.bed")
    echo "  ${base}: aging_open=${n_ao}  aging_close=${n_ac}"

    # -- Relaxed aging BED (padj<0.2) --------------------------
    awk -F'\t' '
      NR>1 && $6~/^[0-9]/ && $6+0 < 0.2 && $2+0 > 0 {
        gsub(/:/, "\t", $7); gsub(/-/, "\t", $7); print $7
      }' "${aging_tsv}" | sort -k1,1 -k2,2n \
      > "${TMP_DIR}/${base}_aging_open_relaxed.bed"

    awk -F'\t' '
      NR>1 && $6~/^[0-9]/ && $6+0 < 0.2 && $2+0 < 0 {
        gsub(/:/, "\t", $7); gsub(/-/, "\t", $7); print $7
      }' "${aging_tsv}" | sort -k1,1 -k2,2n \
      > "${TMP_DIR}/${base}_aging_close_relaxed.bed"

    n_aor=$(wc -l < "${TMP_DIR}/${base}_aging_open_relaxed.bed")
    n_acr=$(wc -l < "${TMP_DIR}/${base}_aging_close_relaxed.bed")
    echo "  ${base}: aging_open_relaxed=${n_aor}  aging_close_relaxed=${n_acr}"

    # -- Disease BED (our peaks: chr1-start-end) --------------
    # Columns: peak(1=chr-start-end) baseMean(2) log2FC(3) ... padj(7)
    # Filter: padj<0.05 only (no LFC threshold, matching Lu 2026 / Patrick 2024 aging atlases)
    mapped_cts="${CT_MAP[${base}]:-${ct}}"
    tmp_combined="${TMP_DIR}/${base}_disease_combined.tsv"
    first_file=true
    for mapped_ct in ${mapped_cts}; do
        f="${DAR_TABLES}/${mapped_ct}__${CONTRAST}__005_DESeq2_all.tsv"
        if [ -f "${f}" ]; then
            if ${first_file}; then
                cat "${f}" > "${tmp_combined}"; first_file=false
            else
                tail -n +2 "${f}" >> "${tmp_combined}"
            fi
        fi
    done
    if ${first_file}; then
        echo "  SKIP no disease TSVs found for ${base} (mapped: ${mapped_cts})"; continue
    fi
    dis_tsv="${tmp_combined}"

    awk -F'\t' -v cut="${PADJ_CUT}" '
      NR>1 && $7~/^[0-9]/ && $7+0 < cut && $3+0 > 0 {
        split($1, a, "-"); print a[1]"\t"a[2]"\t"a[3]
      }' "${dis_tsv}" | sort -k1,1 -k2,2n \
      > "${TMP_DIR}/${base}_disease_open.bed"

    awk -F'\t' -v cut="${PADJ_CUT}" '
      NR>1 && $7~/^[0-9]/ && $7+0 < cut && $3+0 < 0 {
        split($1, a, "-"); print a[1]"\t"a[2]"\t"a[3]
      }' "${dis_tsv}" | sort -k1,1 -k2,2n \
      > "${TMP_DIR}/${base}_disease_close.bed"

    n_do=$(wc -l < "${TMP_DIR}/${base}_disease_open.bed")
    n_dc=$(wc -l < "${TMP_DIR}/${base}_disease_close.bed")
    echo "  ${base}: disease_open=${n_do}  disease_close=${n_dc}"
done

### Step 2: bedtools intersect ###
echo ""
echo "$(date)  Step 2: bedtools intersect..."

PYTHON="/home/s4869245/.conda/envs/scanpy_env/bin/python"

${PYTHON} - <<'PYEOF'
import os, subprocess
import pandas as pd

TMP_DIR   = "/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_science_comparison/tmp_beds_Lung_5movs21mo"
OUT_DIR   = "/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_science_comparison/results_Lung_5movs21mo"
AGING_DIR = "/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_science_comparison/aging_DARs_5movs21mo"

def count_bed(f):
    if not os.path.exists(f): return 0
    with open(f) as fh:
        return sum(1 for l in fh if l.strip())

def bedtools_overlap(a, b):
    if not os.path.exists(a) or not os.path.exists(b): return 0
    if count_bed(a) == 0 or count_bed(b) == 0: return 0
    res = subprocess.run(
        ["bedtools", "intersect", "-a", a, "-b", b, "-u"],
        capture_output=True, text=True, check=True
    )
    return len([l for l in res.stdout.strip().split('\n') if l.strip()])

rows = []
for f in sorted(os.listdir(TMP_DIR)):
    if not f.endswith("_aging_open.bed"): continue
    base = f.replace("_aging_open.bed", "")
    tissue = base.split("_")[0]
    ct     = base[len(tissue)+1:]

    aging_tsv = f"{AGING_DIR}/{base}_Aged_vs_Adult_DAR.tsv"
    if os.path.exists(aging_tsv):
        with open(aging_tsv) as fh:
            n_tested = sum(1 for line in fh if line.strip()) - 1
    else:
        n_tested = None

    ao  = f"{TMP_DIR}/{base}_aging_open.bed"
    ac  = f"{TMP_DIR}/{base}_aging_close.bed"
    aor = f"{TMP_DIR}/{base}_aging_open_relaxed.bed"
    acr = f"{TMP_DIR}/{base}_aging_close_relaxed.bed"
    do  = f"{TMP_DIR}/{base}_disease_open.bed"
    dc  = f"{TMP_DIR}/{base}_disease_close.bed"

    n_ao  = count_bed(ao);  n_ac  = count_bed(ac)
    n_aor = count_bed(aor); n_acr = count_bed(acr)
    n_do  = count_bed(do);  n_dc  = count_bed(dc)

    if n_do == 0 and n_dc == 0:
        print(f"  SKIP (no disease DARs): {base}"); continue

    oo   = bedtools_overlap(do, ao)
    cc   = bedtools_overlap(dc, ac)
    oc   = bedtools_overlap(do, ac)
    co   = bedtools_overlap(dc, ao)
    oo_r = bedtools_overlap(do, aor)
    cc_r = bedtools_overlap(dc, acr)

    row = dict(
        tissue=tissue, cell_type=ct,
        disease_opening=n_do, disease_closing=n_dc,
        aging_opening=n_ao,   aging_closing=n_ac,
        aging_opening_relaxed=n_aor, aging_closing_relaxed=n_acr,
        overlap_open_open=oo, overlap_close_close=cc,
        overlap_open_close=oc, overlap_close_open=co,
        overlap_open_open_relaxed=oo_r, overlap_close_close_relaxed=cc_r,
        n_tested=n_tested,
        pct_disease_open_in_aging_open=    oo/n_do*100   if n_do else 0,
        pct_disease_close_in_aging_close=  cc/n_dc*100   if n_dc else 0,
        pct_disease_open_in_aging_open_relaxed=  oo_r/n_do*100 if n_do else 0,
        pct_disease_close_in_aging_close_relaxed=cc_r/n_dc*100 if n_dc else 0,
    )
    rows.append(row)
    print(f"  {base}: "
          f"dis_open={n_do} aging_open={n_ao} oo={oo} ({row['pct_disease_open_in_aging_open']:.1f}%) | "
          f"dis_close={n_dc} aging_close={n_ac} cc={cc} ({row['pct_disease_close_in_aging_close']:.1f}%)")

df = pd.DataFrame(rows)
out_csv = f"{OUT_DIR}/DAR_overlap_summary.csv"
df.to_csv(out_csv, index=False)
print(f"\nSaved: {out_csv}")
print(df.to_string())
PYEOF

### Step 3: Enrichment plots (Fisher's exact test) ###
echo ""
echo "$(date)  Step 3: Plotting..."

Rscript - <<'REOF'
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

OUT_DIR <- "/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_science_comparison/results_Lung_5movs21mo"

overlap <- read.csv(file.path(OUT_DIR, "DAR_overlap_summary.csv"),
                    stringsAsFactors = FALSE)

fisher_p <- function(x, K, n, N) {
  if (x == 0 | K == 0 | n == 0 | is.na(N) | N <= 0) return(1)
  a <- x;      b <- n - x
  c <- K - x;  d <- N - n - K + x
  if (b < 0 | c < 0 | d < 0) return(1)
  fisher.test(matrix(c(a, b, c, d), nrow = 2),
              alternative = "greater")$p.value
}

fisher_or <- function(x, K, n, N) {
  if (x == 0 | K == 0 | n == 0 | is.na(N) | N <= 0) return(NA_real_)
  a <- x;      b <- n - x
  c <- K - x;  d <- N - n - K + x
  if (b <= 0 | c <= 0 | d <= 0) return(NA_real_)
  (a * d) / (b * c)
}

sig_label <- function(p) {
  dplyr::case_when(p < 0.001 ~ "***", p < 0.01 ~ "**",
                   p < 0.05  ~ "*",   TRUE       ~ "ns")
}

overlap <- overlap %>%
  rowwise() %>%
  mutate(
    N_use      = n_tested,
    pval_open  = fisher_p(overlap_open_open_relaxed,  aging_opening_relaxed, disease_opening,  N_use),
    pval_close = fisher_p(overlap_close_close_relaxed, aging_closing_relaxed, disease_closing, N_use),
    fold_open  = fisher_or(overlap_open_open_relaxed,  aging_opening_relaxed, disease_opening,  N_use),
    fold_close = fisher_or(overlap_close_close_relaxed, aging_closing_relaxed, disease_closing, N_use),
    sig_open   = sig_label(pval_open),
    sig_close  = sig_label(pval_close),
    label      = paste0(tissue, "\n", cell_type)
  ) %>%
  ungroup()

write.csv(
  overlap %>% select(tissue, cell_type,
    disease_opening, aging_opening_relaxed, overlap_open_open_relaxed,
    fold_open,  pval_open,  sig_open,
    disease_closing, aging_closing_relaxed, overlap_close_close_relaxed,
    fold_close, pval_close, sig_close),
  file.path(OUT_DIR, "overlap_stats.csv"), row.names = FALSE
)
message("Saved: overlap_stats.csv")

KEEP <- c("AT2", "B", "EC-vasc", "Mac-alv", "T")
overlap <- overlap %>% filter(cell_type %in% KEEP)
message("Cell types retained: ", nrow(overlap))

ct_order <- overlap %>% arrange(cell_type) %>% pull(label) %>% unique()
overlap$label <- factor(overlap$label, levels = ct_order)

dir_colors <- c(Opening = "#B2182B", Closing = "#2166AC")

fmt_pval <- function(p) {
  dplyr::case_when(
    is.na(p) | p >= 1   ~ "",
    p == 0              ~ "p<1e-300",
    p < 0.001           ~ paste0("p=", formatC(p, format = "e", digits = 1)),
    TRUE                ~ paste0("p=", round(p, 3))
  )
}

plot_df <- overlap %>%
  select(label, tissue, fold_open, fold_close, sig_open, sig_close,
         pval_open, pval_close) %>%
  pivot_longer(cols = c(fold_open, fold_close),
               names_to = "direction", values_to = "fold") %>%
  mutate(
    direction  = ifelse(direction == "fold_open", "Opening", "Closing"),
    sig        = ifelse(direction == "Opening", sig_open, sig_close),
    pval       = ifelse(direction == "Opening", pval_open, pval_close),
    pval_label = fmt_pval(pval),
    fold_plot  = pmin(fold, 20, na.rm = TRUE)
  )

p_bar <- ggplot(plot_df, aes(x = label, y = fold_plot, fill = direction)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed",
             color = "grey50", linewidth = 0.5) +
  geom_text(aes(label = sig, group = direction),
            position = position_dodge(0.8),
            vjust = -0.3, size = 4.5, fontface = "bold") +
  geom_text(aes(label = pval_label, group = direction),
            position = position_dodge(0.8),
            vjust = -1.6, size = 2.8, color = "grey30") +
  scale_fill_manual(values = dir_colors) +
  labs(
    title    = "Lung - Disease DAR overlap with PURE aging DARs (Aged 21mo vs Adult 5mo)",
    subtitle = "Odds ratio from Fisher's exact test (one-sided, greater)\n* p<0.05  ** p<0.01  *** p<0.001",
    x = NULL, y = "Odds Ratio (OR)", fill = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x   = element_text(angle = 45, hjust = 1, size = 9),
    plot.title    = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey40"),
    legend.position = "top"
  )

ggsave(file.path(OUT_DIR, "figures", "DAR_overlap_enrichment.pdf"),
       p_bar, width = 8, height = 6)
ggsave(file.path(OUT_DIR, "figures", "DAR_overlap_enrichment.png"),
       p_bar, width = 8, height = 6, dpi = 300)
message("Saved: DAR_overlap_enrichment")

p_bar_pval <- ggplot(plot_df, aes(x = label, y = fold_plot, fill = direction)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed",
             color = "grey50", linewidth = 0.5) +
  geom_text(aes(label = pval_label, group = direction),
            position = position_dodge(0.8),
            vjust = -0.4, size = 2.6, color = "grey20") +
  scale_fill_manual(values = dir_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.20))) +
  labs(
    title    = "Lung - Disease DAR overlap with PURE aging DARs (Aged 21mo vs Adult 5mo)",
    subtitle = "Odds ratio from Fisher's exact test (one-sided, greater) | p-values shown",
    x = NULL, y = "Odds Ratio (OR)", fill = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x   = element_text(angle = 45, hjust = 1, size = 9),
    plot.title    = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey40"),
    legend.position = "top"
  )

ggsave(file.path(OUT_DIR, "figures", "DAR_overlap_enrichment_pval.pdf"),
       p_bar_pval, width = 8, height = 6)
ggsave(file.path(OUT_DIR, "figures", "DAR_overlap_enrichment_pval.png"),
       p_bar_pval, width = 8, height = 6, dpi = 300)
message("Saved: DAR_overlap_enrichment_pval")
REOF

### Step 4: Scatter plot - disease-significant peaks vs aging landscape ###
echo ""
echo "$(date)  Step 4: Generating scatter plots..."

PAIRED_DIR="${OUT_DIR}/paired_lfc"
mkdir -p "${PAIRED_DIR}"

for aging_tsv in "${AGING_DIR}"/*_Aged_vs_Adult_DAR.tsv; do
    base=$(basename "${aging_tsv}" _Aged_vs_Adult_DAR.tsv)
    tissue="${base%%_*}"
    [[ "${tissue}" != "${TISSUE}" ]] && continue
    ct="${base#*_}"

    mapped_cts="${CT_MAP[${base}]:-${ct}}"
    tmp_combined="${TMP_DIR}/${base}_disease_combined.tsv"
    if [ ! -f "${tmp_combined}" ]; then
        first_file=true
        for mapped_ct in ${mapped_cts}; do
            f="${DAR_TABLES}/${mapped_ct}__${CONTRAST}__005_DESeq2_all.tsv"
            if [ -f "${f}" ]; then
                if ${first_file}; then
                    cat "${f}" > "${tmp_combined}"; first_file=false
                else
                    tail -n +2 "${f}" >> "${tmp_combined}"
                fi
            fi
        done
        ${first_file} && continue
    fi
    dis_tsv="${tmp_combined}"
    if [ ! -f "${dis_tsv}" ]; then continue; fi

    awk -F'\t' 'NR>1 && $7~/^[0-9eE.+-]/ && $7+0 < 0.05 && $3+0 != 0 {
        split($1, a, "-"); print a[1]"\t"a[2]"\t"a[3]"\t"$3"\t"$7
    }' "${dis_tsv}" | sort -k1,1 -k2,2n \
      > "${TMP_DIR}/${base}_disease_sig_lfc.bed"

    awk -F'\t' 'NR>1 && $2~/^[0-9eE.+-]/ {
        gsub(/:/, "\t", $7); gsub(/-/, "\t", $7)
        print $7"\t"$2"\t"($6~/^[0-9]/ ? $6 : "1")
    }' "${aging_tsv}" | sort -k1,1 -k2,2n \
      > "${TMP_DIR}/${base}_aging_landscape_lfc.bed"

    n_sig=$(wc -l < "${TMP_DIR}/${base}_disease_sig_lfc.bed")
    n_land=$(wc -l < "${TMP_DIR}/${base}_aging_landscape_lfc.bed")
    if [ "${n_sig}" -eq 0 ] || [ "${n_land}" -eq 0 ]; then continue; fi

    bedtools intersect \
        -a "${TMP_DIR}/${base}_disease_sig_lfc.bed" \
        -b "${TMP_DIR}/${base}_aging_landscape_lfc.bed" \
        -wa -wb \
      > "${PAIRED_DIR}/${base}_paired_sig.bed"

    n_paired=$(wc -l < "${PAIRED_DIR}/${base}_paired_sig.bed")
    echo "  ${base}: disease_sig=${n_sig}  aging_landscape=${n_land}  paired=${n_paired}"
done

Rscript - <<'REOF2'
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

PAIRED_DIR <- "/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_science_comparison/results_Lung_5movs21mo/paired_lfc"
OUT_DIR    <- "/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_science_comparison/results_Lung_5movs21mo"

KEEP <- c("AT2", "B", "EC-vasc", "Mac-alv", "T")

files <- list.files(PAIRED_DIR, pattern = "_paired_sig\\.bed$", full.names = TRUE)
if (length(files) == 0) {
  message("No paired_sig BED files found - skipping scatter plots.")
  quit(status = 0)
}

plots   <- list()
r_stats <- list()

for (f in sort(files)) {
  base   <- sub("_paired_sig\\.bed$", "", basename(f))
  tissue <- sub("_.*", "", base)
  ct     <- sub("^[^_]+_", "", base)

  if (!(ct %in% KEEP)) { message("  SKIP (not in whitelist): ", base); next }

  dat <- tryCatch(
    read.table(f, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
               col.names = c("d_chr","d_start","d_end","d_lfc","d_padj",
                             "a_chr","a_start","a_end","a_lfc","a_padj")),
    error = function(e) NULL
  )
  if (is.null(dat) || nrow(dat) < 10) {
    message("  SKIP (too few pairs): ", base); next
  }

  dat <- dat %>%
    mutate(across(c(d_lfc, d_padj, a_lfc, a_padj), as.numeric)) %>%
    filter(!is.na(d_lfc), !is.na(a_lfc))

  if (nrow(dat) < 10) { message("  SKIP (too few valid pairs): ", base); next }

  dat <- dat %>%
    mutate(
      aging_sig = case_when(
        !is.na(a_padj) & a_padj < 0.05  ~ "Aging sig (padj<0.05)",
        !is.na(a_padj) & a_padj < 0.2   ~ "Aging trending (padj<0.2)",
        TRUE                             ~ "Not sig in aging"
      ),
      # Factor order: Not sig FIRST so it plots underneath as grey background;
      # Sig LAST so it plots on top.
      aging_sig = factor(aging_sig, levels = c("Not sig in aging",
                                               "Aging trending (padj<0.2)",
                                               "Aging sig (padj<0.05)"))
    ) %>%
    arrange(aging_sig)

  r_all   <- cor(dat$d_lfc, dat$a_lfc, method = "pearson")
  aging_s <- filter(dat, aging_sig == "Aging sig (padj<0.05)")
  # Relaxed set: sig + trending - used for regression when sig alone too small
  aging_r <- filter(dat, aging_sig %in% c("Aging sig (padj<0.05)",
                                           "Aging trending (padj<0.2)"))

  if (nrow(aging_s) < 5) {
    message("  Note: aging-sig n=", nrow(aging_s), " < 5, r_aging_sig set to NA: ", base)
    r_aging <- NA_real_
  } else {
    r_aging <- cor(aging_s$d_lfc, aging_s$a_lfc, method = "pearson")
  }
  if (nrow(aging_r) < 5) {
    r_aging_relaxed <- NA_real_
  } else {
    r_aging_relaxed <- cor(aging_r$d_lfc, aging_r$a_lfc, method = "pearson")
  }

  dis_open  <- filter(dat, d_lfc > 0)
  dis_close <- filter(dat, d_lfc < 0)
  pct_open_concordant  <- if (nrow(dis_open)  > 0) mean(dis_open$a_lfc  > 0)*100 else NA
  pct_close_concordant <- if (nrow(dis_close) > 0) mean(dis_close$a_lfc < 0)*100 else NA

  label_r <- sprintf(
    paste0("r = %.3f (all, n=%d)\n",
           "r = %s (aging sig, n=%d)\n",
           "r = %s (sig+trend, n=%d)\n",
           "Concordant: open %.0f%% | close %.0f%%"),
    r_all, nrow(dat),
    if (is.na(r_aging))         "NA" else sprintf("%.3f", r_aging), nrow(aging_s),
    if (is.na(r_aging_relaxed)) "NA" else sprintf("%.3f", r_aging_relaxed), nrow(aging_r),
    ifelse(is.na(pct_open_concordant), 0, pct_open_concordant),
    ifelse(is.na(pct_close_concordant), 0, pct_close_concordant)
  )

  sig_colors <- c("Not sig in aging"          = "grey80",
                  "Aging trending (padj<0.2)" = "#F4A582",
                  "Aging sig (padj<0.05)"     = "#B2182B")
  sig_alpha  <- c("Not sig in aging"          = 0.25,
                  "Aging trending (padj<0.2)" = 0.7,
                  "Aging sig (padj<0.05)"     = 0.9)

  p <- ggplot(dat, aes(x = d_lfc, y = a_lfc, color = aging_sig,
                       alpha = aging_sig)) +
    geom_point(size = 0.9) +
    { if (nrow(aging_s) >= 5)
        geom_smooth(data = aging_s, aes(x = d_lfc, y = a_lfc),
                    method = "lm", se = TRUE,
                    color = "#B2182B", fill = "#FDDBC7",
                    linewidth = 0.9, inherit.aes = FALSE)
      else if (nrow(aging_r) >= 5)
        geom_smooth(data = aging_r, aes(x = d_lfc, y = a_lfc),
                    method = "lm", se = TRUE,
                    color = "#F4A582", fill = "#FDDBC7",
                    linewidth = 0.9, linetype = "longdash",
                    inherit.aes = FALSE) } +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    annotate("text", x = -Inf, y = Inf, label = label_r,
             hjust = -0.05, vjust = 1.1, size = 2.3, color = "grey20") +
    scale_color_manual(values = sig_colors, drop = FALSE) +
    scale_alpha_manual(values = sig_alpha, guide = "none") +
    labs(title = paste0(tissue, " - ", ct),
         x = "Disease log2FC (padj<0.05)", y = "Aging log2FC (Aged 21mo vs Adult 5mo, all peaks)",
         color = NULL) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
          legend.position = "bottom", legend.text = element_text(size = 7),
          legend.key.size = unit(0.4, "cm"))

  plots[[base]] <- p
  r_stats[[base]] <- data.frame(
    tissue = tissue, cell_type = ct,
    r_all = r_all, n_all = nrow(dat),
    r_aging_sig = r_aging, n_aging_sig = nrow(aging_s),
    r_aging_relaxed = r_aging_relaxed, n_aging_relaxed = nrow(aging_r),
    pct_open_concordant  = ifelse(is.na(pct_open_concordant),  NA_real_, pct_open_concordant),
    pct_close_concordant = ifelse(is.na(pct_close_concordant), NA_real_, pct_close_concordant)
  )
}

empty_r <- data.frame(tissue=character(), cell_type=character(),
                       r_all=numeric(), n_all=integer(),
                       r_aging_sig=numeric(), n_aging_sig=integer(),
                       r_aging_relaxed=numeric(), n_aging_relaxed=integer(),
                       pct_open_concordant=numeric(), pct_close_concordant=numeric())
r_stats_df <- if (length(r_stats) > 0) do.call(rbind, r_stats) else empty_r
write.csv(r_stats_df, file.path(OUT_DIR, "scatter_r_stats.csv"), row.names = FALSE)
message("Saved: scatter_r_stats.csv")

if (length(plots) == 0) { message("No scatter plots generated."); quit(status = 0) }

for (nm in names(plots)) {
  ggsave(file.path(OUT_DIR, "figures", paste0("scatter_", nm, ".pdf")),
         plots[[nm]], width = 5, height = 5)
}

combined <- wrap_plots(plots, ncol = 2) +
  plot_annotation(
    title    = "Lung - Disease-significant DARs vs PURE aging accessibility (Aged 21mo vs Adult 5mo)",
    subtitle = "X: disease log2FC (padj<0.05) | Y: aging log2FC (all peaks)\nRed=aging sig | Orange=trending | Grey=ns",
    theme = theme(plot.title    = element_text(face = "bold", hjust = 0.5, size = 13),
                  plot.subtitle = element_text(hjust = 0.5, size = 8, color = "grey40"))
  )

ggsave(file.path(OUT_DIR, "figures", "DAR_scatter_combined.pdf"),
       combined, width = 10, height = 10)
ggsave(file.path(OUT_DIR, "figures", "DAR_scatter_combined.png"),
       combined, width = 10, height = 10, dpi = 300)
message("Saved: DAR_scatter_combined (", length(plots), " panels)")
REOF2

### Step 5: Aim 3 summary figure (Fisher OR + Pearson r) ###
echo ""
echo "$(date)  Step 5: Aim 3 summary figure..."

Rscript - <<'REOF3'
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

OUT_DIR <- "/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_science_comparison/results_Lung_5movs21mo"
FIG_DIR <- file.path(OUT_DIR, "figures")

fisher <- read.csv(file.path(OUT_DIR, "overlap_stats.csv"), stringsAsFactors = FALSE)

r_csv <- file.path(OUT_DIR, "scatter_r_stats.csv")
r_df  <- if (file.exists(r_csv) && nrow(read.csv(r_csv)) > 0)
           read.csv(r_csv, stringsAsFactors = FALSE) else NULL

KEEP <- c("AT2", "B", "EC-vasc", "Mac-alv", "T")
fisher <- fisher %>% filter(cell_type %in% KEEP)
if (!is.null(r_df)) r_df <- r_df %>% filter(cell_type %in% KEEP)

sig_label <- function(p) {
  dplyr::case_when(p < 0.001 ~ "***", p < 0.01 ~ "**",
                   p < 0.05  ~ "*",   TRUE       ~ "ns")
}

ct_order <- paste0("Lung\n", c("AT2", "B", "EC-vasc", "Mac-alv", "T"))

fisher_long <- fisher %>%
  mutate(label     = paste0(tissue, "\n", cell_type),
         sig_open  = sig_label(pval_open),
         sig_close = sig_label(pval_close)) %>%
  select(label, tissue, fold_open, fold_close, sig_open, sig_close) %>%
  pivot_longer(cols = c(fold_open, fold_close),
               names_to = "direction", values_to = "OR") %>%
  mutate(direction = ifelse(direction == "fold_open", "Opening", "Closing"),
         sig       = ifelse(direction == "Opening", sig_open, sig_close),
         OR_plot   = pmin(OR, 12, na.rm = TRUE),
         label     = factor(label, levels = ct_order))

dir_colors <- c("Opening" = "#B2182B", "Closing" = "#2166AC")

pA <- ggplot(fisher_long, aes(x = label, y = OR_plot, fill = direction)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_text(aes(label = sig, group = direction),
            position = position_dodge(0.8), vjust = -0.3, size = 4, fontface = "bold") +
  scale_fill_manual(values = dir_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "A  Enrichment of disease DARs in PURE aging DARs (Lung, 5mo->21mo)",
       x = NULL, y = "Enrichment fold", fill = NULL) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        plot.title  = element_text(face = "bold", size = 11),
        legend.position = "top")

if (!is.null(r_df) && nrow(r_df) > 0) {
  r_long <- r_df %>%
    mutate(label = factor(paste0(tissue, "\n", cell_type), levels = ct_order)) %>%
    select(label, tissue, r_all, n_all)

  pB <- ggplot(r_long, aes(x = label, y = r_all)) +
    geom_col(fill = "#4D4D4D", width = 0.55) +
    geom_text(aes(label = sprintf("n=%d", n_all)),
              vjust = -0.4, size = 3, color = "grey30") +
    geom_hline(yintercept = 0, linetype = "solid", color = "grey50", linewidth = 0.4) +
    scale_y_continuous(limits = c(-0.1, 0.65),
                       expand = expansion(mult = c(0.05, 0.12))) +
    labs(title = "B  Pearson r - disease vs PURE aging (5mo->21mo) log2FC (Lung)",
         x = NULL, y = "Pearson r (all disease-sig peaks)") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
          plot.title  = element_text(face = "bold", size = 11))

  combined <- pA / pB
} else {
  message("No scatter r stats available - saving Panel A only")
  combined <- pA
}

combined <- combined +
  plot_annotation(
    title    = "Lung: Overlap between disease DARs and aging DARs",
    subtitle = "Fisher's exact test enrichment | Pearson r on aging-significant peaks",
    theme = theme(plot.title    = element_text(face = "bold", hjust = 0.5, size = 13),
                  plot.subtitle = element_text(hjust = 0.5, size = 8, color = "grey40"))
  )

ggsave(file.path(FIG_DIR, "Aim3_DAR_overlap_summary_Lung.pdf"),
       combined, width = 8, height = 9)
ggsave(file.path(FIG_DIR, "Aim3_DAR_overlap_summary_Lung.png"),
       combined, width = 8, height = 9, dpi = 300)
message("Saved: Aim3_DAR_overlap_summary_Lung")
REOF3

echo "$(date)  Lung all done."

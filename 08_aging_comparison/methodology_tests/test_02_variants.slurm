#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --job-name=test_02_variants
#SBATCH --time=03:00:00
#SBATCH --partition=general
#SBATCH --account=a_imb_ccbcd
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

# Test two changes to step 02 on aging DAR counts:
#   Variant A: add sex covariate (design = ~ sex + condition)
#   Variant C: drop CPM top 25% pre-filter
# Cell types: Kidney PC, TAL; Lung AT2, Mac-alv
#   V0: baseline (current 02 settings)
#   V1: +sex covariate            (A)
#   V2: -CPM top 25%              (C)
#   V3: +sex AND -CPM             (A + C)
# Output: CSV comparing DAR counts across the 4 variants

set -euo pipefail

PYTHON="/home/s4869245/.conda/envs/scanpy_env/bin/python"

${PYTHON} - <<'PYEOF'
import scanpy as sc
import pandas as pd
import numpy as np
import os, warnings, re
warnings.filterwarnings('ignore')

SCIENCE_DIR = "/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_science_comparison"
OUT_CSV     = f"{SCIENCE_DIR}/test_02_variants_results.csv"

### Config - cell types to test ###
TESTS = [
    ("Kidney", "kidney_processed/GSM8774007_Kidney_peak_count.h5ad", "PC",
     ["Principal cells"]),
    ("Kidney", "kidney_processed/GSM8774007_Kidney_peak_count.h5ad", "TAL",
     ["Thick ascending limb of LOH cells"]),
    ("Lung",   "lung_processed/GSM8774006_Lung_peak_count.h5ad", "AT2",
     ["Type II alveolar epithelial cells"]),
    ("Lung",   "lung_processed/GSM8774006_Lung_peak_count.h5ad", "Mac-alv",
     ["Myeloid cells_Alveolar macrophages"]),
]

YOUNG_LABEL = "Adult"       # 5mo
AGED_LABEL  = "Aged"        # 21mo
MIN_CELLS_PER_SAMPLE = 10
MIN_SAMPLES          = 4
MIN_PEAK_FRAC        = 0.05
MIN_MEAN_COUNT       = 1

### Generic pseudo-bulk function ###
def pseudobulk(adata, science_labels, use_cpm_filter=True):
    mask = adata.obs['Main_cell_type'].isin(science_labels)
    ct   = adata[mask].copy()
    ct   = ct[ct.obs['Age'].isin([YOUNG_LABEL, AGED_LABEL])].copy()
    if ct.shape[0] < 50: return None, None

    pb_counts, pb_meta = {}, {}
    for samp in ct.obs['Sample'].unique():
        idx = ct.obs['Sample'] == samp
        n = idx.sum()
        if n < MIN_CELLS_PER_SAMPLE: continue
        counts = np.asarray(ct.X[idx.values].sum(axis=0)).flatten()
        age = ct.obs.loc[idx, 'Age'].iloc[0]
        # Extract sex from the Sample name (e.g. "Adult_Female_3" -> "Female")
        sex_match = re.search(r'(Female|Male)', samp)
        sex = sex_match.group(1) if sex_match else "Unknown"
        pb_counts[samp] = counts
        pb_meta[samp] = {'Age': age, 'sex': sex, 'n_cells': n}

    if len(pb_counts) < MIN_SAMPLES: return None, None
    counts_df = pd.DataFrame(pb_counts, index=ct.var_names).T.astype(int)
    meta_df = pd.DataFrame(pb_meta).T

    # Same frac + mean filters
    min_s = max(2, int(MIN_PEAK_FRAC * len(pb_counts)))
    keep_frac = (counts_df > 0).sum(axis=0) >= min_s
    keep_mean = counts_df.mean(axis=0) >= MIN_MEAN_COUNT

    # CPM top 25% filter (optional)
    if use_cpm_filter:
        total_counts = counts_df.sum(axis=1)
        cpm_df = counts_df.div(total_counts, axis=0) * 1e6
        max_cpm = cpm_df.max(axis=0)
        cpm_thresh = max_cpm.quantile(0.75)
        keep_cpm = max_cpm >= cpm_thresh
        keep = keep_cpm & keep_frac & keep_mean
    else:
        keep = keep_frac & keep_mean

    counts_df = counts_df.loc[:, keep]
    return counts_df, meta_df


def run_deseq2(counts_df, meta_df, use_sex=False):
    """Run DESeq2, optionally adding a sex covariate"""
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    meta_df = meta_df.copy()
    meta_df['condition'] = meta_df['Age'].astype(str)

    if use_sex:
        # Add sex covariate
        meta_use = meta_df[['condition', 'sex']]
        design_factors = ['sex', 'condition']
        ref_level = [['condition', YOUNG_LABEL], ['sex', 'Female']]
    else:
        meta_use = meta_df[['condition']]
        design_factors = 'condition'
        ref_level = ['condition', YOUNG_LABEL]

    try:
        dds = DeseqDataSet(
            counts=counts_df, metadata=meta_use,
            design_factors=design_factors, ref_level=ref_level,
            refit_cooks=True, n_cpus=8
        )
        dds.deseq2()
        stat = DeseqStats(dds, contrast=["condition", AGED_LABEL, YOUNG_LABEL], n_cpus=8)
        stat.summary()
        return stat.results_df
    except Exception as e:
        print(f"    ERROR: {e}")
        return None


def count_sig(res, padj_cut=0.05):
    if res is None: return (0, 0, 0)
    res = res.dropna(subset=['padj', 'log2FoldChange'])
    n_open  = ((res['padj'] < padj_cut) & (res['log2FoldChange'] > 0)).sum()
    n_close = ((res['padj'] < padj_cut) & (res['log2FoldChange'] < 0)).sum()
    n_tested = len(res)
    return (n_tested, n_open, n_close)


### Main loop ###
results = []
adata_cache = {}

for tissue, fn, ct_label, science_cts in TESTS:
    print(f"  {tissue} - {ct_label}")
    # Load h5ad (cache)
    if tissue not in adata_cache:
        print(f"  Loading {tissue} h5ad...")
        adata_cache[tissue] = sc.read_h5ad(f"{SCIENCE_DIR}/{fn}")
        print(f"  Loaded: {adata_cache[tissue].shape}")
    adata = adata_cache[tissue]

    # 4 variants
    for vname, use_sex, use_cpm in [
        ("V0_baseline",            False, True),
        ("V1_sex_covariate",       True,  True),
        ("V2_no_cpm_filter",       False, False),
        ("V3_sex_AND_no_cpm",      True,  False),
    ]:
        print(f"\n  --- {vname} ---")
        counts_df, meta_df = pseudobulk(adata, science_cts, use_cpm_filter=use_cpm)
        if counts_df is None:
            print(f"    SKIP: insufficient samples/cells")
            continue
        print(f"    Pseudo-bulk: {counts_df.shape[0]} samples x {counts_df.shape[1]} peaks")

        # Check sex balance
        n_female = (meta_df['sex'] == 'Female').sum()
        n_male   = (meta_df['sex'] == 'Male').sum()
        print(f"    Sex: Female={n_female}, Male={n_male}")
        # If sex is too imbalanced for V1/V3, fall back to V0 (sex covariate would fail)
        if use_sex and (n_female < 2 or n_male < 2):
            print(f"    WARNING: sex too imbalanced for covariate model, skipping {vname}")
            continue

        res = run_deseq2(counts_df, meta_df, use_sex=use_sex)
        n_tested, n_open, n_close = count_sig(res, padj_cut=0.05)

        results.append({
            'tissue': tissue, 'cell_type': ct_label,
            'variant': vname,
            'use_sex_covariate': use_sex,
            'use_cpm_filter': use_cpm,
            'n_samples': counts_df.shape[0],
            'n_peaks_tested': counts_df.shape[1],
            'n_open_padj05': n_open,
            'n_close_padj05': n_close,
            'n_total_sig': n_open + n_close,
        })
        print(f"    Result: tested={n_tested}, open={n_open}, close={n_close}, total_sig={n_open + n_close}")

### Write comparison table ###
df = pd.DataFrame(results)
df.to_csv(OUT_CSV, index=False)

print(f"  Summary comparison table (saved to {OUT_CSV})")
print(df.to_string(index=False))

# Fold improvement of each cell type relative to baseline
print(f"  Fold improvement relative to V0_baseline")
for (tissue, ct), grp in df.groupby(['tissue', 'cell_type']):
    print(f"\n{tissue} - {ct}:")
    baseline = grp[grp['variant'] == 'V0_baseline']
    if len(baseline) == 0:
        print("  (no baseline)")
        continue
    b_total = baseline['n_total_sig'].iloc[0]
    if b_total == 0:
        print(f"  Baseline n_sig = 0, can't compute ratio")
        continue
    for _, row in grp.iterrows():
        ratio = row['n_total_sig'] / b_total
        print(f"  {row['variant']:25s} n_sig={row['n_total_sig']:>6}   ratio={ratio:.2f}x")
PYEOF

echo ""
echo "$(date)  Test complete. CSV saved to: ${SCIENCE_DIR}/test_02_variants_results.csv"

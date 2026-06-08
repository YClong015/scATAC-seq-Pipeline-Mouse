#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --job-name=test_02_sex
#SBATCH --time=01:00:00
#SBATCH --partition=general
#SBATCH --account=a_imb_ccbcd
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

# Test variant A (sex covariate) vs baseline only; 4 cell types x 2 variants

set -euo pipefail
PYTHON="/home/s4869245/.conda/envs/scanpy_env/bin/python"

${PYTHON} - <<'PYEOF'
import scanpy as sc
import pandas as pd
import numpy as np
import os, warnings, re
warnings.filterwarnings('ignore')

SCIENCE_DIR = "/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_science_comparison"
OUT_CSV     = f"{SCIENCE_DIR}/test_02_sex_only_results.csv"

TESTS = [
    ("Kidney", "kidney_processed/GSM8774007_Kidney_peak_count.h5ad", "PC",
     ["Principal cells"]),
    ("Kidney", "kidney_processed/GSM8774007_Kidney_peak_count.h5ad", "TAL",
     ["Thick ascending limb of LOH cells"]),
    ("Lung", "lung_processed/GSM8774006_Lung_peak_count.h5ad", "AT2",
     ["Type II alveolar epithelial cells"]),
    ("Lung", "lung_processed/GSM8774006_Lung_peak_count.h5ad", "Mac-alv",
     ["Myeloid cells_Alveolar macrophages"]),
]

YOUNG_LABEL, AGED_LABEL = "Adult", "Aged"
MIN_CELLS_PER_SAMPLE, MIN_SAMPLES = 10, 4
MIN_PEAK_FRAC, MIN_MEAN_COUNT = 0.05, 1

def pseudobulk(adata, science_labels):
    """Keep CPM top 25% filter (same as the original 02)"""
    mask = adata.obs['Main_cell_type'].isin(science_labels)
    ct = adata[mask].copy()
    ct = ct[ct.obs['Age'].isin([YOUNG_LABEL, AGED_LABEL])].copy()
    if ct.shape[0] < 50: return None, None

    pb_counts, pb_meta = {}, {}
    for samp in ct.obs['Sample'].unique():
        idx = ct.obs['Sample'] == samp
        n = idx.sum()
        if n < MIN_CELLS_PER_SAMPLE: continue
        counts = np.asarray(ct.X[idx.values].sum(axis=0)).flatten()
        age = ct.obs.loc[idx, 'Age'].iloc[0]
        sex_match = re.search(r'(Female|Male)', samp)
        sex = sex_match.group(1) if sex_match else "Unknown"
        pb_counts[samp] = counts
        pb_meta[samp] = {'Age': age, 'sex': sex, 'n_cells': n}

    if len(pb_counts) < MIN_SAMPLES: return None, None
    counts_df = pd.DataFrame(pb_counts, index=ct.var_names).T.astype(int)
    meta_df = pd.DataFrame(pb_meta).T

    # Keep CPM top 25%
    total_counts = counts_df.sum(axis=1)
    cpm_df = counts_df.div(total_counts, axis=0) * 1e6
    max_cpm = cpm_df.max(axis=0)
    cpm_thresh = max_cpm.quantile(0.75)
    keep_cpm = max_cpm >= cpm_thresh
    min_s = max(2, int(MIN_PEAK_FRAC * len(pb_counts)))
    keep_frac = (counts_df > 0).sum(axis=0) >= min_s
    keep_mean = counts_df.mean(axis=0) >= MIN_MEAN_COUNT
    counts_df = counts_df.loc[:, keep_cpm & keep_frac & keep_mean]
    return counts_df, meta_df

def run_deseq2(counts_df, meta_df, use_sex=False):
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    meta_df = meta_df.copy()
    meta_df['condition'] = meta_df['Age'].astype(str)
    if use_sex:
        meta_use = meta_df[['condition', 'sex']]
        design_factors = ['sex', 'condition']
        ref_level = [['condition', YOUNG_LABEL], ['sex', 'Female']]
    else:
        meta_use = meta_df[['condition']]
        design_factors = 'condition'
        ref_level = ['condition', YOUNG_LABEL]
    try:
        dds = DeseqDataSet(counts=counts_df, metadata=meta_use,
                          design_factors=design_factors, ref_level=ref_level,
                          refit_cooks=True, n_cpus=8)
        dds.deseq2()
        stat = DeseqStats(dds, contrast=["condition", AGED_LABEL, YOUNG_LABEL], n_cpus=8)
        stat.summary()
        return stat.results_df
    except Exception as e:
        print(f"    ERROR: {e}"); return None

results = []
adata_cache = {}

for tissue, fn, ct_label, science_cts in TESTS:
    print(f"\n\n  {tissue} - {ct_label}\n")
    if tissue not in adata_cache:
        print(f"  Loading {tissue} h5ad...")
        adata_cache[tissue] = sc.read_h5ad(f"{SCIENCE_DIR}/{fn}")
        print(f"  Loaded: {adata_cache[tissue].shape}")
    adata = adata_cache[tissue]
    
    counts_df, meta_df = pseudobulk(adata, science_cts)
    if counts_df is None:
        print(f"  SKIP: insufficient samples/cells"); continue
    print(f"  Pseudo-bulk: {counts_df.shape[0]} samples x {counts_df.shape[1]} peaks")
    n_female = (meta_df['sex'] == 'Female').sum()
    n_male = (meta_df['sex'] == 'Male').sum()
    print(f"  Sex: Female={n_female}, Male={n_male}")

    for vname, use_sex in [("V0_baseline", False), ("V1_sex_covariate", True)]:
        if use_sex and (n_female < 2 or n_male < 2):
            print(f"  --- {vname}: SKIPPED (sex too imbalanced)"); continue
        print(f"\n  --- {vname} ---")
        res = run_deseq2(counts_df, meta_df, use_sex=use_sex)
        if res is None: continue
        res = res.dropna(subset=['padj', 'log2FoldChange'])
        n_open = ((res['padj'] < 0.05) & (res['log2FoldChange'] > 0)).sum()
        n_close = ((res['padj'] < 0.05) & (res['log2FoldChange'] < 0)).sum()
        results.append({
            'tissue': tissue, 'cell_type': ct_label, 'variant': vname,
            'use_sex_covariate': use_sex,
            'n_samples': counts_df.shape[0],
            'n_peaks_tested': counts_df.shape[1],
            'n_open': n_open, 'n_close': n_close, 'n_total_sig': n_open + n_close,
        })
        print(f"    Result: open={n_open}, close={n_close}, total_sig={n_open + n_close}")

df = pd.DataFrame(results)
df.to_csv(OUT_CSV, index=False)
print(f"\n\n  Saved to {OUT_CSV}\n\n")
print(df.to_string(index=False))

print(f"\n\n  V1 vs V0 fold improvement\n\n")
for (tissue, ct), grp in df.groupby(['tissue', 'cell_type']):
    v0 = grp[grp['variant'] == 'V0_baseline']
    v1 = grp[grp['variant'] == 'V1_sex_covariate']
    if len(v0) == 0 or len(v1) == 0: continue
    b = v0['n_total_sig'].iloc[0]
    s = v1['n_total_sig'].iloc[0]
    ratio = s/b if b > 0 else float('nan')
    print(f"  {tissue:8s} {ct:10s}  V0={b:>5}  V1={s:>5}  ratio={ratio:.2f}x")
PYEOF

echo "$(date)  Done."

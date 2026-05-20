#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --job-name=pseudobulk_Kidney
#SBATCH --time=24:00:00
#SBATCH --partition=general
#SBATCH --account=a_imb_ccbcd
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

set -euo pipefail

PYTHON="/home/s4869245/.conda/envs/scanpy_env/bin/python"

${PYTHON} - <<'PYEOF'
import scanpy as sc
import pandas as pd
import numpy as np
import os, warnings
warnings.filterwarnings('ignore')

SCIENCE_DIR = "/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_science_comparison"
OUT_DIR     = f"{SCIENCE_DIR}/aging_DARs"
os.makedirs(OUT_DIR, exist_ok=True)

MIN_CELLS_PER_SAMPLE = 10
MIN_SAMPLES          = 4
MIN_PEAK_FRAC        = 0.05
MIN_MEAN_COUNT       = 1

kidney_map = {
    "PT":          ["Proximal tubule cells", "Proximal tubule cells_S3T2"],
    "TAL":         ["Thick ascending limb of LOH cells"],
    "DCT":         ["Distal convoluted tubule cells"],
    "PC":          ["Principal cells"],
    "Macrophages": ["Myeloid cells_Macrophages"],
}

def pseudobulk(adata, ct_label, science_labels):
    mask = adata.obs['Main_cell_type'].isin(science_labels)
    ct   = adata[mask].copy()
    ct   = ct[ct.obs['Age'].isin(['Young', 'Aged'])].copy()

    print(f"    Cells (Young+Aged): {ct.shape[0]}")
    if ct.shape[0] < 50:
        print("    SKIP: too few cells"); return None, None

    pb_counts = {}
    pb_meta   = {}
    for samp in ct.obs['Sample'].unique():
        idx = ct.obs['Sample'] == samp
        n   = idx.sum()
        if n < MIN_CELLS_PER_SAMPLE:
            continue
        counts = np.asarray(ct.X[idx.values].sum(axis=0)).flatten()
        age    = ct.obs.loc[idx, 'Age'].iloc[0]
        pb_counts[samp] = counts
        pb_meta[samp]   = {'Age': age, 'n_cells': n}

    if len(pb_counts) < MIN_SAMPLES:
        print(f"    SKIP: only {len(pb_counts)} pseudo-bulk samples"); return None, None

    counts_df = pd.DataFrame(pb_counts, index=ct.var_names).T.astype(int)
    meta_df   = pd.DataFrame(pb_meta).T

    n_young = (meta_df['Age'] == 'Young').sum()
    n_aged  = (meta_df['Age'] == 'Aged').sum()
    if n_young == 0 or n_aged == 0:
        print("    SKIP: missing one age group"); return None, None

    print(f"    Pseudo-bulk: Young={n_young}, Aged={n_aged}, Peaks={counts_df.shape[1]}")

    # Step 1: CPM pre-filter — keep top 25% most accessible peaks per cell type
    total_counts = counts_df.sum(axis=1)
    cpm_df       = counts_df.div(total_counts, axis=0) * 1e6
    max_cpm      = cpm_df.max(axis=0)
    cpm_thresh   = max_cpm.quantile(0.75)
    keep_cpm     = max_cpm >= cpm_thresh

    # Step 2: fraction filter (non-zero in >=5% of samples)
    min_s      = max(2, int(MIN_PEAK_FRAC * len(pb_counts)))
    keep_frac  = (counts_df > 0).sum(axis=0) >= min_s

    # Step 3: mean count filter
    keep_mean  = counts_df.mean(axis=0) >= MIN_MEAN_COUNT

    keep       = keep_cpm & keep_frac & keep_mean
    counts_df  = counts_df.loc[:, keep]
    print(f"    Peaks after CPM-top25% + frac>={min_s} + mean>={MIN_MEAN_COUNT}: {counts_df.shape[1]}")

    if counts_df.shape[1] < 100:
        print("    SKIP: too few peaks"); return None, None

    return counts_df, meta_df

def run_deseq2(counts_df, meta_df, tissue, ct_label, out_dir):
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds  import DeseqStats

    meta_df = meta_df.copy()
    meta_df['condition'] = meta_df['Age'].map({'Young': 'Young', 'Aged': 'Aged'})

    try:
        dds = DeseqDataSet(
            counts         = counts_df,
            metadata       = meta_df[['condition']],
            design_factors = "condition",
            ref_level      = ["condition", "Young"],
            refit_cooks    = True,
            n_cpus         = 8
        )
        dds.deseq2()

        stat = DeseqStats(dds, contrast=["condition", "Aged", "Young"], n_cpus=8)
        stat.summary()
        res = stat.results_df.copy()
        res['peak']      = res.index
        res['cell_type'] = ct_label
        res['tissue']    = tissue

        out_file = os.path.join(out_dir, f"{tissue}_{ct_label}_Aged_vs_Young_DAR.tsv")
        res.to_csv(out_file, sep='\t', index=False)

        n_open  = ((res['padj'] < 0.05) & (res['log2FoldChange'] > 0)).sum()
        n_close = ((res['padj'] < 0.05) & (res['log2FoldChange'] < 0)).sum()
        print(f"    Saved: {out_file}")
        print(f"    Significant (padj<0.05): opening={n_open}  closing={n_close}")
        return res

    except Exception as e:
        print(f"    ERROR in DESeq2: {e}"); return None

# ── Kidney ───────────────────────────────────────────────────
print("\n" + "="*60)
print("Loading Kidney h5ad...")
adata_k = sc.read_h5ad(
    f"{SCIENCE_DIR}/kidney_processed/GSM8774007_Kidney_peak_count.h5ad"
)
print(f"Shape: {adata_k.shape}")

print("\nAll Kidney Main_cell_type labels:")
for ct in sorted(adata_k.obs['Main_cell_type'].unique()):
    n = (adata_k.obs['Main_cell_type'] == ct).sum()
    print(f"  {n:>7,}  {ct}")
print("\nKidney Age breakdown:")
print(adata_k.obs['Age'].value_counts().to_string())

for our_ct, science_cts in kidney_map.items():
    print(f"\n--- Kidney: {our_ct} ---")
    counts_df, meta_df = pseudobulk(adata_k, our_ct, science_cts)
    if counts_df is not None:
        run_deseq2(counts_df, meta_df, "Kidney", our_ct, OUT_DIR)

print("\nKidney done.")
PYEOF

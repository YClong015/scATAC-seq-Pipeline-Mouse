#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --job-name=explore_h5ad
#SBATCH --time=00:30:00
#SBATCH --partition=general
#SBATCH --account=a_imb_ccbcd
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

set -euo pipefail

PYTHON="/home/s4869245/.conda/envs/scanpy_env/bin/python"

${PYTHON} - <<'PYEOF'
import scanpy as sc
import numpy as np

SCIENCE_DIR = "/QRISdata/Q8448/Mouse_disease_data/DAR/DAR_science_comparison"

for tissue, path in [
    ("Kidney", f"{SCIENCE_DIR}/kidney_processed/GSM8774007_Kidney_peak_count.h5ad"),
    ("Lung",   f"{SCIENCE_DIR}/lung_processed/GSM8774006_Lung_peak_count.h5ad"),
]:
    print(f"  {tissue}: {path}")
    adata = sc.read_h5ad(path)
    print(f"\nShape: {adata.shape[0]} cells x {adata.shape[1]} peaks")

    print(f"\n--- obs columns ---")
    for col in adata.obs.columns:
        n_unique = adata.obs[col].nunique()
        vals = adata.obs[col].unique()[:5]
        print(f"  {col:<30} ({n_unique} unique)  e.g. {list(vals)}")

    print(f"\n--- Peak ID format (first 5 var_names) ---")
    print(f"  {list(adata.var_names[:5])}")

    print(f"\n--- X matrix type ---")
    print(f"  {type(adata.X)}  dtype={adata.X.dtype}")

    del adata

print("\nDone.")
PYEOF

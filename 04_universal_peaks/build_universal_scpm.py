#!/usr/bin/env python3
"""
Build universal peak set from narrowPeak files using pycisTopic
get_consensus_peaks with proper SCPM scoring — two-stage merge:
  Stage 1: per-tissue consensus (cell types → tissue peaks)
  Stage 2: cross-tissue universal (tissue peaks → universal peaks)
  Stage 3: SCPM > 1 filter

Usage (set env vars before running):
  export CHROMSIZES_FILE=/path/to/mm10.chrom.sizes
  export BLACKLIST=/scratch/user/s4869245/pycisTopic/blacklist/mm10-blacklist.v2.bed
  export OUTDIR=/QRISdata/Q8448/Mouse_disease_data/universal_peaks_v4
  python build_universal_scpm.py
"""

import os, sys
import pandas as pd
import pyranges as pr
from pycisTopic.iterative_peak_calling import get_consensus_peaks

BLACKLIST   = os.environ["BLACKLIST"]
OUTDIR      = os.environ["OUTDIR"]
PEAK_HW     = int(os.environ.get("PEAK_HALF_WIDTH", "250"))
SCPM_CUTOFF = float(os.environ.get("SCPM_CUTOFF", "1.0"))

os.makedirs(OUTDIR, exist_ok=True)

# ── Chromsizes ────────────────────────────────────────────────
chrom_local = os.environ.get("CHROMSIZES_FILE", "")
if chrom_local and os.path.exists(chrom_local):
    chromsizes = pd.read_table(chrom_local, header=None,
                               names=["Chromosome", "End"])
else:
    chromsizes = pd.read_table(
        "http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes",
        header=None, names=["Chromosome", "End"])
chromsizes.insert(1, "Start", 0)
valid_chroms = set(chromsizes["Chromosome"].values)
print(f"[OK] Loaded chromsizes: {len(valid_chroms)} chromosomes")

# ── NarrowPeak files per tissue ───────────────────────────────
PEAKS = {
    "Kidney": {
        "IC":       "/QRISdata/Q8448/Mouse_disease_data/Kidney/peaks/IC/sample_peak/IC_peaks.narrowPeak",
        "PC_URO":   "/QRISdata/Q8448/Mouse_disease_data/Kidney/peaks/PC_URO/sample_peak/PC_URO_peaks.narrowPeak",
        "DCT_CNT":  "/QRISdata/Q8448/Mouse_disease_data/Kidney/peaks/DCT_CNT/sample_peak/DCT_CNT_peaks.narrowPeak",
        "LEUK":     "/QRISdata/Q8448/Mouse_disease_data/Kidney/peaks/LEUK/sample_peak/LEUK_peaks.narrowPeak",
        "DTL_ATL":  "/QRISdata/Q8448/Mouse_disease_data/Kidney/peaks/DTL_ATL/sample_peak/DTL_ATL_peaks.narrowPeak",
        "PCT":      "/QRISdata/Q8448/Mouse_disease_data/Kidney/peaks/PCT/sample_peak/PCT_peaks.narrowPeak",
        "PST":      "/QRISdata/Q8448/Mouse_disease_data/Kidney/peaks/PST/sample_peak/PST_peaks.narrowPeak",
        "FIB":      "/QRISdata/Q8448/Mouse_disease_data/Kidney/peaks/FIB/sample_peak/FIB_peaks.narrowPeak",
        "EC":       "/QRISdata/Q8448/Mouse_disease_data/Kidney/peaks/EC/sample_peak/EC_peaks.narrowPeak",
        "PODO_PEC": "/QRISdata/Q8448/Mouse_disease_data/Kidney/peaks/PODO_PEC/sample_peak/PODO_PEC_peaks.narrowPeak",
    },
    "Lung": {
        "AT2":              "/QRISdata/Q8448/Mouse_disease_data/Lung/peaks/AT2/sample_peak/AT2_peaks.narrowPeak",
        "Ciliated":         "/QRISdata/Q8448/Mouse_disease_data/Lung/peaks/Ciliated/sample_peak/Ciliated_peaks.narrowPeak",
        "Mo_Ly6c+":         "/QRISdata/Q8448/Mouse_disease_data/Lung/peaks/Mo_Ly6c+/sample_peak/Mo_Ly6c+_peaks.narrowPeak",
        "Endothelial_cells":"/QRISdata/Q8448/Mouse_disease_data/Lung/peaks/Endothelial_cells/sample_peak/Endothelial_cells_peaks.narrowPeak",
        "Mac_inter":        "/QRISdata/Q8448/Mouse_disease_data/Lung/peaks/Mac_inter/sample_peak/Mac_inter_peaks.narrowPeak",
        "Fib":              "/QRISdata/Q8448/Mouse_disease_data/Lung/peaks/Fib/sample_peak/Fib_peaks.narrowPeak",
        "Pen":              "/QRISdata/Q8448/Mouse_disease_data/Lung/peaks/Pen/sample_peak/Pen_peaks.narrowPeak",
        "SMCs":             "/QRISdata/Q8448/Mouse_disease_data/Lung/peaks/SMCs/sample_peak/SMCs_peaks.narrowPeak",
        "Mac_alv":          "/QRISdata/Q8448/Mouse_disease_data/Lung/peaks/Mac_alv/sample_peak/Mac_alv_peaks.narrowPeak",
        "B":                "/QRISdata/Q8448/Mouse_disease_data/Lung/peaks/B/sample_peak/B_peaks.narrowPeak",
    },
    "Aorta": {
        "Endothelial": "/QRISdata/Q8448/Mouse_disease_data/Aorta/peaks/Endothelial/sample_peak/Endothelial_peaks.narrowPeak",
        "Pericyte":    "/QRISdata/Q8448/Mouse_disease_data/Aorta/peaks/Pericyte/sample_peak/Pericyte_peaks.narrowPeak",
        "SMC":         "/QRISdata/Q8448/Mouse_disease_data/Aorta/peaks/SMC/sample_peak/SMC_peaks.narrowPeak",
        "T-cell":      "/QRISdata/Q8448/Mouse_disease_data/Aorta/peaks/T-cell/sample_peak/T-cell_peaks.narrowPeak",
        "Mac":         "/QRISdata/Q8448/Mouse_disease_data/Aorta/peaks/Mac/sample_peak/Mac_peaks.narrowPeak",
        "Fibroblast":  "/QRISdata/Q8448/Mouse_disease_data/Aorta/peaks/Fibroblast/sample_peak/Fibroblast_peaks.narrowPeak",
    },
    "Tcells": {
        "Naive_T":         "/QRISdata/Q8448/Mouse_disease_data/Tcells/peaks/Naive_T/sample_peak/Naive_T_peaks.narrowPeak",
        "Naive_CD8_T":     "/QRISdata/Q8448/Mouse_disease_data/Tcells/peaks/Naive_CD8_T/sample_peak/Naive_CD8_T_peaks.narrowPeak",
        "Effector_CD8_T":  "/QRISdata/Q8448/Mouse_disease_data/Tcells/peaks/Effector_CD8_T/sample_peak/Effector_CD8_T_peaks.narrowPeak",
        "Cytotoxic_CD8_T": "/QRISdata/Q8448/Mouse_disease_data/Tcells/peaks/Cytotoxic_CD8_T/sample_peak/Cytotoxic_CD8_T_peaks.narrowPeak",
        "CD8_Eff":         "/QRISdata/Q8448/Mouse_disease_data/Tcells/peaks/CD8_Eff/sample_peak/CD8_Eff_peaks.narrowPeak",
        "Memory_CD8_T":    "/QRISdata/Q8448/Mouse_disease_data/Tcells/peaks/Memory_CD8_T/sample_peak/Memory_CD8_T_peaks.narrowPeak",
        "Treg":            "/QRISdata/Q8448/Mouse_disease_data/Tcells/peaks/Treg/sample_peak/Treg_peaks.narrowPeak",
        "Tfh_like_T":      "/QRISdata/Q8448/Mouse_disease_data/Tcells/peaks/Tfh_like_T/sample_peak/Tfh_like_T_peaks.narrowPeak",
        "NK":              "/QRISdata/Q8448/Mouse_disease_data/Tcells/peaks/NK/sample_peak/NK_peaks.narrowPeak",
    },
}

# ── Helper: load narrowPeak ───────────────────────────────────
def load_narrowpeak(path):
    df = pd.read_csv(path, sep="\t", header=None,
                     usecols=range(10),
                     names=["Chromosome","Start","End","Name","Score",
                            "Strand","FC_summit","-log10_pval",
                            "-log10_qval","Summit"])
    df = df[df["Chromosome"].isin(valid_chroms)].copy()
    df["Start"]  = df["Start"].astype(int)
    df["End"]    = df["End"].astype(int)
    df["Summit"] = df["Summit"].astype(int)
    # Use -log10(qval) as the score for SCPM (same as pycisTopic convention)
    df["Score"]  = df["-log10_qval"].clip(lower=0)
    return pr.PyRanges(df)

# ── Stage 1: per-tissue consensus ─────────────────────────────
print("\n=== Stage 1: Per-tissue consensus ===")
tissue_consensus = {}

for tissue, ct_dict in PEAKS.items():
    cached = os.path.join(OUTDIR, f"{tissue}_consensus.bed.gz")
    if os.path.exists(cached):
        print(f"  [{tissue}] Loading cached consensus: {cached}")
        cons_gr = pr.read_bed(cached)
        df_c = cons_gr.df.copy()
        df_c = df_c[df_c["Chromosome"].isin(valid_chroms)].copy()
        df_c["Start"] = df_c["Start"].astype(int)
        df_c["End"]   = df_c["End"].astype(int)
        if "Score" not in df_c.columns:
            df_c["Score"] = 1.0
        tissue_consensus[tissue] = pr.PyRanges(df_c)
        print(f"    → {len(tissue_consensus[tissue])} peaks")
        continue

    print(f"\n  [{tissue}] Loading {len(ct_dict)} cell types...")
    narrow_peaks = {}
    for ct, path in ct_dict.items():
        if not os.path.exists(path):
            print(f"    SKIP (not found): {path}"); continue
        pr_obj = load_narrowpeak(path)
        if len(pr_obj) == 0:
            print(f"    SKIP (empty): {ct}"); continue
        narrow_peaks[ct] = pr_obj
        print(f"    {ct}: {len(pr_obj)} peaks")

    if not narrow_peaks:
        print(f"  SKIP {tissue}: no peaks found"); continue

    consensus = get_consensus_peaks(
        narrow_peaks_dict   = narrow_peaks,
        peak_half_width     = PEAK_HW,
        chromsizes          = chromsizes,
        path_to_blacklist   = BLACKLIST,
    )
    tissue_consensus[tissue] = consensus
    consensus.to_bed(path=cached, keep=True, compression="gzip", chain=False)
    print(f"  [{tissue}] → {len(consensus)} consensus peaks → {cached}")

# ── Stage 2: cross-tissue universal merge ─────────────────────
# Use simple cluster-based merge (equivalent to bedtools merge -d 0).
# get_consensus_peaks is intentionally NOT used here because we have only one
# consensus per tissue — its iterative overlap-removal would discard
# tissue-specific peaks (biologically meaningful signal).
print("\n=== Stage 2: Cross-tissue universal merge (cluster merge) ===")

# Concatenate all tissue consensus peaks, preserving Score.
# Narrow each 500bp peak to ±125bp around its centre before cross-tissue merge
# to reduce false sharing caused by the wide Stage-1 peaks.
NARROW_HW = int(os.environ.get("NARROW_HALF_WIDTH", "125"))
print(f"  Narrowing tissue peaks to ±{NARROW_HW}bp before merge")

df_parts = []
tissue_narrow_gr = {}   # narrowed peaks per tissue — reused for annotation
for tissue, cons_gr in tissue_consensus.items():
    df_tmp = cons_gr.df[["Chromosome", "Start", "End", "Score"]].copy()
    mid = ((df_tmp["Start"] + df_tmp["End"]) / 2).astype(int)
    df_tmp["Start"] = (mid - NARROW_HW).clip(lower=0)
    df_tmp["End"]   = mid + NARROW_HW
    df_tmp["tissue"] = tissue
    tissue_narrow_gr[tissue] = pr.PyRanges(df_tmp[["Chromosome","Start","End"]].copy())
    df_parts.append(df_tmp)

df_all = pd.concat(df_parts, ignore_index=True)
all_gr  = pr.PyRanges(df_all)

# Remove blacklist
bl_gr  = pr.read_bed(BLACKLIST)
all_gr = all_gr.subtract(bl_gr)

# Clean up after subtract: drop NaN coords and degenerate intervals
df_all2 = all_gr.df.copy()
df_all2 = df_all2.dropna(subset=["Start", "End"])
df_all2["Start"] = df_all2["Start"].astype(int)
df_all2["End"]   = df_all2["End"].astype(int)
df_all2 = df_all2[df_all2["Start"] < df_all2["End"]].copy()
all_gr  = pr.PyRanges(df_all2)

# Cluster overlapping peaks; for each cluster keep merged coords + max Score
clustered = all_gr.cluster(slack=0)
df_cl = clustered.df.copy()
df_cl = df_cl.dropna(subset=["Start", "End", "Cluster"])
df_merged = (
    df_cl.groupby(["Chromosome", "Cluster"], sort=False)
    .agg(Start=("Start", "min"), End=("End", "max"), Score=("Score", "max"))
    .reset_index()
    .drop(columns="Cluster")
)
df_merged = df_merged.dropna(subset=["Start", "End"])
df_merged["Start"] = df_merged["Start"].astype(int)
df_merged["End"]   = df_merged["End"].astype(int)
universal = pr.PyRanges(df_merged)
print(f"  Universal peaks (before SCPM filter): {len(universal)}")

out_pre = os.path.join(OUTDIR, "universal_peaks_prescpm.bed.gz")
universal.to_bed(path=out_pre, keep=True, compression="gzip", chain=False)

# ── Stage 3: SCPM > cutoff filter ─────────────────────────────
print(f"\n=== Stage 3: SCPM distribution & filter (cutoff={SCPM_CUTOFF}) ===")

df_univ = universal.df.copy()

if "Score" in df_univ.columns:
    total_score = df_univ["Score"].sum()
    df_univ["SCPM"] = df_univ["Score"] / total_score * 1e6

    # Print distribution to help choose cutoff
    import numpy as np
    scpm = df_univ["SCPM"]
    print(f"  Total peaks:  {len(scpm)}")
    print(f"  SCPM min:     {scpm.min():.4f}")
    print(f"  SCPM median:  {scpm.median():.4f}")
    print(f"  SCPM mean:    {scpm.mean():.4f}")
    print(f"  SCPM max:     {scpm.max():.4f}")
    print(f"  Percentiles:")
    for p in [10, 25, 50, 75, 90, 95, 99]:
        print(f"    {p:3d}th: {np.percentile(scpm, p):.4f}  "
              f"→ retains {(scpm > np.percentile(scpm, p)).sum()} peaks "
              f"({100-p}%)")
    print(f"\n  Peaks retained at different cutoffs:")
    for cutoff in [0.5, 1.0, 2.0, 5.0]:
        n = (scpm > cutoff).sum()
        print(f"    SCPM>{cutoff:.1f}: {n} peaks ({n/len(scpm)*100:.1f}%)")

    df_filt = df_univ[df_univ["SCPM"] > SCPM_CUTOFF].copy()
    print(f"\n  --> Using SCPM>{SCPM_CUTOFF}: {len(df_filt)} peaks retained "
          f"({len(df_filt)/len(df_univ)*100:.1f}%)")
else:
    print("  WARNING: No Score column — skipping SCPM filter")
    df_filt = df_univ

# ── Save final universal peak set ─────────────────────────────
df_out = df_filt[["Chromosome","Start","End"]].copy()
df_out["Name"] = (df_out["Chromosome"].astype(str) + ":" +
                  df_out["Start"].astype(str) + "-" +
                  df_out["End"].astype(str))

gr_out = pr.PyRanges(df_filt)
out_final = os.path.join(OUTDIR, "consensus_regions_v5.bed.gz")
gr_out.to_bed(path=out_final, keep=True, compression="gzip", chain=False)
print(f"\n  Final universal peak set: {out_final}")

# ── Tissue membership for UpSet plot ──────────────────────────
print("\n=== Annotating tissue membership ===")

def count_overlaps(univ_gr, tissue_gr):
    if tissue_gr is None or len(tissue_gr) == 0:
        return [0] * len(univ_gr.df)
    joined = univ_gr.count_overlaps(tissue_gr)
    return (joined.df["NumberOverlaps"] > 0).astype(int).tolist()

univ_gr = pr.PyRanges(df_filt[["Chromosome","Start","End"]].copy())

membership = {}
for tissue, cons in tissue_narrow_gr.items():
    membership[tissue] = count_overlaps(univ_gr, cons)
    print(f"  {tissue}: {sum(membership[tissue])} peaks")

upset_df = df_filt[["Chromosome","Start","End"]].copy()
for tissue in ["Kidney","Lung","Aorta","Tcells"]:
    upset_df[tissue] = membership.get(tissue, [0]*len(upset_df))
upset_df["n_tissues"] = upset_df[["Kidney","Lung","Aorta","Tcells"]].sum(axis=1)

upset_csv = os.path.join(OUTDIR, "upset_input.csv")
upset_df.rename(columns={"Chromosome":"chr","Start":"start","End":"end"}) \
        .to_csv(upset_csv, index=False)
print(f"\n  UpSet CSV: {upset_csv}")

# ── Summary statistics ─────────────────────────────────────────
print("\n=== Summary ===")
for n in [1,2,3,4]:
    print(f"  Shared {n} tissue(s): {(upset_df['n_tissues']==n).sum()}")
    if n == 1:
        for t in ["Kidney","Lung","Aorta","Tcells"]:
            n1 = ((upset_df[t]==1) & (upset_df["n_tissues"]==1)).sum()
            print(f"    {t}-only: {n1}")

print("\n[DONE]")

import os
import sys
import pandas as pd
import pyranges as pr
from pycisTopic.iterative_peak_calling import get_consensus_peaks

def load_chromsizes_mm10() -> pd.DataFrame:
    chrom_local = os.environ.get("CHROMSIZES_FILE", "")
    if chrom_local and os.path.exists(chrom_local):
        chromsizes = pd.read_table(
            chrom_local,
            header=None,
            names=["Chromosome", "End"],
        )
    else:
        url = (
            "http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/"
            "mm10.chrom.sizes"
        )
        chromsizes = pd.read_table(
            url,
            header=None,
            names=["Chromosome", "End"],
        )
    chromsizes.insert(1, "Start", 0)
    return chromsizes

def load_consensus_bed(path: str) -> pr.PyRanges:
    try:
        df = pd.read_csv(
            path,
            sep="\t",
            header=None,
            comment="#",
            compression="gzip",
            usecols=[0, 1, 2],
            names=["Chromosome", "Start", "End"],
            dtype={"Chromosome": str},
        )
    except Exception:
        df = pd.read_csv(
            path,
            sep="\t",
            header=None,
            comment="#",
            usecols=[0, 1, 2],
            names=["Chromosome", "Start", "End"],
            dtype={"Chromosome": str},
        )

    df = df.dropna(subset=["Chromosome", "Start", "End"])
    df = df[df["Chromosome"] != ""]
    df = df[~df["Chromosome"].str.startswith(("track", "browser"))]

    df["Start"] = pd.to_numeric(df["Start"], errors="coerce").fillna(-1).astype(int)
    df["End"] = pd.to_numeric(df["End"], errors="coerce").fillna(-1).astype(int)
    df = df[(df["Start"] >= 0) & (df["End"] > df["Start"])]
    df = df.drop_duplicates(subset=["Chromosome", "Start", "End"]).reset_index(drop=True)

    if df.empty:
        raise ValueError(f"No valid peaks found in: {path}")

    width = df["End"] - df["Start"]
    summit = (width // 2).astype(int)

    df["Name"] = df["Chromosome"] + ":" + df["Start"].astype(str) + "-" + df["End"].astype(str)
    df["Score"] = width.astype(float) + 100.0
    df["Strand"] = "."
    df["SignalValue"] = 10.0
    df["pValue"] = 1.0
    df["qValue"] = 1.0
    df["Summit"] = summit

    return pr.PyRanges(df)

def main() -> None:
    try:
        peakpath_file = os.environ["PEAKPATH"]
        outdir = os.environ["OUTDIR"]
        peak_half_width = int(os.environ.get("PEAK_HALF_WIDTH", "250"))
        blacklist_path = os.environ.get("BLACKLIST")
    except KeyError as e:
        sys.exit(1)
    
    os.makedirs(outdir, exist_ok=True)

    chromsizes = load_chromsizes_mm10()
    valid_chroms = set(chromsizes["Chromosome"].values)

    peaks = {}
    with open(peakpath_file, "r") as f:
        for line in f:
            if not line.strip(): continue
            parts = line.strip().split(":", 1)
            if len(parts) != 2: continue
            k, v = parts[0].strip(), parts[1].strip()
            
            if not os.path.exists(v):
                raise FileNotFoundError(f"Missing: {v}")
            
            pr_obj = load_consensus_bed(v)
            pr_obj = pr_obj[pr_obj.Chromosome.isin(valid_chroms)]
            
            if len(pr_obj) == 0:
                continue
                
            peaks[k] = pr_obj
            print(f"[OK] Loaded {k}: {len(peaks[k])} peaks")

    if not peaks:
        sys.exit(1)

    consensus = get_consensus_peaks(
        narrow_peaks_dict=peaks,
        peak_half_width=peak_half_width,
        chromsizes=chromsizes,
        path_to_blacklist=blacklist_path,
    )

    out_path = os.path.join(outdir, "consensus_regions.bed.gz")
    consensus.to_bed(
        path=out_path,
        keep=True,
        compression="gzip",
        chain=False,
    )
    print(f"[SUCCESS] Done! Wrote: {out_path}")

if __name__ == "__main__":
    main()

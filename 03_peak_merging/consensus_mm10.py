import os
import pandas as pd
import pyranges as pr
from pycisTopic.iterative_peak_calling import get_consensus_peaks


def load_narrowpeak(path: str) -> pr.PyRanges:
    df = pd.read_csv(path, sep="\t", header=None)
    df = df.iloc[:, :10].copy()
    df.columns = [
        "Chromosome",
        "Start",
        "End",
        "Name",
        "Score",
        "Strand",
        "FC_summit",
        "-log10_pval",
        "-log10_qval",
        "Summit",
    ]
    return pr.PyRanges(df)


def load_chromsizes_mm10() -> pd.DataFrame:
    chrom_local = os.environ.get("CHROMSIZES_FILE", "")
    if chrom_local and os.path.exists(chrom_local):
        chromsizes = pd.read_table(
            chrom_local,
            header=None,
            names=["Chromosome", "End"],
        )
    else:
        chromsizes = pd.read_table(
            "http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/"
            "mm10.chrom.sizes",
            header=None,
            names=["Chromosome", "End"],
        )
    chromsizes.insert(1, "Start", 0)
    return chromsizes


def main() -> None:
    peakpath = os.environ["PEAKPATH"]
    outdir = os.environ["OUTDIR"]
    os.makedirs(outdir, exist_ok=True)

    peak_half_width = int(os.environ.get("PEAK_HALF_WIDTH", "250"))
    blacklist = os.environ["BLACKLIST"]

    narrow_peak_dict = {}
    with open(peakpath, "r") as f:
        for line in f:
            k, v = line.strip().split(": ", 1)
            narrow_peak_dict[k] = load_narrowpeak(v)

    chromsizes = load_chromsizes_mm10()

    consensus_peaks = get_consensus_peaks(
        narrow_peaks_dict=narrow_peak_dict,
        peak_half_width=peak_half_width,
        chromsizes=chromsizes,
        path_to_blacklist=blacklist,
    )

    out_path = os.path.join(outdir, "consensus_regions.bed.gz")
    consensus_peaks.to_bed(
        path=out_path,
        keep=True,
        compression="gzip",
        chain=False,
    )

    print(f"[OK] wrote: {out_path}")


if __name__ == "__main__":
    main()


import sys
import pandas as pd
import pyranges as pr
import os
from pycisTopic.pseudobulk_peak_calling import peak_calling
from pycisTopic.iterative_peak_calling import get_consensus_peaks

out_dir="/QRISdata/Q8448/Mouse_disease_data/Kidney/peaks/SRR27367347"

bed_paths = {}
with open(os.path.join(out_dir, "bed_path.tsv")) as f:
    for line in f:
        v, p = line.strip().split("\t")
        bed_paths.update({v: p})

macs_path = "/home/s4869245/.conda/envs/scenicplus/bin/macs2"

os.makedirs(os.path.join(out_dir, "consensus_peak_calling/MACS"), exist_ok = True)

temp_dir = os.environ.get('TMPDIR', '/tmp')
os.makedirs(temp_dir, exist_ok=True)

narrow_peak_dict = peak_calling(
    macs_path = macs_path,
    bed_paths = bed_paths,
    outdir = os.path.join(os.path.join(out_dir, "consensus_peak_calling/MACS")),
    genome_size = '2.7e9',
    n_cpu = 10,
    input_format = 'BED',
    keep_dup = 'all',
    q_value = 0.05,
    _temp_dir = temp_dir
)
# Other parameters
chromsizes = pd.read_table(
    "http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes"   ,
    header = None,
    names = ["Chromosome", "End"]
)
chromsizes.insert(1, "Start", 0)
peak_half_width=250
path_to_blacklist="/scratch/user/s4869245/pycisTopic/blacklist/mm10-blacklist.v2.bed"
# Get consensus peaks
consensus_peaks = get_consensus_peaks(
    narrow_peaks_dict = narrow_peak_dict,
    peak_half_width = peak_half_width,
    chromsizes = chromsizes,
    path_to_blacklist = path_to_blacklist)

consensus_peaks.to_bed(
    path = os.path.join(out_dir, "consensus_peak_calling/consensus_regions.bed"),
    keep =True,
    compression = 'gzip',
    chain = False)

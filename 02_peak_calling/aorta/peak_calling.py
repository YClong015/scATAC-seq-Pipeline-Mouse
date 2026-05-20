import os
import pycisTopic

bed_paths = {}
with open(os.environ['BEDPATH']) as f:
    for line in f:
        v, p = line.strip().split(": ")
        bed_paths.update({v: p})

from pycisTopic.pseudobulk_peak_calling import peak_calling
macs_path = "/sw/local/rocky8/noarch/qcif/software/miniconda3/envs/macs2_2.2.9.1/bin/macs2"

narrow_peak_dict = peak_calling(
    macs_path = macs_path,
    bed_paths = bed_paths,
    outdir = os.environ['OUTDIR'],
    genome_size = 'mm',
    input_format = 'BEDPE',
    n_cpu = 10,
    shift = 73,
    ext_size = 146,
    keep_dup = 'all',
    q_value = 0.05,
    _temp_dir = os.environ['TMPDIR']
)


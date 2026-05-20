# Lung COPD raw data download — CNGBdb CNP0004399

Source paper: **Zhang Q. et al., 2025 PLOS ONE** (`paper/Lung_mice_paper.pdf`)
Project: https://db.cngb.org/search/project/CNP0004399/

## Sample mapping

| Group | Sample ID | CNGB run ID | CNGB BGI library |
|---|---|---|---|
| Control | Control_F2 | CL100168054_L01 | CNX0739876/CNR0841132 |
| Control | Control_M1 | CL100167942_L01 | CNX0739878/CNR0841134 |
| Case | Case_F1 | CL100168078_L02 | CNX0739875/CNR0841131 |
| Case | Case_F3 | CL100168054_L02 | CNX0739877/CNR0841133 |
| Case | Case_M2 | CL100167942_L02 | CNX0739879/CNR0841135 |
| Case | Case_M3 | CL100168078_L01 | CNX0739880/CNR0841136 |

## Two download scripts

### `download_full_directory.sh` (recursive — recommended)

Recursive `wget` of a whole CNGB sample directory (gets R1 + R2 + metadata in one shot). Edit `FTP_DIR_URL` for each sample.

```bash
sbatch 00_download/lung_cngb/download_full_directory.sh
```

### `download_per_sample.sh` (single URL, fallback)

Targeted single-file `wget` if recursive download fails (e.g. CNGB rate-limiting). Edit the `URLS` array.

```bash
sbatch 00_download/lung_cngb/download_per_sample.sh
```

## Output directory structure

After download, FASTQ files land at:
```
${DATA_ROOT}/Lung/Lung_raw_data/CNX0739{875..880}/CNR0841{131..136}/
  CL10016XXXX_LXX_read_1.fq.gz
  CL10016XXXX_LXX_read_2.fq.gz
```

## Next step: dnbc4tools alignment

This is **MGI** sequencing, NOT 10x — so no Cell Ranger. Use `dnbc4tools atac run` per sample:

```bash
# Build mm10 reference first (one-time)
sbatch 00_download/mkref_MGI.sh

# Per-sample alignment (6 jobs, run in parallel)
for i in 75 76 77 78 79 80; do
  sbatch 01_preprocessing/lung/dnbc4tools_per_sample/Run_MGI_${i}.slurm
done
```

Outputs land at `${DATA_ROOT}/Lung/Lung_cellatac/{CL_id}/outs/` — equivalent to Cell Ranger's `outs/` directory (fragments.tsv.gz, peaks.bed, filter_peak_matrix/, singlecell.csv).

Then proceed to `01_preprocessing/lung/atac_Lung.Rmd` to build the merged Seurat object.

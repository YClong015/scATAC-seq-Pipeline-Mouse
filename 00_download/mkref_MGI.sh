#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16         
#SBATCH --mem=64G                  
#SBATCH --job-name=dnbc4tools_mkref_mm10  
#SBATCH --time=2:00:00             
#SBATCH --partition=general
#SBATCH --account=a_nefzger
#SBATCH --output=DWLD_%j.out 
#SBATCH --error=DWLD_%j.err 
#SBATCH --array=1-1

echo "--- Job Started ---"
echo "Building mm10 reference for dnbc4tools..."

/scratch/user/s4869245/dnbc4tools/dnbc4tools2.1.3/dnbc4tools atac mkref \
    --fasta /scratch/user/s4869245/mouse_mm10_ref/GRCm38.primary_assembly.genome.fa \
    --ingtf /scratch/user/s4869245/mouse_mm10_ref/gencode.vM25.annotation.gtf \
    --species mm10 \
    --genomeDir /scratch/user/s4869245/mouse_mm10_ref

echo "--- Job Finished ---"

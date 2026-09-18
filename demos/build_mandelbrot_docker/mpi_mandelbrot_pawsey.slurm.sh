#!/bin/bash --login
#SBATCH --job-name=mpi-mandelbrot
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --time=00:05:00

set -euo pipefail
module load singularity/4.1.0-mpi

SINGULARITY_MPI_IMAGE="${MYSOFTWARE}/singularity/images/mandelbrot-mpi--2026.09.sif"
OUTPUT_PPM="$PWD/mandelbrot.ppm"
OUTPUT_PNG="$PWD/mandelbrot.png"

srun \
    -N "$SLURM_JOB_NUM_NODES" \
    -n "$SLURM_NTASKS" \
    -c "$SLURM_CPUS_PER_TASK" \
    singularity exec "$SINGULARITY_MPI_IMAGE" \
    mpi-mandelbrot \
        --width 1200 \
        --height 800 \
        --iterations 500 \
        --output "$OUTPUT_PPM"

singularity exec "$SINGULARITY_MPI_IMAGE" \
    convert "$OUTPUT_PPM" "$OUTPUT_PNG"

rm -f "$OUTPUT_PPM"
printf 'Created %s\n' "$OUTPUT_PNG"

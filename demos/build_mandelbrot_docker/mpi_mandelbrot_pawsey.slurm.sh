#!/bin/bash --login
#SBATCH --job-name=mpi-mandelbrot
#SBATCH --partition=work
#SBATCH --reservation=ContainersTraining
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --time=00:05:00

#--- Error handling
set -euo pipefail

#--- Load modules and define image
module load singularity/4.1.0-mpi
SINGULARITY_MPI_IMAGE="${MYSOFTWARE}/singularity/images/mandelbrot-mpi--2026.09.sif"

#--- Set variables
OUTPUT_PPM="$PWD/mandelbrot.ppm"
OUTPUT_PNG="$PWD/mandelbrot.png"

#--- Run the MPI Mandelbrot program in parallel using srun and the Singularity image
srun \
    -N "$SLURM_JOB_NUM_NODES" \
    -n "$SLURM_NTASKS" \
    -c 1 \
    singularity exec "$SINGULARITY_MPI_IMAGE" \
    mpi-mandelbrot \
        --width 1200 \
        --height 800 \
        --iterations 500 \
        --output "$OUTPUT_PPM"

#--- Convert the output PPM file to PNG format using ImageMagick's convert command inside the Singularity image
singularity exec "$SINGULARITY_MPI_IMAGE" convert "$OUTPUT_PPM" "$OUTPUT_PNG"
rm -f "$OUTPUT_PPM"
printf 'Created %s\n' "$OUTPUT_PNG"

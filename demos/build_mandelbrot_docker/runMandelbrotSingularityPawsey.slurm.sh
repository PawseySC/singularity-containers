#!/bin/bash --login
#SBATCH --job-name=mpi-mandelbrot
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --time=00:05:00

#--- Error handling
set -euo pipefail

#--- Load modules and define image
module load singularity/4.1.0-mpi
MANDEL_IMAGE="${MYSOFTWARE}/singularity/images/mandelbrot-mpi--2026.09.sif"

#--- Set environment and default values for the Mandelbrot workload and view
WIDTH="${WIDTH:-6000}"
HEIGHT="${HEIGHT:-4000}"
ITERATIONS="${ITERATIONS:-2000}"
CENTRE_REAL="${CENTRE_REAL:--0.743643887037151}"
CENTRE_IMAGINARY="${CENTRE_IMAGINARY:-0.131825904205330}"
SCALE="${SCALE:-0.002}"

#--- Set output files
OUTPUT_DIR="$PWD/output"
mkdir -p "$OUTPUT_DIR"
FILE_PPM="mandelbrot.singularity.setonix.ppm"
FILE_PNG="mandelbrot.singularity.setonix.png"

#--- Report the selected Mandelbrot view
printf 'Image size: %s x %s\n' "$WIDTH" "$HEIGHT"
printf 'Maximum iterations: %s\n' "$ITERATIONS"
printf 'Centre: (%s, %s)\n' "$CENTRE_REAL" "$CENTRE_IMAGINARY"
printf 'Scale: %s\n' "$SCALE"

#--- Run the MPI Mandelbrot program in parallel using srun and the Singularity image
srun \
    -N "$SLURM_JOB_NUM_NODES" \
    -n "$SLURM_NTASKS" \
    -c 1 \
    singularity exec "$MANDEL_IMAGE" \
    mpi-mandelbrot \
        --width "$WIDTH" \
        --height "$HEIGHT" \
        --iterations "$ITERATIONS" \
        --centre-real "$CENTRE_REAL" \
        --centre-imaginary "$CENTRE_IMAGINARY" \
        --scale "$SCALE" \
        --output "$OUTPUT_DIR/$FILE_PPM"

#--- Convert the output PPM file to PNG format using ImageMagick's convert command inside the Singularity image
singularity exec "$MANDEL_IMAGE" \
   convert "$OUTPUT_DIR/$FILE_PPM" "$OUTPUT_DIR/$FILE_PNG"

rm -f "$OUTPUT_DIR/$FILE_PPM"
printf 'Created %s\n' "$OUTPUT_DIR/$FILE_PNG"

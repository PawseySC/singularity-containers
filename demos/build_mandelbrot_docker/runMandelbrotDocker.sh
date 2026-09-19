#!/usr/bin/env bash

#--- Error handling
set -euo pipefail

#--- Set environment and default values for variables
MPI_IMAGE="${MPI_IMAGE:-mandelbrot-mpi:2026.09}"
MPI_PROCESSES="${MPI_PROCESSES:-4}"
WIDTH="${WIDTH:-1200}"
HEIGHT="${HEIGHT:-800}"
ITERATIONS="${ITERATIONS:-500}"
CENTRE_REAL="${CENTRE_REAL:--0.5}"
CENTRE_IMAGINARY="${CENTRE_IMAGINARY:-0.0}"
SCALE="${SCALE:-3.0}"

#--- Set output files
OUTPUT_DIR="${OUTPUT_DIR:-$PWD/output}"
mkdir -p "$OUTPUT_DIR"
FILE_PPM="mandelbrot.ppm"
FILE_PNG="mandelbrot.png"

#--- Report the selected Mandelbrot view
printf 'Image size: %s x %s\n' "$WIDTH" "$HEIGHT"
printf 'Maximum iterations: %s\n' "$ITERATIONS"
printf 'Centre: (%s, %s)\n' "$CENTRE_REAL" "$CENTRE_IMAGINARY"
printf 'Scale: %s\n' "$SCALE"

#--- Run the MPI Mandelbrot program in parallel using internal mpiexec
docker run \
    --rm \
    --platform linux/amd64 \
    --mount type=bind,source="$OUTPUT_DIR",target=/output \
    "$MPI_IMAGE" \
    mpiexec -n "$MPI_PROCESSES" mpi-mandelbrot \
        --width "$WIDTH" \
        --height "$HEIGHT" \
        --iterations "$ITERATIONS" \
        --centre-real "$CENTRE_REAL" \
        --centre-imaginary "$CENTRE_IMAGINARY" \
        --scale "$SCALE" \
        --output "/output/$FILE_PPM"

#--- Convert the output PPM file to PNG format using ImageMagick's convert command inside the Docker image
docker run \
    --rm \
    --platform linux/amd64 \
    --mount type=bind,source="$OUTPUT_DIR",target=/output \
    "$MPI_IMAGE" \
    convert "/output/$FILE_PPM" "/output/$FILE_PNG"

rm -f "$OUTPUT_DIR/$FILE_PPM"
printf 'Created %s\n' "$OUTPUT_DIR/$FILE_PNG"

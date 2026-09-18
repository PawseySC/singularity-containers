#!/usr/bin/env bash
set -euo pipefail

MPI_IMAGE="${MPI_IMAGE:-mandelbrot-mpi:2026.09}"
MPI_PROCESSES="${MPI_PROCESSES:-4}"
OUTPUT_DIR="${OUTPUT_DIR:-$PWD/output}"

mkdir -p "$OUTPUT_DIR"

docker run \
    --rm \
    --platform linux/amd64 \
    --mount type=bind,source="$OUTPUT_DIR",target=/output \
    "$MPI_IMAGE" \
    mpiexec -n "$MPI_PROCESSES" mpi-mandelbrot \
        --width 1200 \
        --height 800 \
        --iterations 500 \
        --output /output/mandelbrot.ppm

docker run \
    --rm \
    --platform linux/amd64 \
    --mount type=bind,source="$OUTPUT_DIR",target=/output \
    "$MPI_IMAGE" \
    convert /output/mandelbrot.ppm /output/mandelbrot.png

rm -f "$OUTPUT_DIR/mandelbrot.ppm"
printf 'Created %s\n' "$OUTPUT_DIR/mandelbrot.png"

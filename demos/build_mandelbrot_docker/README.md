# MPI Mandelbrot Docker build example

This directory contains the build context for an MPI Mandelbrot renderer built
on Pawsey's MPICH base image.

Files:

- `mandelbrot_mpi.dockerfile`: Docker/OCI build recipe
- `mpi-mandelbrot.cpp`: attributed MPI C++ application
- `render-mandelbrot`: wrapper that produces PNG output
- `THIRD_PARTY_NOTICES.md`: upstream acknowledgement and MIT licence

Build:

    docker build --platform linux/amd64 \
        --file mandelbrot_mpi.dockerfile \
        --tag mandelbrot-mpi:2026.09 .

Run from this directory:

    docker run --rm --platform linux/amd64 \
        --mount type=bind,source="$PWD",target=/work \
        mandelbrot-mpi:2026.09 \
        render-mandelbrot --processes 4 --output /work/mandelbrot.png

# Build the application on Pawsey's Setonix-compatible MPICH base image
FROM quay.io/pawsey/mpich-base:3.4.3_ubuntu24.04

# Record standard OCI image metadata
LABEL org.opencontainers.image.title="MPI Mandelbrot renderer" \
      org.opencontainers.image.description="MPI training application built on Pawsey's MPICH base image" \
      org.opencontainers.image.vendor="Pawsey Supercomputing Research Centre" \
      org.opencontainers.image.licenses="MIT"

# Install the utility used to convert the PPM result to PNG
RUN set -eux; \
    export DEBIAN_FRONTEND=noninteractive; \
    apt-get update; \
    apt-get install -y --no-install-recommends imagemagick; \
    apt-get clean; \
    rm -rf /var/lib/apt/lists/*

# Copy and compile the MPI application with the compiler from the base image
COPY mpi-mandelbrot.cpp /tmp/mpi-mandelbrot.cpp
RUN mpic++ \
        -std=c++17 \
        -O3 \
        -Wall \
        -Wextra \
        -Wpedantic \
        -o /usr/local/bin/mpi-mandelbrot \
        /tmp/mpi-mandelbrot.cpp \
    && rm -f /tmp/mpi-mandelbrot.cpp

# Install the wrapper that launches MPI and converts the result to PNG
COPY render-mandelbrot /usr/local/bin/render-mandelbrot
RUN chmod 0755 /usr/local/bin/render-mandelbrot

# Preserve third-party acknowledgements and licence information
COPY THIRD_PARTY_NOTICES.md \
    /usr/local/share/doc/mpi-mandelbrot/THIRD_PARTY_NOTICES.md

# Display the wrapper help when no other command is supplied
CMD ["render-mandelbrot", "--help"]

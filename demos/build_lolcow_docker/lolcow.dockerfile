# Start from a versioned Ubuntu image on Docker Hub
FROM docker.io/ubuntu:24.04

# Install the applications and remove package-manager cache files
RUN set -eux; \
    export DEBIAN_FRONTEND=noninteractive; \
    apt-get update; \
    apt-get install -y --no-install-recommends \
        cowsay \
        lolcat; \
    apt-get clean; \
    rm -rf /var/lib/apt/lists/*

# Make the installed commands available by name at runtime
ENV PATH="/usr/games:${PATH}"

# Copy the message displayed by the default container action
COPY lolcow-message.txt /usr/local/share/lolcow/message.txt

# Define the default action for docker run and singularity run
CMD ["bash", "-c", "cowsay < /usr/local/share/lolcow/message.txt | lolcat"]

# Preserve the recipe and build input files inside the image
ARG IMAGE_BUILD_INFO_DIR="/opt/build-info-and-recipes/lolcow"
RUN mkdir -p "${IMAGE_BUILD_INFO_DIR}"
COPY lolcow.dockerfile \
     lolcow-message.txt \
     "${IMAGE_BUILD_INFO_DIR}/"

# Record standard OCI image metadata and the build-information location
LABEL org.opencontainers.image.title="lolcow training image" \
      org.opencontainers.image.description="Small image used to teach Docker builds for HPC" \
      org.opencontainers.image.vendor="Pawsey Supercomputing Research Centre" \
      org.opencontainers.image.source="https://github.com/PawseySC/singularity-containers" \
      au.org.pawsey.image.build-info-dir="${IMAGE_BUILD_INFO_DIR}"

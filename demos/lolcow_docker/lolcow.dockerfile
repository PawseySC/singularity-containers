FROM docker.io/ubuntu:24.04

LABEL org.opencontainers.image.title="lolcow training image" \
      org.opencontainers.image.description="Small image used to teach Docker builds for HPC" \
      org.opencontainers.image.vendor="Pawsey Supercomputing Research Centre"

RUN set -eux; \
    export DEBIAN_FRONTEND=noninteractive; \
    apt-get update; \
    apt-get install -y --no-install-recommends \
        cowsay \
        fortune-mod \
        lolcat; \
    apt-get clean; \
    rm -rf /var/lib/apt/lists/*

ENV PATH="/usr/games:${PATH}"

RUN useradd --create-home --uid 1000 training
USER training
WORKDIR /home/training

CMD ["bash", "-c", "fortune | cowsay | lolcat"]

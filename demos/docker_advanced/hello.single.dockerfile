FROM docker.io/ubuntu:24.04

RUN set -eux; \
    export DEBIAN_FRONTEND=noninteractive; \
    apt-get update; \
    apt-get install -y --no-install-recommends g++; \
    apt-get clean; \
    rm -rf /var/lib/apt/lists/*

COPY hello.cpp /tmp/hello.cpp
RUN g++ -O2 -Wall -Wextra -Wpedantic \
        -o /usr/local/bin/hello.exe \
        /tmp/hello.cpp

ENTRYPOINT ["hello.exe"]
CMD ["Hello from the single-stage image"]

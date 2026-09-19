FROM docker.io/ubuntu:24.04 AS build

RUN set -eux; \
    export DEBIAN_FRONTEND=noninteractive; \
    apt-get update; \
    apt-get install -y --no-install-recommends g++; \
    rm -rf /var/lib/apt/lists/*

COPY hello.cpp /tmp/hello.cpp
RUN mkdir -p /out \
    && g++ -O2 -Wall -Wextra -Wpedantic \
        -o /out/hello.exe \
        /tmp/hello.cpp

FROM docker.io/ubuntu:24.04

RUN set -eux; \
    export DEBIAN_FRONTEND=noninteractive; \
    apt-get update; \
    apt-get install -y --no-install-recommends file; \
    rm -rf /var/lib/apt/lists/*

COPY --from=build /out/hello.exe /usr/local/bin/hello.exe

ENTRYPOINT ["hello.exe"]
CMD ["Hello from the multi-stage image"]

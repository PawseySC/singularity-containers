# Build stage: install the compiler and compile the application
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

# Runtime stage: start again from Ubuntu without the build tools
FROM docker.io/ubuntu:24.04

# Install only the utility required for the later inspection exercise
RUN set -eux; \
    export DEBIAN_FRONTEND=noninteractive; \
    apt-get update; \
    apt-get install -y --no-install-recommends file; \
    rm -rf /var/lib/apt/lists/*

# Copy only the compiled executable from the build stage
COPY --from=build /out/hello.exe /usr/local/bin/hello.exe

# Run the application with a default message
ENTRYPOINT ["hello.exe"]
CMD ["Hello from the multi-stage image"]

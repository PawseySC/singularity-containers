# Use Ubuntu for both compilation and runtime
FROM docker.io/ubuntu:24.04

# Install the compiler and remove temporary package data in the same layer
RUN set -eux; \
    export DEBIAN_FRONTEND=noninteractive; \
    apt-get update; \
    apt-get install -y --no-install-recommends g++; \
    apt-get clean; \
    rm -rf /var/lib/apt/lists/*

# Copy and compile the application inside the final image
COPY hello.cpp /tmp/hello.cpp
RUN g++ -O2 -Wall -Wextra -Wpedantic \
        -o /usr/local/bin/hello.exe \
        /tmp/hello.cpp

# Run the application with a default message
ENTRYPOINT ["hello.exe"]
CMD ["Hello from the single-stage image"]

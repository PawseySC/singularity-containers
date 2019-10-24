FROM ubuntu:18.04

MAINTAINER Your Name <youremail@yourdomain>

RUN apt-get update && \
    apt-get install -y \
        build-essential \
        git \
        wget \
    && apt-get clean all \
    && rm -rf /var/lib/apt/lists/*

CMD ["/bin/bash"]

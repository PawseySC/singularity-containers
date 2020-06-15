FROM ubuntu:18.04

LABEL maintainer="Pawsey Supercomputing Centre"
LABEL version="v0.0.1"

RUN apt-get -y update && \
  apt-get -y install fortune cowsay lolcat

ENV PATH=/usr/games:$PATH

CMD fortune | cowsay | lolcat

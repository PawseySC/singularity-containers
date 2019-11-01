FROM ubuntu:18.04

MAINTAINER Pawsey Supercomputing Centre

RUN apt-get -y update && \
  apt-get -y install fortune cowsay lolcat

ENV PATH=/usr/games:$PATH

VOLUME /data
WORKDIR /data

CMD ["/bin/bash"]

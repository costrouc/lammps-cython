FROM costrouc/lammps:stable_16Mar2018-debian-mpi-none
MAINTAINER Chris Ostrouchov

ARG VERSION=v0.2
ARG USERNAME=costrouc
ARG PROJECT=lammps-cython

RUN apt update && \
    apt install python3.5 python3-pip -y && \
    rm -rf /var/lib/apt/lists/*

RUN pip3 install cython numpy mpi4py

RUN pip3 install --no-cache-dir https://gitlab.com/$USERNAME/$PROJECT/repository/$VERSION/archive.tar.gz

FROM firedrakeproject/firedrake

LABEL maintainer="Nguyen Chien <nchien@ed.ac.uk>" 
LABEL version='main'

# Making use of https://hub.docker.com/r/firedrakeproject/firedrake/dockerfile

USER root

RUN mkdir /model_data
RUN mkdir /model
WORKDIR /model
RUN mkdir inputs
RUN mkdir tools
RUN mkdir share

COPY ./*  ./
COPY ./inputs/* ./inputs/
COPY ./tools/* ./tools/

# Set up user `thetis`, avoid runing as root
RUN useradd -m -s /bin/bash -G sudo thetis && \
    echo "thetis:docker" | chpasswd && \
    echo "thetis ALL=(ALL) NOPASSWD: ALL" >> /etc/sudoers && \
    ldconfig

ENV OMP_NUM_THREADS 1
ENV OPENBLAS_NUM_THREADS 1

CMD . /home/firedrake/firedrake/bin/activate &&  \
    pip3 install -r requirements.txt && \
    apt-get update && apt-get -y install mpich && \
    sh sim.sh && \
    /bin/bash

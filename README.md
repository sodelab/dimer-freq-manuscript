# Builds a Docker image for
# https://github.com/sodelab/dimer-freq-manuscript.git
# docker build . -t <MYIMAGE>
# docker run -it <MYIMAGE> bash

FROM continuumio/miniconda:4.4.10
LABEL maintainer "Murat Keceli <keceli@gmail.com>"

RUN apt-get update && \
   apt-get install -y --no-install-recommends build-essential \
                    cmake \
                    wget \
                    gcc \
                    gfortran \
                    libopenblas-dev \
                    liblapack-dev \
                    libgsl2 \
                    vim \
                    libgsl0-dev && \
   apt-get clean


RUN mkdir /container && \
    cd /container && \
    git clone https://github.com/sodelab/dimer-freq-manuscript.git 

RUN cd /container && \
    cd dimer-freq-manuscript/tools && \
    wget --content-disposition 'https://sourceforge.net/projects/nitrogen-downloads/files/nitrogen_v1.9.tar.gz/download' && \
    tar -xvzf nitrogen_v1.9.tar.gz && \
    rm -rf nitrogen_v1.9.tar.gz && \
    cd nitro* && \
    make 

RUN cd /container && \
    cd dimer-freq-manuscript/tools && \
    git clone https://github.com/keceli/MaVi.git && \
    cd MaVi/Debug && \
    make

#RUN cd /container && \ 
#    cd dimer-freq-manuscript/tools && \
#    git clone https://github.com/sodelab/sindo-Yagi.git && \
#    cd sindo
#    cd MaVi/Debug && \
#    make

RUN cd /container && \ 
    cd dimer-freq-manuscript/tools/library && \
    make && \
    cd ../optimize && \
    make && \
    cd ../anharmonic && \
    make && \
    cd ../evaluate && \
    make

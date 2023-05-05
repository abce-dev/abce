# This Dockerfile builds, installs, and tests ABCE in a clean environment
FROM ubuntu:20.04

SHELL ["/bin/bash", "-c"]

WORKDIR /home

## Commands that need to be run from root
# Install packages where apt requires root
RUN apt update && \
    apt -y upgrade && \
    apt install -y \
        git \
        wget \
        python3 \
        python3-distutils \
        python3-dev \
        python3-pip

RUN pip install python-config

# Install tzdata
RUN DEBIAN_FRONTEND=noninteractive TZ=America/Chicago apt -y install tzdata

# Install julia
RUN wget https://julialang-s3.julialang.org/bin/linux/x64/1.8/julia-1.8.2-linux-x86_64.tar.gz && \
    tar zxvf julia-1.8.2-linux-x86_64.tar.gz

# Set up environment variables
ENV PATH="$PATH:/home/julia-1.8.2/bin"

# Set up a directory for the ABCE environment files
RUN mkdir /home/abce_env
WORKDIR /home/abce_env

# Bring the python pip requirements file into the docker image scope
COPY ./requirements.txt .

# Install required python packages
RUN pip install -r ./requirements.txt

# Bring the julia environment setup script into the docker image scope
COPY ./make_julia_environment.jl .

# Run the julia installation script to set up the julia package environment
RUN julia make_julia_environment.jl

# Set up environment variables
ENV ABCE_DIR="/home/abce"
ENV ALEAF_DIR=""
ENV ABCE_ENV="/home/abce_env"

WORKDIR /home

LABEL org.opencontainers.image.source=https://github.com/abce-dev/abce

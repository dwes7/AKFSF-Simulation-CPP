FROM ubuntu:focal-20220826

ENV DEBIAN_FRONTEND=noninteractive
WORKDIR /workspace

RUN apt update -y && apt upgrade -y && apt install  -y \
    tzdata \
    cmake \
    git \ 
    vim \
    libeigen3-dev \
    libsdl2-dev \ 
    libsdl2-ttf-dev \
    gdb \ 
    build-essential


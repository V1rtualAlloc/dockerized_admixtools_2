FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y \
        build-essential \
        gfortran \
        git \
        libgsl-dev \
        liblapack-dev \
        liblapacke-dev \
        libopenblas-dev && \
    rm -rf /var/lib/apt/lists/*

RUN git clone --depth 1 https://github.com/DReichLab/AdmixTools /opt/AdmixTools && \
    cd /opt/AdmixTools/src && \
    make clobber && \
    make LDLIBS="-llapacke" install && \
    find /opt/AdmixTools/bin -maxdepth 1 -type f -exec cp {} /usr/local/bin/ \;

RUN apt-get update && \
    apt-get install -y wget unzip && \
    rm -rf /var/lib/apt/lists/* && \
    wget -q https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231211.zip \
        -O /tmp/plink.zip && \
    unzip /tmp/plink.zip plink -d /usr/local/bin/ && \
    rm /tmp/plink.zip

WORKDIR /data

ENTRYPOINT ["convertf"]

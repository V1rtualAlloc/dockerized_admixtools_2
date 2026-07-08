FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y \
        build-essential \
        git \
        libgsl-dev \
        liblapack-dev \
        libopenblas-dev \
        libfftw3-dev && \
    rm -rf /var/lib/apt/lists/*

RUN git clone --depth 1 https://github.com/priyamoorjani/DATES /opt/DATES

WORKDIR /opt/DATES

# The repo ships a pre-built binary linked against libgsl.so.19.
# Ubuntu 22.04 provides libgsl.so.27 — create a compat symlink.
RUN GSL_SO=$(find /usr/lib -name 'libgsl.so.*' -not -name '*.la' | sort -V | tail -1) && \
    echo "System GSL: $GSL_SO" && \
    ln -sf "$GSL_SO" /usr/lib/x86_64-linux-gnu/libgsl.so.19 && \
    ldconfig

# Find and install the dates binary; print all deps so missing libs surface at build time
RUN find /opt/DATES -maxdepth 3 -type f -executable \
        \( -name 'dates' -o -name 'DATES' -o -name 'dates_main' \) \
    | head -1 \
    | xargs -I{} sh -c 'echo "Using binary: $1" && ldd "$1" && cp "$1" /usr/local/bin/dates' _ {}

WORKDIR /data
ENTRYPOINT ["dates"]

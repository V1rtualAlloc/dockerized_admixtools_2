FROM rocker/tidyverse:4.4

RUN apt-get update && \
    apt-get install -y \
        libgsl-dev \
        libopenblas-dev && \
    rm -rf /var/lib/apt/lists/*

RUN Rscript -e " \
    install.packages(c('remotes', 'igraph', 'foreach', 'doParallel', 'data.table'), \
                     repos = 'https://cloud.r-project.org'); \
    remotes::install_github('uqrmaie1/admixtools', upgrade = 'never')"

WORKDIR /data

ENTRYPOINT ["Rscript"]

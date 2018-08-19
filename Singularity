BootStrap: docker
From: continuumio/miniconda:latest
IncludeCmd: yes

%post
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# this will install all necessary packages and prepare the container
    apt-get -y update
    apt-get -y install make gcc zlib1g-dev libncurses5-dev
    wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 \
        && tar -xjf samtools-1.9.tar.bz2 \
        && cd samtools-1.9 \
        && make \
        && make prefix=/usr/local install
    export PATH=/opt/conda/bin:$PATH
    conda install --yes -c intel \
        intelpython2_core \
        cython \
        scipy
    conda install --yes -c bioconda \
        pysam \
        matplotlib
    conda clean --index-cache --tarballs --packages --yes
    mkdir -p /athena /scratchLocal


%runscript
#!/bin/bash
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# this text will will run whenever the container is called as an executable
function usage() {
    cat <<EOF

NAME
    atacflow - atac-seq dependencies signularity container 

SYNOPSIS
    rnaseq tool [tool options]
    rnaseq list
    rnaseq help

DESCRIPTION
    Singularity container with nucleoatac. 

EOF
}

function tools() {
    echo "conda: $(which conda)"
    echo "---------------------------------------------------------------"
    conda list
    echo "---------------------------------------------------------------"
    echo "samtools: $(samtools --version | head -n1)"
    echo "nucleoatac: $( --version | head -n1)"
}

arg="${1:-none}"

case "$arg" in
    none) usage; exit 1;;
    help) usage; exit 0;;
    list) tools; exit 0;;
    # just try to execute it then
    *)    $@;;
esac

%environment
export PATH="/opt/conda/bin:/usr/local/bin:/usr/bin:/bin:"
unset CONDA_DEFAULT_ENV
export ANACONDA_HOME=/opt/conda


FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="bc19c784b0bfbaeee5e888b2adc543a25aebe0a403b47afcfffe58e02781a355"

# now I add all my non conda stuff

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y \
git build-essential zlib1g-dev libbz2-dev \
liblzma-dev samtools bedtools subread wget \
python3 pip unzip openjdk-11-jre-headless

WORKDIR /data
RUN wget --output-document sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
RUN tar -vxzf sratoolkit.tar.gz
ENV PATH="$PATH:/data/sratoolkit.3.0.1-ubuntu64/bin"


RUN git clone --recursive https://github.com/mhammell-laboratory/TElocal.git
WORKDIR TElocal
RUN pip install setuptools
RUN python3 setup.py install


WORKDIR /data
RUN wget https://github.com/alexdobin/STAR/archive/refs/tags/2.7.1a.tar.gz && tar -xvzf 2.7.1a.tar.gz
RUN cd STAR-2.7.1a/source; make STAR
ENV PATH="$PATH:/data/STAR-2.7.1a/source"

# ENV PATH="/root/miniconda3/bin:${PATH}"
# ARG PATH="/root/miniconda3/bin:${PATH}"
# RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
#     && mkdir /root/.conda \
#     && bash Miniconda3-latest-Linux-x86_64.sh -b \
#     && rm -f Miniconda3-latest-Linux-x86_64.sh \
#     && echo "Running $(conda --version)" && \
#     conda init bash && \
#     . /root/.bashrc

RUN wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip \
    && unzip Trimmomatic-0.39.zip


# Now I add all my environments

# Step 1: Retrieve conda environments

# Conda environment:
#   source: envs/deeptools.yml
#   prefix: /conda-envs/2b6fafbb07fcb29c8cfbb32b1fb919f1
#   name: deeptools
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#     - r
#   dependencies:
#     - bedtools
#     - pygenometracks
#     - deeptools
#     - ipykernel
#     - ipython
#     - python=3.9
#   prefix: /home/mk/anaconda3/envs/deeptools
RUN mkdir -p /conda-envs/2b6fafbb07fcb29c8cfbb32b1fb919f1
COPY envs/deeptools.yml /conda-envs/2b6fafbb07fcb29c8cfbb32b1fb919f1/environment.yaml

# Conda environment:
#   source: envs/renv.yml
#   prefix: /conda-envs/59c06ea203069ce2f424aefba5e2efc5
#   name: renv
#   channels:
#    - conda-forge
#    - bioconda
#    - defaults
#   dependencies:
#    - r-base=4.1.1
#    - r-tibble
#    - r-readr
#    - r-dplyr
#    - r-ggplot2
#    - r-pheatmap
#    - r-cowplot
#    - bioconductor-DESeq2
#    - bioconductor-GenomicFeatures
RUN mkdir -p /conda-envs/59c06ea203069ce2f424aefba5e2efc5
COPY envs/renv.yml /conda-envs/59c06ea203069ce2f424aefba5e2efc5/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/2b6fafbb07fcb29c8cfbb32b1fb919f1 --file /conda-envs/2b6fafbb07fcb29c8cfbb32b1fb919f1/environment.yaml && \
    mamba env create --prefix /conda-envs/59c06ea203069ce2f424aefba5e2efc5 --file /conda-envs/59c06ea203069ce2f424aefba5e2efc5/environment.yaml && \
    mamba clean --all -y



ENTRYPOINT [ "/bin/bash" ]


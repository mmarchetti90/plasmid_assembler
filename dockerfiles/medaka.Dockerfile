
### SOFTWARE VERSIONS ---------------------- ###

ARG CUDA_IMAGE=nvidia/cuda:12.5.1-cudnn-devel-ubuntu20.04
ARG CONDA_IMAGE=continuumio/miniconda3:25.3.1-1
ARG DV_GPU_BUILD=1
ARG HTSLIB=1.22
ARG MEDAKA=2.1.0
ARG MINIMAP2=2.29
ARG SAMTOOLS=1.22

### GET CONDA PACKAGES --------------------- ###

FROM ${CONDA_IMAGE} as conda_img

# Update conda
RUN conda update -y conda

# Add channels
RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge

# Install mamba
RUN conda install -y \
    mamba \
    python=3.11

# Installing software
RUN mamba install -y \
    htslib=${HTSLIB} \
    medaka=${MEDAKA} \
    minimap2=${MINIMAP2} \
    samtools=${SAMTOOLS} && \
    conda clean -afty

### COPY CONDA TO CUDA IMAGE --------------- ###

FROM ${CUDA_IMAGE} as cuda_img

COPY --from=conda_img /opt/conda /opt/conda

ENV DV_GPU_BUILD=${DV_GPU_BUILD}

ENV PATH $PATH:/opt/conda/bin

### SETTING WORKING ENVIRONMENT ------------ ###

# Set workdir to /home/
WORKDIR /home/

# Launch bash automatically
CMD ["/bin/bash"]

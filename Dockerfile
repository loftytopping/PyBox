# ================================================================================
#
# This Dockerfile sets up a container based on the Ubuntu 16:04 image to run PyBox
# 
FROM ubuntu:16.04
MAINTAINER David Topping (david.topping@manchester.ac.uk)

# ========================== Basic configuration =================================
RUN apt-get update 
RUN apt-get install -y build-essential \
    apt-utils \
    wget \
    libgl1-mesa-glx \
    gfortran \
    linux-headers-generic \
    cmake \
    vim \
    unzip \
    pkgconf \
    libpng-dev \
    libfreetype6-dev \
    libfontconfig1 \
    libxrender1 \
    xauth \
    git \
    subversion 
# Install BLAS and LAPACK
RUN apt-get install -y \
    libblas-dev \
    liblapack-dev \
    libblas-doc \
    liblapack-doc 

# ========================== Create directory structure ==========================
RUN mkdir -p /Code
RUN mkdir -p /Code/Anaconda
RUN mkdir -p /Code/Git_repos
RUN mkdir -p /Downloads

# ========================== Install Anaconda Python =============================
WORKDIR /Code/Anaconda
RUN wget -q https://repo.continuum.io/archive/Anaconda3-5.1.0-Linux-x86_64.sh
RUN printf '\nyes\n\nyes\nno\n' | bash Anaconda3-5.1.0-Linux-x86_64.sh 
RUN rm Anaconda3-5.1.0-Linux-x86_64.sh 
# - add anaconda python to path
ENV PATH="/root/anaconda3/bin/:${PATH}"
# - add relevant conda channels and install modules for UManSysProp and PyBox
RUN conda config --append channels conda-forge
RUN conda install -c openbabel openbabel 
RUN conda install flask-wtf
RUN conda install -c chria assimulo

# =========================== Clone UManSysProp ==================================
WORKDIR /Code/Git_repos
RUN git clone https://github.com/loftytopping/UManSysProp_public.git
# ============================= Clone PyBox ======================================
RUN git clone https://github.com/loftytopping/PyBox.git

WORKDIR /Code



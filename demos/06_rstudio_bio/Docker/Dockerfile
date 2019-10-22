FROM rocker/tidyverse:3.5

RUN apt-get update -qq && apt-get -y --no-install-recommends install \
      autoconf \
      automake \
      g++ \
      gcc \
      gfortran \
      make \
      && apt-get clean all \
      && rm -rf /var/lib/apt/lists/*

RUN mkdir -p $HOME/.R
COPY Makevars /root/.R/Makevars

RUN Rscript -e "library('devtools')" \
      -e "install_github('Rdatatable/data.table', build_vignettes=FALSE)" \
      -e "install.packages('reshape2')" \
      -e "install.packages('fields')" \
      -e "install.packages('ggbeeswarm')" \
      -e "install.packages('gridExtra')" \
      -e "install.packages('dynamicTreeCut')" \
      -e "install.packages('DEoptimR')" \
      -e "install.packages('http://cran.r-project.org/src/contrib/Archive/robustbase/robustbase_0.90-2.tar.gz', repos=NULL, type='source')" \
      -e "install.packages('dendextend')" \
      -e "install.packages('RColorBrewer')" \
      -e "install.packages('locfit')" \
      -e "install.packages('KernSmooth')" \
      -e "install.packages('BiocManager')" \
      -e "source('http://bioconductor.org/biocLite.R')" \
      -e "biocLite('Biobase')" \
      -e "biocLite('BioGenerics')" \
      -e "biocLite('BiocParallel')" \
      -e "biocLite('SingleCellExperiment')" \
      -e "biocLite('GenomeInfoDb')" \
      -e "biocLite('GenomeInfgoDbData')" \
      -e "biocLite('DESeq')" \
      -e "biocLite('DESeq2')" \
      -e "BiocManager::install(c('scater', 'scran'))" \
      -e "library('devtools')" \
      -e "install_github('IMB-Computational-Genomics-Lab/ascend', ref = 'devel')" \
      && rm -rf /tmp/downloaded_packages

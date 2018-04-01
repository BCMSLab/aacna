FROM bioconductor/release_base2:latest

RUN apt-get update
RUN apt-get -y upgrade

RUN apt-get -y install libpng-dev
RUN apt-get -y install libpcre3-dev
RUN apt-get -y install liblzma-dev
RUN apt-get -y install libbz2-dev
RUN apt-get -y install libudunits2-dev

RUN Rscript -e "update.packages(ask=FALSE)"
RUN Rscript -e "BiocInstaller::biocLite('Biobase', quiet = TRUE)"
RUN Rscript -e "BiocInstaller::biocLite('AnnotationDbi', quiet = TRUE)"
RUN Rscript -e "BiocInstaller::biocLite('GO.db', quiet = TRUE)"
RUN Rscript -e "BiocInstaller::biocLite('org.Mm.eg.db', quiet = TRUE)"
RUN Rscript -e "BiocInstaller::biocLite('GEOquery', quiet = TRUE)"
RUN Rscript -e "BiocInstaller::biocLite('limma', quiet = TRUE)"
RUN Rscript -e "BiocInstaller::biocLite('clusterProfiler', quiet = TRUE)"
RUN Rscript -e "BiocInstaller::biocLite('STRINGdb', quiet = TRUE)"
RUN Rscript -e "BiocInstaller::biocLite('WGCNA', quiet = TRUE)"
RUN Rscript -e "BiocInstaller::biocLite('mogene10sttranscriptcluster.db', quiet = TRUE)"

RUN Rscript -e "install.packages('devtools', quiet = TRUE)"
RUN Rscript -e "install.packages('RSQLite', quiet = TRUE)"
RUN Rscript -e "install.packages('tidyverse', quiet = TRUE)"
RUN Rscript -e "install.packages('reshape2', quiet = TRUE)"
RUN Rscript -e "install.packages('igraph', quiet = TRUE)"
RUN Rscript -e "install.packages('xtable', quiet = TRUE)"
RUN Rscript -e "install.packages('cowplot', quiet = TRUE)"
RUN Rscript -e "install.packages('tidygraph', quiet = TRUE)"
RUN Rscript -e "install.packages('ggraph', quiet = TRUE)"
RUN Rscript -e "install.packages('covr', quiet = TRUE)"
RUN Rscript -e "install.packages('ggdendro', quiet = TRUE)"
RUN Rscript -e "install.packages('testthat', quiet = TRUE)"
RUN Rscript -e "install.packages('dynamicTreeCut', quiet = TRUE)"

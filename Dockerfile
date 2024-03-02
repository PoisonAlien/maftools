FROM ubuntu:xenial

RUN apt-get update && apt-get install -y \
	libcurl4-openssl-dev \
	libxml2-dev \
	apt-transport-https \
	libssl-dev \
	r-base-dev
	

RUN Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("ComplexHeatmap") ; biocLite("VariantAnnotation"); biocLite("Biostrings")'
#RUN Rscript -e 'install.packages("kinitr", repos="http://cran.uk.r-project.org")'
RUN Rscript -e 'install.packages("devtools", repos="http://cran.uk.r-project.org", dependencies=TRUE)'

ADD . /opt/maftools/
WORKDIR /opt/maftools

# We need to build the script locally
RUN Rscript -e 'library("devtools"); devtools::install()'

FROM ubuntu:latest
RUN apt-get update -qq -y --fix-missing
RUN apt-get install -y apt-utils
RUN apt-get install -y build-essential git autoconf wget lbzip2
RUN apt-get install -y zlib1g-dev libbz2-dev gsl-bin libgsl0-dev
RUN apt-get install -y libncurses5-dev liblzma-dev libssl-dev libcurl4-openssl-dev
RUN cd /usr/local/ \
    && wget http://research-pub.gene.com/gmap/src/gmap-gsnap-2020-04-08.tar.gz \
    && tar xvf gmap-gsnap-2020-04-08.tar.gz && cd gmap-2020-04-08 && ./configure && make && make install
RUN wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.masked.gz && gzip -d hg38.fa.masked 
    && gmap_build -d hg38 hg38.fa.masked
RUN wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.masked.gz && gzip -d hg19.fa.masked 
    && gmap_build -d hg19 hg19.fa.masked
RUN wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr19.fa.gz && gzip -d chr19.fa.gz \
    && gmap_build -d hg19_chr19 chr19.fa 
    # chr19 only for testing. 
CMD ["/bin/bash"] 

# docker image for gemBS binaries (gem3-mapper and bs_call)
FROM ubuntu:latest
RUN apt-get update -qq -y --fix-missing
RUN apt-get install -y apt-utils
RUN apt-get install -y build-essential git autoconf wget lbzip2
RUN apt-get install -y zlib1g-dev libbz2-dev gsl-bin libgsl0-dev
RUN apt-get install -y libncurses5-dev liblzma-dev libssl-dev libcurl4-openssl-dev
RUN cd /usr/local/ \
    && git clone --recursive https://github.com/smarco/gem3-mapper.git gem3-mapper \
    && cd gem3-mapper && ./configure --disable-cuda \
    && make 
RUN ln -s /usr/local/gem3-mapper/bin/* /usr/bin
RUN cd /usr/local/ \
    && wget https://github.com/samtools/htslib/releases/download/1.10.2/htslib-1.10.2.tar.bz2 \
    && tar jxvf htslib-1.10.2.tar.bz2 \
    && cd htslib-1.10.2 && ./configure && make && make install
RUN cd /usr/local/ \
    && git clone --recursive https://github.com/heathsc/bs_call bs_call \
    && cd bs_call && ./configure && make all \
    && cp bin/bs_call /usr/bin 
CMD ["/bin/bash"] 

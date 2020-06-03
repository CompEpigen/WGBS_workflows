FROM ubuntu:latest
RUN apt-get update 
RUN apt-get install -y default-jre wget 
RUN cd /usr/local/ \
    && wget https://sourceforge.net/projects/bissnp/files/BisSNP-0.71/BisSNP-0.71.jar
CMD ["/bin/bash"] 

#########################################################
# Docker file for PHOTOMOD
# Base on Ubuntu 16.04
#########################################################

FROM ubuntu:16.04

#Maintainer/Author
MAINTAINER Apiwat Sangphukieo

ADD / /working_directory

RUN apt-get update \
    && apt-get install -y wget \
    && apt-get install -y software-properties-common \
    && apt-get update \
    && apt-get install -y openjdk-8-jdk \
    && apt-get install -y r-base \
    && apt-get install -y ncbi-blast+-legacy

RUN apt-get install -y python-pip \
    && python -m pip install pandas \
    && python -m pip install scipy==0.16 \
    && python -m pip install rpy2==2.7.7 \
    && python -m pip install Biopython

RUN tar xvjf /working_directory/call_neighbor_3evalue/call_neighbor.tar.bz2 -C /working_directory/call_neighbor_3evalue/

WORKDIR /working_directory

#sudo docker build -t photomod .
#sudo docker run -ti -v `pwd`/INPUT:/working_directory/INPUT -v `pwd`/OUTPUT:/working_directory/OUTPUT -d --name photomod_con1 photomod
#sudo docker exec -t -i photomod_con1 /bin/bash -c "python run.py test_protein.fasta"

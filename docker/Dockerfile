# Standard ops
FROM ubuntu:trusty
MAINTAINER Rail Suleymanov
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update

# Install applications and python tools
RUN apt-get install -yq git wget supervisor
RUN apt-get install -yq python python-dev python-distribute python-pip
RUN wget https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_2_x/ViennaRNA-2.2.9.tar.gz
RUN tar -xvzf ViennaRNA-2.2.9.tar.gz
RUN cd ViennaRNA-2.2.9/ && ./configure --without-perl && make && make install
RUN pip install supervisor-stdout
RUN pip install gunicorn

# SSH mess and sources download


# Prepare services
RUN service supervisor stop
ADD . /opt/peptideprops
ADD supervisord.conf /etc/supervisor/conf.d/
ADD logging.conf /opt/peptideprops/

# Install requirements
RUN cd /opt/peptideprops/ && pip install -r requirements.txt

# Set the default command
ENTRYPOINT supervisord & bash


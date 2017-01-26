FROM continuumio/miniconda

# Adds the current working directory to a root ParmEd directory
COPY . /ParmEd
ARG PYTHON_VERSION=local
ENV PYTHON_VERSION ${PYTHON_VERSION}
RUN /bin/bash -c "apt-get install -y gcc g++ gfortran gromacs"
RUN /bin/bash -c "cd /ParmEd && devtools/ci/jenkins/install.sh && devtools/ci/jenkins/runtest.sh"
